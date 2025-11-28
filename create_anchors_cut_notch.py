import os
import glob
import random
import traceback
import math
import sys
import json


# ==============================================================================
# Parameter passing
# ==============================================================================
if len(sys.argv) < 2:
    raise ValueError("Usage: python create_anchors_cut_notch.py <config_filename.json>")

config_filename = sys.argv[1]
if not config_filename.endswith('.json'):
    config_filename += '.json'

with open(config_filename, 'r') as f:
    params = json.load(f)

domain_min_coords = params['domain_min_coords']
domain_max_coords = params['domain_max_coords']
circular_burrito_radius = params.get('circular_burrito_radius')
pbc_direction = params.get('pbc_direction')
bridge_thickness = params['bridge_thickness']
default_mesh_size = params['default_mesh_size']

mesh_md_region = params['mesh_md_region']
mesh_cutout_region = params['mesh_cutout_region']
have_iface = params['have_iface']
have_padding = params['have_padding']
tol = params['tol']

dpd_thickness = params['dpd_thickness']
anchor_placement_probability = params['anchor_placement_probability']
ANCHOR_ATOM_MASS_VAL = params['ANCHOR_ATOM_MASS_VAL']
truncate_chains = params['truncate_chains']
remove_bridge_interactions = params['remove_bridge_interactions']

md_remove_shape = params['md_remove_shape']

cut_thickness = params['cut_thickness']
cut_shift = params['cut_shift']

notch_thickness = params['notch_thickness']
notch_normal = params['notch_normal']
notch_plane = params['notch_plane']
notch_shift = params['notch_shift']
notch_depth = params['notch_depth']

lammps_file = params['lammps_file']
output_file = params['output_file']

# ==============================================================================
# find_lammps_data_file
# ==============================================================================
def find_lammps_data_file(script_dir, specific_name=None):
    file_path = None
    if specific_name:
        potential_path = os.path.join(script_dir, specific_name)
        if os.path.exists(potential_path): file_path = potential_path
        else: raise IOError("Specified LAMMPS file not found: {}".format(potential_path))
    else:
        data_files = sorted(glob.glob(os.path.join(script_dir, '*.data')))
        if not data_files: raise IOError("No '.data' file found in directory: {}".format(script_dir))
        if len(data_files) > 1: print("WARNING: Multiple '.data' files found.")
        file_path = data_files[0]
    base_name = os.path.splitext(os.path.basename(file_path))[0].replace(' ', '_').replace('-', '_')
    if not base_name: base_name = "DefaultModel"
    if base_name and base_name[0].isdigit(): base_name = "_" + base_name
    print("Using LAMMPS file: {}".format(file_path))
    return file_path, base_name

# ==============================================================================
# Atom style and line parsing
# ==============================================================================
def _determine_atom_style_info(atom_section_header_line):
    style_comment_parts = atom_section_header_line.split('#', 1)
    style_comment = style_comment_parts[1].strip().lower() if len(style_comment_parts) > 1 else "atomic"
    info = {'style_name': style_comment, 'id_idx': 0, 'first_extra_col_idx': 0}
    if 'molecular' in info['style_name'] or 'full' in info['style_name']:
        info['mol_id_idx'] = 1; info['type_idx'] = 2; current_idx = 3
    else:
        info['mol_id_idx'] = None; info['type_idx'] = 1; current_idx = 2
    if 'charge' in info['style_name'] or 'full' in info['style_name']:
        info['q_idx'] = current_idx; current_idx += 1
    else:
        info['q_idx'] = None
    info['x_idx'] = current_idx; info['y_idx'] = current_idx + 1; info['z_idx'] = current_idx + 2
    info['first_extra_col_idx'] = info['z_idx'] + 1
    return info

def _parse_atom_line(line_str, style_info):
    parts = line_str.split()
    try:
        atom = {'id': int(parts[style_info['id_idx']]), 'type': int(parts[style_info['type_idx']]),
                'x': float(parts[style_info['x_idx']]), 'y': float(parts[style_info['y_idx']]),
                'z': float(parts[style_info['z_idx']]), 'original_parts': parts}
        if style_info['mol_id_idx'] is not None: atom['mol_id'] = int(parts[style_info['mol_id_idx']])
        if style_info['q_idx'] is not None: atom['q'] = float(parts[style_info['q_idx']])
        return atom
    except: return None

# ==============================================================================
# _create_anchor_atom_parts
# ==============================================================================
def _create_anchor_atom_parts(original_atom_data, new_id, anchor_type_id_val, style_info):
    num_original_cols = len(original_atom_data['original_parts'])
    parts = ['0'] * num_original_cols
    parts[style_info['id_idx']] = str(new_id)
    parts[style_info['type_idx']] = str(anchor_type_id_val)
    parts[style_info['x_idx']] = "{:.15g}".format(original_atom_data['x'])
    parts[style_info['y_idx']] = "{:.15g}".format(original_atom_data['y'])
    parts[style_info['z_idx']] = "{:.15g}".format(original_atom_data['z'])
    if style_info['mol_id_idx'] is not None: parts[style_info['mol_id_idx']] = str(original_atom_data.get('mol_id', 0))
    if style_info['q_idx'] is not None: parts[style_info['q_idx']] = '0.0'
    first_extra_idx = style_info['first_extra_col_idx']
    if num_original_cols > first_extra_idx:
        for i in range(first_extra_idx, num_original_cols):
            parts[i] = original_atom_data['original_parts'][i]
    return parts

# ==============================================================================
# Main data parser
# ==============================================================================
def _parse_lammps_data_simplified(filepath):
    data = {
        'initial_comment': "", 'ordered_header_count_keys': [],
        'header_counts': {}, 'header_count_lines_map': {},
        'box_bounds_lines': [], 'box_dims': {},
        'section_headers': {}, 'sections_data_lines': {}
    }
    header_keywords = ["atoms", "bonds", "angles", "dihedrals", "impropers",
                       "atom types", "bond types", "angle types", "dihedral types", "improper types"]
    main_data_sections = ["Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"]
    coeff_section_names = ["Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", "PairIJ Coeffs"]
    all_known_section_headers = main_data_sections + coeff_section_names

    current_section_name = None
    lines = []
    with open(filepath, 'r') as f: lines = [line.strip() for line in f.readlines()]

    line_iter = iter(lines)
    try: data['initial_comment'] = next(line_iter)
    except StopIteration: return data

    for line in line_iter:
        if not line: continue
        words = line.split()
        if not words: continue

        is_processed = False
        if len(words) > 1 and words[0].isdigit():
            potential_kw = " ".join(words[1:])
            if potential_kw in header_keywords:
                if potential_kw not in data['header_counts']:
                    data['ordered_header_count_keys'].append(potential_kw)
                data['header_counts'][potential_kw] = int(words[0])
                data['header_count_lines_map'][potential_kw] = line
                current_section_name = None; is_processed = True
        if is_processed: continue

        if len(words) >= 3 and any(val in line for val in ["xlo xhi", "ylo yhi", "zlo zhi"]):
            try:
                lo, hi = float(words[0]), float(words[1])
                if "xlo xhi" in line: data['box_dims']['x'] = [lo, hi]
                elif "ylo yhi" in line: data['box_dims']['y'] = [lo, hi]
                elif "zlo zhi" in line: data['box_dims']['z'] = [lo, hi]
                data['box_bounds_lines'].append(line)
                current_section_name = None; is_processed = True
            except ValueError: pass
        if is_processed: continue

        section_header_candidate = line.split('#')[0].strip()
        if section_header_candidate in all_known_section_headers:
            current_section_name = section_header_candidate
            if current_section_name not in coeff_section_names: # Store header and prepare data list
                data['section_headers'][current_section_name] = line
                data['sections_data_lines'][current_section_name] = []
            continue

        if current_section_name and current_section_name not in coeff_section_names:
            data['sections_data_lines'][current_section_name].append(line)

    data['atom_style_info'] = None; data['atoms'] = []
    if 'Atoms' in data['section_headers'] and 'Atoms' in data['sections_data_lines']:
        data['atom_style_info'] = _determine_atom_style_info(data['section_headers']['Atoms'])
        for atom_line_str in data['sections_data_lines']['Atoms']:
            atom_obj = _parse_atom_line(atom_line_str, data['atom_style_info'])
            if atom_obj: data['atoms'].append(atom_obj)
        data['atoms'].sort(key=lambda x: x['id']) 

    for section_key, num_int_fields in [('Velocities', 1), ('Bonds', 4), ('Angles', 5), ('Dihedrals', 6), ('Impropers', 6)]:
        data[section_key.lower()] = [] 
        if section_key in data['section_headers'] and section_key in data['sections_data_lines']:
            for line_str in data['sections_data_lines'][section_key]:
                parts = line_str.split()
                try:
                    if section_key == 'Velocities': 
                         data[section_key.lower()].append([int(parts[0])] + [float(p) for p in parts[1:4]])
                    else: 
                         data[section_key.lower()].append([int(p) for p in parts[:num_int_fields]])
                except (IndexError, ValueError):
                    print("WARNING: Could not parse line in {}: {}".format(section_key, line_str))
    return data


# ==============================================================================
# Overall MD Bounds Calculation
# ==============================================================================
def _calculate_overall_md_bounds(lammps_box_dims, user_domain_min, user_domain_max, tol_val):
    overall_bounds = {}
    axes = ['x', 'y', 'z']
    for i, axis in enumerate(axes):
        l_min, l_max = lammps_box_dims[axis][0], lammps_box_dims[axis][1]
        overall_bounds[axis + '_min'] = user_domain_min[i] if user_domain_min[i] is not None else l_min
        overall_bounds[axis + '_max'] = user_domain_max[i] if user_domain_max[i] is not None else l_max
        if overall_bounds[axis + '_min'] > overall_bounds[axis + '_max'] - tol_val:
            raise ValueError("Overall MD bound {}_min ({}) >= {}_max ({}). Check domain_min/max_coords.".format(
                axis, overall_bounds[axis + '_min'], axis, overall_bounds[axis + '_max']))
    return overall_bounds

# ==============================================================================
# Atom Deletion Logic for Cut/Notch
# ==============================================================================
def _is_close_custom(a, b, rel_tol=1e-9, abs_tol=0.0): # Custom isclose for Python 2.7
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def _identify_atoms_for_removal(atoms_list, overall_md_bounds, lammps_data_box_dims,
                                in_cut_thickness, in_cut_shift,
                                in_notch_thickness, in_notch_normal, in_notch_plane,
                                in_notch_shift, in_notch_depth,
                                remove_shape, tol_val):
    atom_ids_to_delete = set()
    axes_map = {0: 'x', 1: 'y', 2: 'z'}
    npg = {} # notch_params_geom dictionary

    is_cut_active_for_md = False
    num_active_cut_dims = 0
    cut_params_by_dim = {}

    if any(ct is not None for ct in in_cut_thickness):
        for i_ax in range(3):
            if in_cut_thickness[i_ax] is not None:
                num_active_cut_dims += 1
                ct_val = abs(in_cut_thickness[i_ax])
                if ct_val < tol_val: ct_val = tol_val
                cs_val = in_cut_shift[i_ax] if in_cut_shift[i_ax] is not None else 0.0
                axis_char = axes_map[i_ax]
                center_overall_cut = (overall_md_bounds[axis_char + '_min'] + overall_md_bounds[axis_char + '_max']) / 2.0
                min_coord_cut = center_overall_cut + cs_val - ct_val / 2.0
                max_coord_cut = center_overall_cut + cs_val + ct_val / 2.0
                cut_params_by_dim[axis_char] = {'min': min_coord_cut, 'max': max_coord_cut, 'thick': ct_val, 'center': center_overall_cut + cs_val}
        if num_active_cut_dims >= 2: 
            is_cut_active_for_md = True
        elif num_active_cut_dims == 1:
            print("INFO: Cut defined in only 1D. MD atoms will not be removed for this cut.")

    is_notch_active_for_md = False
    if in_notch_thickness is not None:
        is_notch_active_for_md = True
        nt_val = abs(in_notch_thickness)
        if nt_val < tol_val: nt_val = tol_val
        
        nn_idx = in_notch_normal - 1
        np_param = in_notch_plane
        ns_val = in_notch_shift if in_notch_shift is not None else 0.0
        nd_user = in_notch_depth

        normal_ax_char = axes_map[nn_idx]
        plane_param_ax_char = axes_map[abs(np_param) - 1]
        is_plane_pos_side = np_param > 0
        
        if normal_ax_char == plane_param_ax_char:
            raise ValueError("Notch normal and plane axes cannot be the same for MD removal.")

        all_axes_chars = ['x', 'y', 'z']
        third_axis_char = [ax for ax in all_axes_chars if ax not in [normal_ax_char, plane_param_ax_char]][0]

        npg['normal_axis'] = normal_ax_char
        npg['plane_axis'] = plane_param_ax_char
        npg['third_axis'] = third_axis_char
        npg['is_plane_pos_side'] = is_plane_pos_side
        
        center_N_overall_box = (overall_md_bounds[normal_ax_char + '_min'] + overall_md_bounds[normal_ax_char + '_max']) / 2.0
        npg['shifted_center_N'] = center_N_overall_box + ns_val
        npg['half_thick_N'] = nt_val / 2.0

        surface_P_md = lammps_data_box_dims[plane_param_ax_char][1] if is_plane_pos_side else lammps_data_box_dims[plane_param_ax_char][0]
        dim_P_md = lammps_data_box_dims[plane_param_ax_char][1] - lammps_data_box_dims[plane_param_ax_char][0]
        actual_depth_for_md_removal = 0

        if nd_user is None:
            actual_depth_for_md_removal = dim_P_md / 2.0
        else:
            surface_P_overall = overall_md_bounds[plane_param_ax_char + '_max'] if is_plane_pos_side else overall_md_bounds[plane_param_ax_char + '_min']
            dim_P_overall = overall_md_bounds[plane_param_ax_char + '_max'] - overall_md_bounds[plane_param_ax_char + '_min']
            
            depth_value_from_overall_surface = abs(nd_user)
            depth_value_from_overall_surface = min(depth_value_from_overall_surface, dim_P_overall) 
            notch_tip_coord_overall = (surface_P_overall - depth_value_from_overall_surface) if is_plane_pos_side else (surface_P_overall + depth_value_from_overall_surface)

            if is_plane_pos_side:
                actual_depth_for_md_removal = surface_P_md - notch_tip_coord_overall
                actual_depth_for_md_removal = min(actual_depth_for_md_removal, surface_P_md - lammps_data_box_dims[plane_param_ax_char][0]) 
            else: 
                actual_depth_for_md_removal = notch_tip_coord_overall - surface_P_md
                actual_depth_for_md_removal = min(actual_depth_for_md_removal, lammps_data_box_dims[plane_param_ax_char][1] - surface_P_md) 
            actual_depth_for_md_removal = max(0, actual_depth_for_md_removal) 

        npg['use_rect_lead_in'] = False
        if remove_shape == 'ellipse':
            ellipse_depth = max(actual_depth_for_md_removal, tol_val * 1e-2)
            ellipse_base = surface_P_md
            
            if bridge_thickness > 0:
                npg['use_rect_lead_in'] = True
                if is_plane_pos_side:
                    npg['rect_lead_in_max_P'] = surface_P_md
                    npg['rect_lead_in_min_P'] = surface_P_md - bridge_thickness
                    ellipse_base = npg['rect_lead_in_min_P']
                else:
                    npg['rect_lead_in_min_P'] = surface_P_md
                    npg['rect_lead_in_max_P'] = surface_P_md + bridge_thickness
                    ellipse_base = npg['rect_lead_in_max_P']
                ellipse_depth = max(0, ellipse_depth - bridge_thickness)
            
            npg['ellipse_semi_axis_P'] = ellipse_depth
            npg['ellipse_center_P_on_surface'] = ellipse_base
        
        if is_plane_pos_side:
            npg['cube_plane_min'] = surface_P_md - actual_depth_for_md_removal
            npg['cube_plane_max'] = surface_P_md
        else:
            npg['cube_plane_min'] = surface_P_md
            npg['cube_plane_max'] = surface_P_md + actual_depth_for_md_removal
        
        npg['cube_normal_min'] = npg['shifted_center_N'] - npg['half_thick_N']
        npg['cube_normal_max'] = npg['shifted_center_N'] + npg['half_thick_N']
        
        npg['third_axis_min_extrude'] = lammps_data_box_dims[third_axis_char][0]
        npg['third_axis_max_extrude'] = lammps_data_box_dims[third_axis_char][1]


    for atom in atoms_list:
        ax, ay, az = atom['x'], atom['y'], atom['z']
        atom_coords_map = {'x': ax, 'y': ay, 'z': az}
        delete_this_atom = False

        if is_cut_active_for_md:
            in_cut_region_current_atom = True 
            if remove_shape == 'cube':
                for dim_char, params in cut_params_by_dim.items():
                    if not (params['min'] - tol_val <= atom_coords_map[dim_char] <= params['max'] + tol_val):
                        in_cut_region_current_atom = False; break
            elif remove_shape == 'ellipse':
                ellipse_sum = 0.0
                active_ellipse_dims = [d for d in ['x', 'y', 'z'] if d in cut_params_by_dim]
                valid_ellipse_cut = True
                for dim_char_check in active_ellipse_dims:
                    if cut_params_by_dim[dim_char_check]['thick'] / 2.0 < tol_val:
                        valid_ellipse_cut = False; break
                
                if not valid_ellipse_cut: 
                    in_cut_region_current_atom = False
                elif len(active_ellipse_dims) >= 2: 
                    for dim_char in active_ellipse_dims:
                        ellipse_sum += math.pow((atom_coords_map[dim_char] - cut_params_by_dim[dim_char]['center']) / (cut_params_by_dim[dim_char]['thick'] / 2.0), 2)
                    if ellipse_sum > 1.0 + tol_val: in_cut_region_current_atom = False
                else: 
                    in_cut_region_current_atom = False
            
            if in_cut_region_current_atom: delete_this_atom = True

        if not delete_this_atom and is_notch_active_for_md:
            in_notch_region_current_atom = False
            coord_N = atom_coords_map[npg['normal_axis']]
            coord_P = atom_coords_map[npg['plane_axis']]
            coord_T = atom_coords_map[npg['third_axis']]

            in_third_axis_bounds = (npg['third_axis_min_extrude'] - tol_val <= coord_T <= npg['third_axis_max_extrude'] + tol_val)

            if in_third_axis_bounds:
                if remove_shape == 'cube':
                    in_N_bounds_cube = (npg['cube_normal_min'] - tol_val <= coord_N <= npg['cube_normal_max'] + tol_val)
                    in_P_bounds_cube = (npg['cube_plane_min'] - tol_val <= coord_P <= npg['cube_plane_max'] + tol_val)
                    if in_N_bounds_cube and in_P_bounds_cube:
                        in_notch_region_current_atom = True

                elif remove_shape == 'ellipse':
                    if npg.get('use_rect_lead_in', False):
                        in_rect_N = (npg['cube_normal_min'] - tol_val <= coord_N <= npg['cube_normal_max'] + tol_val)
                        in_rect_P = (npg['rect_lead_in_min_P'] - tol_val <= coord_P <= npg['rect_lead_in_max_P'] + tol_val)
                        if in_rect_N and in_rect_P:
                            in_notch_region_current_atom = True
                    
                    if not in_notch_region_current_atom:
                        surface_P_md_local = npg['ellipse_center_P_on_surface']
                        semi_axis_P_depth_local = npg['ellipse_semi_axis_P']
                        in_depth_region = False
                        if npg['is_plane_pos_side']:
                            if coord_P <= surface_P_md_local + tol_val: in_depth_region = True
                        else:
                            if coord_P >= surface_P_md_local - tol_val: in_depth_region = True
                        
                        if in_depth_region:
                            semi_N_val = npg['half_thick_N']
                            center_N_val = npg['shifted_center_N']
                            
                            is_degenerate_N = semi_N_val < tol_val
                            is_degenerate_P = semi_axis_P_depth_local < tol_val

                            if is_degenerate_N and is_degenerate_P: 
                                if _is_close_custom(coord_N, center_N_val, abs_tol=tol_val) and \
                                   _is_close_custom(coord_P, surface_P_md_local, abs_tol=tol_val):
                                    in_notch_region_current_atom = True
                            elif is_degenerate_N: 
                                if _is_close_custom(coord_N, center_N_val, abs_tol=tol_val) and \
                                   (math.pow((coord_P - surface_P_md_local) / semi_axis_P_depth_local, 2) <= 1.0 + tol_val):
                                    in_notch_region_current_atom = True
                            elif is_degenerate_P: 
                                if _is_close_custom(coord_P, surface_P_md_local, abs_tol=tol_val) and \
                                   (math.pow((coord_N - center_N_val) / semi_N_val, 2) <= 1.0 + tol_val):
                                    in_notch_region_current_atom = True
                            else: 
                                term_N_sq = math.pow((coord_N - center_N_val) / semi_N_val, 2)
                                term_P_sq = math.pow((coord_P - surface_P_md_local) / semi_axis_P_depth_local, 2)
                                if (term_N_sq + term_P_sq) <= 1.0 + tol_val:
                                    in_notch_region_current_atom = True

            if in_notch_region_current_atom: delete_this_atom = True
            
        if delete_this_atom:
            atom_ids_to_delete.add(atom['id'])
            
    return atom_ids_to_delete

# ==============================================================================
# Interaction Filtering and Re-sequencing
# ==============================================================================
def _filter_interactions(interaction_data_list, atom_ids_to_delete_set, atom_id_indices_in_entry):
    kept_interactions = []
    if not interaction_data_list: return kept_interactions
    for entry in interaction_data_list:
        delete_interaction = False
        for atom_idx_pos in atom_id_indices_in_entry:
            if entry[atom_idx_pos] in atom_ids_to_delete_set:
                delete_interaction = True; break
        if not delete_interaction:
            kept_interactions.append(list(entry)) 
    return kept_interactions

def _resequence_atom_ids_and_update_parts(kept_atoms_list, atom_style_info):
    old_to_new_id_map = {}
    for new_id_minus_1, atom_obj in enumerate(kept_atoms_list):
        new_id = new_id_minus_1 + 1
        old_to_new_id_map[atom_obj['id']] = new_id
        atom_obj['id'] = new_id
        atom_obj['original_parts'][atom_style_info['id_idx']] = str(new_id) 
    return old_to_new_id_map

def _update_interaction_atom_ids(interaction_list, old_to_new_atom_id_map, atom_id_indices_in_entry):
    if not interaction_list: return
    for entry in interaction_list:
        for atom_idx_pos in atom_id_indices_in_entry:
            if entry[atom_idx_pos] in old_to_new_atom_id_map:
                 entry[atom_idx_pos] = old_to_new_atom_id_map[entry[atom_idx_pos]]

def _resequence_interaction_ids(interaction_list):
    if not interaction_list: return
    for new_id_minus_1, entry in enumerate(interaction_list):
        entry[0] = new_id_minus_1 + 1 

def _filter_and_resequence_velocities(velocities_list, atom_ids_to_delete_set, old_to_new_atom_id_map):
    kept_velocities = []
    if not velocities_list: return kept_velocities
    for vel_entry in velocities_list:
        atom_id = vel_entry[0]
        if atom_id not in atom_ids_to_delete_set:
            if atom_id in old_to_new_atom_id_map:
                new_vel_entry = [old_to_new_atom_id_map[atom_id]] + vel_entry[1:]
                kept_velocities.append(new_vel_entry)
    kept_velocities.sort(key=lambda x: x[0]) 
    return kept_velocities


# ==============================================================================
# Bond Cutting Logic
# ==============================================================================
def _cut_long_bonds_in_non_periodic_dims(parsed_data, domain_min_coords_in, domain_max_coords_in):
    atom_id_to_obj_map = {atom['id']: atom for atom in parsed_data['atoms']}
    
    non_periodic_dims_for_cutting = []
    if domain_min_coords_in[0] is not None or domain_max_coords_in[0] is not None: non_periodic_dims_for_cutting.append(0)
    if domain_min_coords_in[1] is not None or domain_max_coords_in[1] is not None: non_periodic_dims_for_cutting.append(1)
    if domain_min_coords_in[2] is not None or domain_max_coords_in[2] is not None: non_periodic_dims_for_cutting.append(2)

    if not non_periodic_dims_for_cutting: return

    box_dims = parsed_data['box_dims']
    if not all(k in box_dims for k in ['x', 'y', 'z']): return

    box_lengths = [
        box_dims['x'][1] - box_dims['x'][0],
        box_dims['y'][1] - box_dims['y'][0],
        box_dims['z'][1] - box_dims['z'][0]
    ]

    def is_pair_long(atom1, atom2):
        if not atom1 or not atom2: return False
        coords1, coords2 = [atom1['x'], atom1['y'], atom1['z']], [atom2['x'], atom2['y'], atom2['z']]
        for dim_idx in non_periodic_dims_for_cutting:
            if abs(coords1[dim_idx] - coords2[dim_idx]) > (box_lengths[dim_idx] / 2.0): return True
        return False

    interaction_specs_for_cutting = [
        {'key': 'bonds', 'indices': [2, 3]},
        {'key': 'angles', 'indices': [2, 3, 4]},
        {'key': 'dihedrals', 'indices': [2, 3, 4, 5]},
        {'key': 'impropers', 'indices': [2, 3, 4, 5]}
    ]

    for spec in interaction_specs_for_cutting:
        key, indices = spec['key'], spec['indices']
        if not parsed_data.get(key): continue

        original_count = len(parsed_data[key])
        kept_interactions = []
        for interaction in parsed_data[key]:
            atom_ids = [interaction[i] for i in indices]
            cut_this = False
            for i in range(len(atom_ids) - 1):
                atom1 = atom_id_to_obj_map.get(atom_ids[i])
                atom2 = atom_id_to_obj_map.get(atom_ids[i+1])
                if is_pair_long(atom1, atom2):
                    cut_this = True
                    break
            if not cut_this:
                kept_interactions.append(interaction)
        
        num_cut = original_count - len(kept_interactions)
        if num_cut > 0:
            print("Interaction cutting ({}): {} entries cut.".format(key.capitalize(), num_cut))
            parsed_data[key] = kept_interactions


# ==============================================================================
# Anchor region logic
# ==============================================================================
def _determine_anchor_regions(box_dims, d_min_coords, d_max_coords, buff_thick, tol_val):
    regions = []
    if not box_dims or not all(k in box_dims for k in ['x', 'y', 'z']): return regions
    if d_min_coords[0] is not None: regions.append({'dim': 0, 'type': 'min', 'range': [box_dims['x'][0] - tol_val, box_dims['x'][0] + buff_thick + tol_val]})
    if d_max_coords[0] is not None: regions.append({'dim': 0, 'type': 'max', 'range': [box_dims['x'][1] - buff_thick - tol_val, box_dims['x'][1] + tol_val]})
    if d_min_coords[1] is not None: regions.append({'dim': 1, 'type': 'min', 'range': [box_dims['y'][0] - tol_val, box_dims['y'][0] + buff_thick + tol_val]})
    if d_max_coords[1] is not None: regions.append({'dim': 1, 'type': 'max', 'range': [box_dims['y'][1] - buff_thick - tol_val, box_dims['y'][1] + tol_val]})
    if d_min_coords[2] is not None: regions.append({'dim': 2, 'type': 'min', 'range': [box_dims['z'][0] - tol_val, box_dims['z'][0] + buff_thick + tol_val]})
    if d_max_coords[2] is not None: regions.append({'dim': 2, 'type': 'max', 'range': [box_dims['z'][1] - buff_thick - tol_val, box_dims['z'][1] + tol_val]})
    return regions

def _is_atom_in_anchor_regions(atom_data, anchor_regions, tol_val):
    coords = [atom_data['x'], atom_data['y'], atom_data['z']]
    for region in anchor_regions:
        dim_coord = coords[region['dim']]
        if dim_coord >= region['range'][0] - tol_val and dim_coord <= region['range'][1] + tol_val: 
            return True
    return False

# ==============================================================================
# File writer
# ==============================================================================
def _write_section_if_present(f, section_name_key, parsed_data_view, data_lines_override=None):
    header = parsed_data_view['section_headers'].get(section_name_key)
    data_to_write = []
    if data_lines_override is not None:
        data_to_write = data_lines_override
    elif section_name_key.lower() in parsed_data_view: 
        for entry in parsed_data_view[section_name_key.lower()]:
            data_to_write.append(" ".join(map(str, entry)))
            
    if not header and not data_to_write: return 

    f.write((header if header else section_name_key) + "\n\n") 
    if data_to_write:
        for line in data_to_write: f.write(line + "\n")
    f.write("\n")


def _write_modified_lammps_file(output_filepath, parsed_data,
                                final_atom_lines, final_bond_lines,
                                final_velocity_lines, 
                                final_header_counts, anchor_atom_type_id_for_mass_sec):
    with open(output_filepath, 'w') as f:
        f.write(parsed_data['initial_comment'] + "\n\n")

        for kw in parsed_data['ordered_header_count_keys']: 
            count_val = final_header_counts.get(kw)
            if count_val is not None: f.write("{} {}\n".format(count_val, kw))
        
        original_header_keys_set = set(parsed_data['ordered_header_count_keys'])
        for kw, count_val in final_header_counts.items():
            if kw not in original_header_keys_set and count_val is not None:
                 f.write("{} {}\n".format(count_val, kw))
        f.write("\n")

        for line in parsed_data['box_bounds_lines']: f.write(line + "\n")
        f.write("\n")

        masses_header = parsed_data['section_headers'].get('Masses')
        original_mass_data_lines = parsed_data['sections_data_lines'].get('Masses', [])
        if masses_header or anchor_atom_type_id_for_mass_sec is not None:
            f.write((masses_header if masses_header else "Masses") + "\n\n")
            for mass_line in original_mass_data_lines: f.write(mass_line + "\n")
            if anchor_atom_type_id_for_mass_sec is not None:
                f.write("{} {}\n".format(anchor_atom_type_id_for_mass_sec, str(ANCHOR_ATOM_MASS_VAL)))
            f.write("\n")

        _write_section_if_present(f, "Atoms", parsed_data, data_lines_override=final_atom_lines)
        _write_section_if_present(f, "Velocities", parsed_data, data_lines_override=final_velocity_lines)
        _write_section_if_present(f, "Bonds", parsed_data, data_lines_override=final_bond_lines)

        topology_section_names = ["Angles", "Dihedrals", "Impropers"]
        for sec_name in topology_section_names:
            _write_section_if_present(f, sec_name, parsed_data)

# ==============================================================================
# Main data modification logic
# ==============================================================================
def run_anchor_cut_process():
    global domain_min_coords, domain_max_coords # Allow modification of globals
    print("="*70); print(" LAMMPS Data Modification Script ".center(70, "=")); print("="*70)
    script_dir = os.getcwd()
    try:
        input_file, _ = find_lammps_data_file(script_dir, lammps_file)
        parsed_data = _parse_lammps_data_simplified(input_file)

        if not parsed_data['box_dims'] or not all(k in parsed_data['box_dims'] for k in ['x', 'y', 'z']):
            raise ValueError("Incomplete box dimensions in LAMMPS data.")
        if not parsed_data['atoms']:
            raise ValueError("No atom data parsed from LAMMPS file.")
        
        if circular_burrito_radius is not None:
            print("INFO: Cylindrical 'burrito' mode detected. Recalculating domain center based on LAMMPS box.")
            if pbc_direction not in [1, 2, 3]: 
                raise ValueError("pbc_direction must be 1 (x), 2 (y), or 3 (z).")
            axes = ['x', 'y', 'z']
            pbc_axis = axes[pbc_direction - 1]
            non_pbc_axes = [axis for axis in axes if axis != pbc_axis]

            # In burrito mode, PBC axis bounds are None (periodic)
            # Non-PBC axes are centered on the MD box
            
            # Update the global coordinate lists (mutable lists)
            # Index map: x=0, y=1, z=2
            pbc_idx = pbc_direction - 1
            domain_min_coords[pbc_idx] = None
            domain_max_coords[pbc_idx] = None
            
            for axis in non_pbc_axes:
                ax_idx = {'x':0, 'y':1, 'z':2}[axis]
                l_min = parsed_data['box_dims'][axis][0]
                l_max = parsed_data['box_dims'][axis][1]
                center = (l_min + l_max) / 2.0
                
                domain_min_coords[ax_idx] = center - circular_burrito_radius
                domain_max_coords[ax_idx] = center + circular_burrito_radius
                
            print("  Updated domain_min_coords: {}".format(domain_min_coords))
            print("  Updated domain_max_coords: {}".format(domain_max_coords))

        print("Initial atoms: {}, bonds: {}, angles: {}, dihedrals: {}, impropers: {}".format(
            len(parsed_data['atoms']), len(parsed_data['bonds']), len(parsed_data.get('angles',[])),
            len(parsed_data.get('dihedrals',[])), len(parsed_data.get('impropers',[]))
        ))

        overall_md_bounds = _calculate_overall_md_bounds(parsed_data['box_dims'], domain_min_coords, domain_max_coords, tol)
        atom_ids_to_delete = _identify_atoms_for_removal(
            parsed_data['atoms'], overall_md_bounds, parsed_data['box_dims'],
            cut_thickness, cut_shift,
            notch_thickness, notch_normal, notch_plane, notch_shift, notch_depth,
            md_remove_shape, tol
        )
        if atom_ids_to_delete:
            print("Marked {} atoms for deletion due to cut/notch.".format(len(atom_ids_to_delete)))

        current_atoms = [atom for atom in parsed_data['atoms'] if atom['id'] not in atom_ids_to_delete]
        old_to_new_atom_id_map = _resequence_atom_ids_and_update_parts(current_atoms, parsed_data['atom_style_info'])
        parsed_data['atoms'] = current_atoms

        interaction_specs = [
            {'key': 'bonds', 'atom_indices': [2, 3]},
            {'key': 'angles', 'atom_indices': [2, 3, 4]},
            {'key': 'dihedrals', 'atom_indices': [2, 3, 4, 5]},
            {'key': 'impropers', 'atom_indices': [2, 3, 4, 5]} 
        ]
        for spec in interaction_specs:
            kept_interactions = _filter_interactions(parsed_data.get(spec['key'], []), atom_ids_to_delete, spec['atom_indices'])
            _update_interaction_atom_ids(kept_interactions, old_to_new_atom_id_map, spec['atom_indices'])
            _resequence_interaction_ids(kept_interactions)
            parsed_data[spec['key']] = kept_interactions 

        parsed_data['velocities'] = _filter_and_resequence_velocities(
            parsed_data.get('velocities', []), atom_ids_to_delete, old_to_new_atom_id_map
        )
        
        print("After cut/notch: atoms: {}, bonds: {}, angles: {}, dihedrals: {}, impropers: {}".format(
            len(parsed_data['atoms']), len(parsed_data['bonds']), len(parsed_data.get('angles',[])),
            len(parsed_data.get('dihedrals',[])), len(parsed_data.get('impropers',[]))
        ))

        if truncate_chains:
            print("truncate_chains is True. Proceeding with interaction cutting across non-periodic boundaries.")
            _cut_long_bonds_in_non_periodic_dims(parsed_data, domain_min_coords, domain_max_coords)
            for spec in interaction_specs:
                if parsed_data.get(spec['key']):
                    _resequence_interaction_ids(parsed_data[spec['key']])
        else:
            print("truncate_chains is False. Skipping interaction cutting.")

        if remove_bridge_interactions:
            print("remove_bridge_interactions is True. Removing interactions within bridge region.")
            bridge_regions_def = _determine_anchor_regions(parsed_data['box_dims'], domain_min_coords, domain_max_coords, bridge_thickness, tol)
            
            bridge_atom_ids = set()
            for atom in parsed_data['atoms']:
                if _is_atom_in_anchor_regions(atom, bridge_regions_def, tol):
                    bridge_atom_ids.add(atom['id'])
            
            if bridge_atom_ids:
                print("Identified {} atoms in bridge region. Checking interactions...".format(len(bridge_atom_ids)))
                for spec in interaction_specs:
                    key = spec['key']
                    if not parsed_data.get(key): continue
                    
                    original_count = len(parsed_data[key])
                    atom_id_indices = spec['atom_indices']
                    kept_interactions = [
                        interaction for interaction in parsed_data[key]
                        if not all(interaction[i] in bridge_atom_ids for i in atom_id_indices)
                    ]

                    if len(kept_interactions) < original_count:
                        print("  - {}: Removed {} of {} entries.".format(
                            key.capitalize(), original_count - len(kept_interactions), original_count))
                        parsed_data[key] = kept_interactions
                        _resequence_interaction_ids(parsed_data[key])
            else:
                print("No atoms found in bridge region. No interactions removed.")
        else:
            print("remove_bridge_interactions is False. Skipping bridge interaction removal.")

        if dpd_thickness is not None and dpd_thickness > 0:
            print("dpd_thickness is active. Re-typing atoms in boundary regions.")
            dpd_regions = _determine_anchor_regions(parsed_data['box_dims'], domain_min_coords, domain_max_coords, dpd_thickness, tol)
            
            atoms_in_dpd_region = []
            original_types_in_dpd = set()
            for atom in parsed_data['atoms']:
                if _is_atom_in_anchor_regions(atom, dpd_regions, tol):
                    atoms_in_dpd_region.append(atom)
                    original_types_in_dpd.add(atom['type'])
            
            if original_types_in_dpd:
                max_existing_type = max(atom['type'] for atom in parsed_data['atoms'])
                
                old_to_new_type_map = {}
                next_new_type = max_existing_type + 1
                for old_type in sorted(list(original_types_in_dpd)):
                    old_to_new_type_map[old_type] = next_new_type
                    next_new_type += 1

                mass_map = {}
                for mass_line in parsed_data['sections_data_lines'].get('Masses', []):
                    parts = mass_line.split()
                    if len(parts) >= 2 and parts[0].isdigit():
                        mass_map[int(parts[0])] = parts[1]
                
                for old_type, new_type in old_to_new_type_map.items():
                    if old_type in mass_map:
                        parsed_data['sections_data_lines']['Masses'].append("{} {}".format(new_type, mass_map[old_type]))
                
                type_idx = parsed_data['atom_style_info']['type_idx']
                for atom in atoms_in_dpd_region:
                    new_type = old_to_new_type_map[atom['type']]
                    atom['type'] = new_type
                    atom['original_parts'][type_idx] = str(new_type)
                
                print("Re-typed {} atoms into {} new types. Mapping: {}".format(
                    len(atoms_in_dpd_region), len(old_to_new_type_map), old_to_new_type_map))
            else:
                print("No atoms found in DPD regions. No types changed.")
        
        if anchor_placement_probability > 0:
            original_style_info = parsed_data['atom_style_info']
            style_name = original_style_info['style_name']
            new_style_name = None
            if 'bond' not in style_name and 'molecular' not in style_name and 'full' not in style_name:
                if 'charge' in style_name: new_style_name = 'full'
                else: new_style_name = 'bond'
            
            if new_style_name:
                print("INFO: Atom style changed from '{}' to '{}' to support anchor bonds.".format(style_name, new_style_name))
                new_atom_header_line = "Atoms # {}".format(new_style_name)
                parsed_data['section_headers']['Atoms'] = new_atom_header_line
                new_style_info = _determine_atom_style_info(new_atom_header_line)
                parsed_data['atom_style_info'] = new_style_info
                for atom in parsed_data['atoms']:
                    old_parts = list(atom['original_parts'])
                    num_old_extra = len(old_parts) - original_style_info['first_extra_col_idx']
                    if num_old_extra < 0: num_old_extra = 0
                    num_new_cols = new_style_info['first_extra_col_idx'] + num_old_extra
                    new_parts = ['0'] * num_new_cols
                    new_parts[new_style_info['id_idx']] = old_parts[original_style_info['id_idx']]
                    new_parts[new_style_info['type_idx']] = old_parts[original_style_info['type_idx']]
                    new_parts[new_style_info['x_idx']] = old_parts[original_style_info['x_idx']]
                    new_parts[new_style_info['y_idx']] = old_parts[original_style_info['y_idx']]
                    new_parts[new_style_info['z_idx']] = old_parts[original_style_info['z_idx']]
                    if new_style_info['mol_id_idx'] is not None:
                        if original_style_info['mol_id_idx'] is not None:
                            new_parts[new_style_info['mol_id_idx']] = old_parts[original_style_info['mol_id_idx']]
                        else:
                            new_mol_id = atom['id']
                            atom['mol_id'] = new_mol_id
                            new_parts[new_style_info['mol_id_idx']] = str(new_mol_id)
                    if new_style_info['q_idx'] is not None:
                        new_parts[new_style_info['q_idx']] = old_parts[original_style_info['q_idx']] if original_style_info['q_idx'] is not None else '0.0'
                    if num_old_extra > 0:
                        for i in range(num_old_extra):
                            new_parts[new_style_info['first_extra_col_idx'] + i] = old_parts[original_style_info['first_extra_col_idx'] + i]
                    atom['original_parts'] = new_parts

        anchor_regions = _determine_anchor_regions(parsed_data['box_dims'], domain_min_coords, domain_max_coords, bridge_thickness, tol)
        new_anchor_atom_objects_for_struct = []
        new_anchor_bonds_data_temp = [] 

        max_atom_id = max(atom['id'] for atom in parsed_data['atoms']) if parsed_data['atoms'] else 0
        current_anchor_atom_id_counter = max_atom_id + 1
        
        max_type_from_atoms_section = max(atom['type'] for atom in parsed_data['atoms']) if parsed_data['atoms'] else 0
        max_type_from_masses_section = 0
        original_mass_lines = parsed_data['sections_data_lines'].get('Masses', [])
        for mass_line in original_mass_lines:
            parts = mass_line.split()
            if parts and parts[0].isdigit():
                try: max_type_from_masses_section = max(max_type_from_masses_section, int(parts[0]))
                except ValueError: pass
        dynamic_anchor_atom_type_id = max(max_type_from_atoms_section, max_type_from_masses_section) + 1
        
        original_declared_bond_types_count = parsed_data['header_counts'].get('bond types', 0)
        dynamic_new_bond_type_id_for_anchors = 1 

        num_anchors_placed = 0
        newly_placed_anchor_data_for_vel = [] 

        for atom_obj in parsed_data['atoms']: 
            if _is_atom_in_anchor_regions(atom_obj, anchor_regions, tol):
                if random.random() < anchor_placement_probability:
                    if num_anchors_placed == 0: 
                        if original_declared_bond_types_count == 0:
                            dynamic_new_bond_type_id_for_anchors = 1
                        else:
                            existing_bond_types = set()
                            if parsed_data['bonds']: 
                                for b_entry in parsed_data['bonds']: existing_bond_types.add(b_entry[1])
                            
                            dynamic_new_bond_type_id_for_anchors = original_declared_bond_types_count + 1
                            while dynamic_new_bond_type_id_for_anchors in existing_bond_types:
                                dynamic_new_bond_type_id_for_anchors +=1
                    
                    anchor_parts_list = _create_anchor_atom_parts(atom_obj, current_anchor_atom_id_counter, dynamic_anchor_atom_type_id, parsed_data['atom_style_info'])
                    new_anchor_data = {
                        'id': current_anchor_atom_id_counter, 'type': dynamic_anchor_atom_type_id,
                        'x': atom_obj['x'], 'y': atom_obj['y'], 'z': atom_obj['z'],
                        'original_parts': anchor_parts_list
                    }
                    if parsed_data['atom_style_info']['mol_id_idx'] is not None:
                        new_anchor_data['mol_id'] = atom_obj.get('mol_id', 0)
                    new_anchor_atom_objects_for_struct.append(new_anchor_data)

                    new_anchor_bonds_data_temp.append([dynamic_new_bond_type_id_for_anchors,
                                                       atom_obj['id'], current_anchor_atom_id_counter])
                    newly_placed_anchor_data_for_vel.append([current_anchor_atom_id_counter, 0, 0, 0])
                    current_anchor_atom_id_counter += 1; num_anchors_placed += 1
        
        parsed_data['atoms'].extend(new_anchor_atom_objects_for_struct)

        current_max_bond_id = max(b[0] for b in parsed_data['bonds']) if parsed_data['bonds'] else 0
        for anchor_bond_type_atoms in new_anchor_bonds_data_temp:
            current_max_bond_id += 1
            parsed_data['bonds'].append([current_max_bond_id] + anchor_bond_type_atoms)
        
        if parsed_data.get('velocities') is not None or newly_placed_anchor_data_for_vel:
             if 'velocities' not in parsed_data or parsed_data['velocities'] is None : parsed_data['velocities'] = []
             parsed_data['velocities'].extend(newly_placed_anchor_data_for_vel)
             parsed_data['velocities'].sort(key=lambda x: x[0]) 

        final_header_counts = parsed_data['header_counts'].copy()
        final_header_counts['atoms'] = len(parsed_data['atoms'])
        final_header_counts['bonds'] = len(parsed_data['bonds'])
        for spec in interaction_specs[1:]: 
            final_header_counts[spec['key']] = len(parsed_data.get(spec['key'],[])) 
        
        max_final_atom_type = max(atom['type'] for atom in parsed_data['atoms']) if parsed_data['atoms'] else 0
        current_atom_types_count = final_header_counts.get('atom types', 0)
        final_header_counts['atom types'] = max(current_atom_types_count, max_final_atom_type)

        anchor_atom_type_for_mass_section = None
        if num_anchors_placed > 0:
            anchor_atom_type_for_mass_section = dynamic_anchor_atom_type_id
            
            max_bond_type_used = 0
            if parsed_data['bonds']:
                 max_bond_type_used = max(b[1] for b in parsed_data['bonds']) if parsed_data['bonds'] else 0 
            
            current_bond_types_count = final_header_counts.get('bond types', 0)
            final_header_counts['bond types'] = max(current_bond_types_count, max_bond_type_used)
        
        final_atom_lines = [" ".join(atom['original_parts']) for atom in sorted(parsed_data['atoms'], key=lambda x: x['id'])]
        final_bond_lines = [" ".join(map(str, bond)) for bond in sorted(parsed_data['bonds'], key=lambda x: x[0])] if parsed_data['bonds'] else []
        final_velocity_lines = []
        if parsed_data.get('velocities'):
            final_velocity_lines = [
                "{} 0 0 0".format(v[0]) if isinstance(v[1], int) else "{} {:.8e} {:.8e} {:.8e}".format(v[0], v[1], v[2], v[3])
                for v in sorted(parsed_data['velocities'], key=lambda x: x[0])
            ]

        output_path = os.path.join(script_dir, output_file + ".data")
        _write_modified_lammps_file(output_path, parsed_data,
                                    final_atom_lines, final_bond_lines,
                                    final_velocity_lines,
                                    final_header_counts, anchor_atom_type_for_mass_section)

        if num_anchors_placed > 0:
            print("Placed {} anchor atoms (type {}) and {} anchor bonds (type {}).".format(
                num_anchors_placed, dynamic_anchor_atom_type_id,
                len(new_anchor_bonds_data_temp), dynamic_new_bond_type_id_for_anchors))
        elif anchor_placement_probability > 0:
            print("No anchor atoms were placed (possibly due to probability or no atoms in region).")
        else:
            print("No anchor atoms were placed (anchor_placement_probability was 0).")

        print("Final atoms: {}, bonds: {}, angles: {}, dihedrals: {}, impropers: {}".format(
            len(parsed_data['atoms']), len(parsed_data['bonds']), len(parsed_data.get('angles',[])),
            len(parsed_data.get('dihedrals',[])), len(parsed_data.get('impropers',[]))
        ))
        print("Output written to: {}".format(output_path))

    except Exception as e:
        print("Error during MD data modification: {}".format(e))
        traceback.print_exc()
    finally:
        print("="*70)
        input("Process complete. Press Enter to close this window...")

# ==============================================================================
# Main execution
# ==============================================================================
if __name__ == "__main__":
    run_anchor_cut_process()