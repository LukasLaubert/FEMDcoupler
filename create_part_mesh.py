from abaqus import *
from abaqusConstants import *
import displayGroupMdbToolset as dgm
import part as part_module
import material
import assembly
import step as abaqus_step_module
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import datum
import shutil

import os
import glob
import math
import traceback
import sys
import json


# ==============================================================================
# Parameter passing
# ==============================================================================
config_filename = None
# Search for the .json file in the command-line arguments,
# as Abaqus may pass its own arguments (e.g., '-cae') first.
if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        if '.json' in arg:
            config_filename = arg
            break

# If no .json file was found in the arguments, fall back to searching the directory
if config_filename is None:
    print "INFO: No .json config file found in arguments."
    print "      Searching for a .json file in the current directory..."
    json_files = glob.glob('*.json')
    if not json_files:
        raise IOError("No .json configuration file found in the current directory.")
    config_filename = json_files[0]
    print "      Using found file: {}".format(config_filename)

with open(config_filename, 'r') as f:
    params = json.load(f)

domain_min_coords = params['domain_min_coords']
domain_max_coords = params['domain_max_coords']
circular_burrito_radius = params['circular_burrito_radius']
pbc_direction = params['pbc_direction']

bridge_thickness = params['bridge_thickness']
default_mesh_size = params['default_mesh_size']
tol = params['tol']

full_elem_bridge = params['full_elem_bridge']
bridge_cellset = params['bridge_cellset']
mesh_md_region = params['mesh_md_region']
mesh_cutout_region = params['mesh_cutout_region']
have_iface = params['have_iface']
have_padding = params['have_padding']

anchor_placement_probability = params['anchor_placement_probability']
ANCHOR_ATOM_MASS_VAL = params['ANCHOR_ATOM_MASS_VAL']
truncate_chains = params['truncate_chains']

md_remove_shape = params['md_remove_shape']

cut_width = params['cut_width']
cut_shift = params['cut_shift']

notch_width = params['notch_width']
notch_normal = params['notch_normal']
notch_plane = params['notch_plane']
notch_shift = params['notch_shift']
notch_depth = params['notch_depth']

# Convert unicode from JSON to standard strings for Abaqus API
output_file = params['output_file'].encode('utf-8')
if params['lammps_file']:
    lammps_file = params['lammps_file'].encode('utf-8')
else:
    lammps_file = None


# ==============================================================================
# find_lammps_data_file
# ==============================================================================
def find_lammps_data_file(script_dir, specific_name=None):
    file_path = None
    if specific_name:
        # If an absolute path is provided, os.path.join will correctly use it
        potential_path = os.path.join(script_dir, specific_name)
        if os.path.exists(potential_path): file_path = potential_path
        else: raise IOError("Specified LAMMPS file not found: {}".format(potential_path))
    else:
        data_files = sorted(glob.glob(os.path.join(script_dir, '*.data')))
        if not data_files: raise IOError("No '.data' file found in directory: {}".format(script_dir))
        if len(data_files) > 1: print "WARNING: Multiple '.data' files found."
        file_path = data_files[0]
    base_name = os.path.splitext(os.path.basename(file_path))[0].replace(' ', '_').replace('-', '_')
    if not base_name: base_name = "DefaultModel"
    if base_name and base_name[0].isdigit(): base_name = "_" + base_name
    print "Using LAMMPS file: {}".format(file_path)
    return file_path, base_name


# ==============================================================================
# ==============================================================================
# Mesh creation part
# ==============================================================================
# ==============================================================================

def parse_lammps_data(file_path):
    box_bounds = {}
    xlo, xhi, ylo, yhi, zlo, zhi = None, None, None, None, None, None
    found_bounds_flags = {'x': False, 'y': False, 'z': False}
    print "Reading LAMMPS file for box boundaries: {}".format(file_path)
    with open(file_path, 'r') as f:
        for line_content in f:
            line_content = line_content.strip()
            if not line_content or line_content.startswith('#'): continue
            parts_bounds = line_content.split()
            if len(parts_bounds) >= 4 and parts_bounds[2] == 'xlo' and parts_bounds[3] == 'xhi':
                try: xlo, xhi = float(parts_bounds[0]), float(parts_bounds[1]); found_bounds_flags['x'] = True
                except ValueError: print "WARNING: Could not parse xlo xhi line: '{}'".format(line_content)
            elif len(parts_bounds) >= 4 and parts_bounds[2] == 'ylo' and parts_bounds[3] == 'yhi':
                try: ylo, yhi = float(parts_bounds[0]), float(parts_bounds[1]); found_bounds_flags['y'] = True
                except ValueError: print "WARNING: Could not parse ylo yhi line: '{}'".format(line_content)
            elif len(parts_bounds) >= 4 and parts_bounds[2] == 'zlo' and parts_bounds[3] == 'zhi':
                try: zlo, zhi = float(parts_bounds[0]), float(parts_bounds[1]); found_bounds_flags['z'] = True
                except ValueError: print "WARNING: Could not parse zlo zhi line: '{}'".format(line_content)
            if all(found_bounds_flags.values()): break
    if not all(found_bounds_flags.values()) or any(b is None for b in [xlo, xhi, ylo, yhi, zlo, zhi]):
        raise ValueError("Could not parse all xlo/xhi, ylo/yhi, zlo/zhi lines.")
    box_bounds = {'xlo': xlo, 'xhi': xhi, 'ylo': ylo, 'yhi': yhi, 'zlo': zlo, 'zhi': zhi}
    print "Original LAMMPS Box Bounds: X:[{:.4f},{:.4f}], Y:[{:.4f},{:.4f}], Z:[{:.4f},{:.4f}]".format(
        box_bounds['xlo'], box_bounds['xhi'], box_bounds['ylo'], box_bounds['yhi'], box_bounds['zlo'], box_bounds['zhi'])
    return box_bounds

def calculate_dimensions(lammps_original_box_bounds, domain_min_user, domain_max_user, buffer_val, tol_calc, full_elem_bridge_flag,
                         in_cut_width, in_cut_shift,
                         in_notch_width, in_notch_normal, in_notch_plane, in_notch_shift, in_notch_depth,
                         padding_flag, padding_thickness):
    print "\nCalculating Part Dimensions and Partitions (full_elem_bridge={}, have_padding={})...".format(full_elem_bridge_flag, padding_flag)
    overall_part_bounds = {}
    md_interface_bounds = {}
    boun_set_bounds = {}

    user_min_map = {'x': domain_min_user[0], 'y': domain_min_user[1], 'z': domain_min_user[2]}
    user_max_map = {'x': domain_max_user[0], 'y': domain_max_user[1], 'z': domain_max_user[2]}

    if circular_burrito_radius is not None:
        print "  Cylindrical 'burrito' mode activated (radius={}). Overriding all domain settings.".format(circular_burrito_radius)
        if pbc_direction not in [1, 2, 3]: raise ValueError("pbc_direction must be 1 (x), 2 (y), or 3 (z).")
        axes = ['x', 'y', 'z']
        pbc_axis = axes[pbc_direction - 1]
        non_pbc_axes = [axis for axis in axes if axis != pbc_axis]

        user_min_map[pbc_axis], user_max_map[pbc_axis] = None, None
        center = {}
        for axis in non_pbc_axes:
            center[axis] = (lammps_original_box_bounds[axis+'lo'] + lammps_original_box_bounds[axis+'hi']) / 2.0
            user_min_map[axis] = center[axis] - circular_burrito_radius
            user_max_map[axis] = center[axis] + circular_burrito_radius

    physical_partition_coords = {'x': [], 'y': [], 'z': []}

    for i, axis in enumerate(['x', 'y', 'z']):
        l_box_min, l_box_max = lammps_original_box_bounds[axis+'lo'], lammps_original_box_bounds[axis+'hi']

        if circular_burrito_radius is None:
            if user_min_map[axis] is not None and user_min_map[axis] > l_box_min - tol_calc:
                raise ValueError("{}-Dir: domain_min_coord ({:.4f}) is INSIDE or EQUAL to LAMMPS box min ({:.4f}). Must be strictly outside or None.".format(axis.upper(), user_min_map[axis], l_box_min))
            if user_max_map[axis] is not None and user_max_map[axis] < l_box_max + tol_calc:
                raise ValueError("{}-Dir: domain_max_coord ({:.4f}) is INSIDE or EQUAL to LAMMPS box max ({:.4f}). Must be strictly outside or None.".format(axis.upper(), user_max_map[axis], l_box_max))

        part_min_coord = user_min_map[axis] if user_min_map[axis] is not None else l_box_min
        part_max_coord = user_max_map[axis] if user_max_map[axis] is not None else l_box_max

        boun_set_bounds[axis + '_min'] = part_min_coord
        boun_set_bounds[axis + '_max'] = part_max_coord

        if padding_flag:
            if user_min_map[axis] is None:
                part_min_coord -= padding_thickness
                physical_partition_coords[axis].append(l_box_min)
            if user_max_map[axis] is None:
                part_max_coord += padding_thickness
                physical_partition_coords[axis].append(l_box_max)

        if part_min_coord > part_max_coord - tol_calc:
             raise ValueError("{}-Dir Final Part Min ({:.4f}) >= Max ({:.4f}). Check logic.".format(axis.upper(), part_min_coord, part_max_coord))

        overall_part_bounds[axis + '_min'] = part_min_coord
        overall_part_bounds[axis + '_max'] = part_max_coord

        md_iface_min = l_box_min
        md_iface_max = l_box_max

        if user_min_map[axis] is not None:
            md_iface_min = l_box_min + buffer_val
        if user_max_map[axis] is not None:
            md_iface_max = l_box_max - buffer_val

        if md_iface_min >= md_iface_max - tol_calc:
            print "WARNING: {}-Dir: Buffer ({:.2f}) made effective MD region invalid (min {:.4f} >= max {:.4f}). MD region zero or negative.".format(
                axis.upper(), buffer_val, md_iface_min, md_iface_max)

        md_interface_bounds[axis+'lo'] = md_iface_min
        md_interface_bounds[axis+'hi'] = md_iface_max

        current_phys_partitions = [part_min_coord, part_max_coord]

        if md_iface_min > part_min_coord + tol_calc and md_iface_min < part_max_coord - tol_calc:
            current_phys_partitions.append(md_iface_min)
        if md_iface_max < part_max_coord - tol_calc and md_iface_max > part_min_coord + tol_calc:
            if not (abs(md_iface_max - md_iface_min) < tol_calc and md_iface_min in current_phys_partitions):
                current_phys_partitions.append(md_iface_max)

        if full_elem_bridge_flag:
            if user_min_map[axis] is not None:
                if l_box_min > part_min_coord + tol_calc and l_box_min < part_max_coord - tol_calc:
                    is_new_coord = True
                    for existing_coord in current_phys_partitions:
                        if abs(l_box_min - existing_coord) < tol_calc: is_new_coord = False; break
                    if is_new_coord: current_phys_partitions.append(l_box_min)

            if user_max_map[axis] is not None:
                if l_box_max < part_max_coord - tol_calc and l_box_max > part_min_coord + tol_calc:
                    is_new_coord = True
                    for existing_coord in current_phys_partitions:
                        if abs(l_box_max - existing_coord) < tol_calc: is_new_coord = False; break
                    if is_new_coord: current_phys_partitions.append(l_box_max)
        physical_partition_coords[axis].extend(current_phys_partitions)

    is_cut_active = any(ct is not None for ct in in_cut_width)
    is_notch_active = in_notch_width is not None

    if is_cut_active:
        print "  Adding cut partitions to physical_partition_coords..."
        for i_axis, axis_char_loop in enumerate(['x', 'y', 'z']):
            if in_cut_width[i_axis] is not None:
                ct_val = abs(in_cut_width[i_axis])
                if ct_val < tol_calc: ct_val = tol_calc

                cs_val = in_cut_shift[i_axis] if in_cut_shift[i_axis] is not None else 0.0

                part_min = overall_part_bounds[axis_char_loop + '_min']
                part_max = overall_part_bounds[axis_char_loop + '_max']
                center = (part_min + part_max) / 2.0

                p1 = center + cs_val - ct_val / 2.0
                p2 = center + cs_val + ct_val / 2.0

                if p1 > part_min + tol_calc and p1 < part_max - tol_calc:
                    physical_partition_coords[axis_char_loop].append(p1)
                if p2 > part_min + tol_calc and p2 < part_max - tol_calc:
                    physical_partition_coords[axis_char_loop].append(p2)

    if is_notch_active:
        print "  Adding notch partitions to physical_partition_coords..."
        n_thick_val = abs(in_notch_width)
        if n_thick_val < tol_calc: n_thick_val = tol_calc
        n_normal_idx = in_notch_normal - 1
        n_plane_param = in_notch_plane
        n_shift_val = in_notch_shift if in_notch_shift is not None else 0.0
        n_depth_user = in_notch_depth
        axes_chars = ['x', 'y', 'z']
        normal_ax_char = axes_chars[n_normal_idx]
        plane_ax_idx = abs(n_plane_param) - 1
        plane_ax_char = axes_chars[plane_ax_idx]
        is_plane_pos_side = n_plane_param > 0

        if normal_ax_char == plane_ax_char:
            raise ValueError("Notch normal ({}) and plane ({}) axes cannot be the same.".format(normal_ax_char, plane_ax_char))

        part_min_n = overall_part_bounds[normal_ax_char + '_min']
        part_max_n = overall_part_bounds[normal_ax_char + '_max']
        center_n = (part_min_n + part_max_n) / 2.0
        t_p1 = center_n + n_shift_val - n_thick_val / 2.0
        t_p2 = center_n + n_shift_val + n_thick_val / 2.0
        if t_p1 > part_min_n + tol_calc and t_p1 < part_max_n - tol_calc:
            physical_partition_coords[normal_ax_char].append(t_p1)
        if t_p2 > part_min_n + tol_calc and t_p2 < part_max_n - tol_calc:
            physical_partition_coords[normal_ax_char].append(t_p2)

        part_min_p = overall_part_bounds[plane_ax_char + '_min']
        part_max_p = overall_part_bounds[plane_ax_char + '_max']
        surface_coord = part_max_p if is_plane_pos_side else part_min_p
        actual_depth = abs(part_max_p - part_min_p) / 2.0 if n_depth_user is None else abs(n_depth_user)
        if actual_depth > abs(part_max_p - part_min_p) - tol_calc : actual_depth = abs(part_max_p - part_min_p)

        depth_coord = surface_coord - actual_depth if is_plane_pos_side else surface_coord + actual_depth

        valid_depth_partition = False
        if is_plane_pos_side and depth_coord < surface_coord - tol_calc and depth_coord > part_min_p + tol_calc:
            valid_depth_partition = True
        elif not is_plane_pos_side and depth_coord > surface_coord + tol_calc and depth_coord < part_max_p - tol_calc:
            valid_depth_partition = True

        if valid_depth_partition:
            physical_partition_coords[plane_ax_char].append(depth_coord)
        else:
            print "    WARNING: Notch depth partition for {}-axis results in non-internal or zero-depth partition. Skipping this depth coord.".format(plane_ax_char.upper())

    for axis in ['x', 'y', 'z']:
        coords_list = physical_partition_coords[axis]
        coords_list.sort()
        unique_coords_axis = [coords_list[0]] if coords_list else []
        for j_coord in range(1, len(coords_list)):
            if abs(coords_list[j_coord] - unique_coords_axis[-1]) > tol_calc:
                unique_coords_axis.append(coords_list[j_coord])
        physical_partition_coords[axis] = unique_coords_axis
        print "  {}-Dir Physical Partitions (incl. cut/notch/pad): {}".format(axis.upper(), ["{:.4f}".format(c) for c in unique_coords_axis])

    print "MD Interface Bounds (for classification, after buffer): X:[{:.4f},{:.4f}], Y:[{:.4f},{:.4f}], Z:[{:.4f},{:.4f}]".format(
        md_interface_bounds['xlo'], md_interface_bounds['xhi'],
        md_interface_bounds['ylo'], md_interface_bounds['yhi'],
        md_interface_bounds['zlo'], md_interface_bounds['zhi'])
    print "Logical BOUN Set Bounding Box: X:[{:.4f},{:.4f}], Y:[{:.4f},{:.4f}], Z:[{:.4f},{:.4f}]".format(
        boun_set_bounds['x_min'], boun_set_bounds['x_max'],
        boun_set_bounds['y_min'], boun_set_bounds['y_max'],
        boun_set_bounds['z_min'], boun_set_bounds['z_max'])
    print "Overall Part Bounding Box: X:[{:.4f},{:.4f}], Y:[{:.4f},{:.4f}], Z:[{:.4f},{:.4f}]".format(
        overall_part_bounds['x_min'], overall_part_bounds['x_max'],
        overall_part_bounds['y_min'], overall_part_bounds['y_max'],
        overall_part_bounds['z_min'], overall_part_bounds['z_max'])
    return overall_part_bounds, physical_partition_coords, md_interface_bounds, boun_set_bounds

def apply_partitions_to_part(part_obj, physical_partition_coords_dict, overall_bounds_dict, tol_geom):
    print "  Applying partitions to part '{}'...".format(part_obj.name)
    all_datum_ids = []
    z_min_part_solid = overall_bounds_dict['z_min']
    z_max_part_solid = overall_bounds_dict['z_max']

    for axis, principal_plane in [('x', YZPLANE), ('y', XZPLANE), ('z', XYPLANE)]:
        for coord in physical_partition_coords_dict[axis][1:-1]:
            if axis == 'z':
                if not (z_min_part_solid + tol_geom < coord < z_max_part_solid - tol_geom):
                    continue
            try:
                d = part_obj.DatumPlaneByPrincipalPlane(principalPlane=principal_plane, offset=coord)
                all_datum_ids.append(d.id)
            except Exception as e_datum: print "    WARNING: Could not create {} Datum Plane at {:.4f}. Err: {}".format(axis.upper(), coord, e_datum)

    if all_datum_ids:
        print "  Partitioning cells using {} Datum Planes...".format(len(all_datum_ids))
        if not part_obj.cells:
            print "    WARNING: Part has no cells to partition, but datum planes were created. Skipping partitioning."
        else:
            for dp_id in all_datum_ids:
                try:
                    part_obj.PartitionCellByDatumPlane(datumPlane=part_obj.datums[dp_id], cells=part_obj.cells)
                except AbaqusException as e_part: print "    INFO: Partition message for Datum ID {}: {}".format(dp_id, e_part)
                except Exception as e_part_gen: print "    ERROR: Unexpected partition error Datum ID {}. Err: {}".format(dp_id, e_part_gen)
    else: print "  No internal datum planes needed for partitioning based on physical_partition_coords."

    print "  Part '{}' has {} cells after partitioning.".format(part_obj.name, len(part_obj.cells))
    if not part_obj.cells and all_datum_ids : raise RuntimeError("Partitioning resulted in zero cells for part '{}'.".format(part_obj.name))

def create_partitioned_part(model_obj, part_name_str, overall_bounds_dict, physical_partition_coords_dict, tol_geom):
    print "\nCreating Abaqus part '{}'...".format(part_name_str)
    if part_name_str in model_obj.parts:
        print "  Part '{}' already exists. Deleting and recreating.".format(part_name_str)
        del model_obj.parts[part_name_str]
    part_obj = model_obj.Part(name=part_name_str, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    if circular_burrito_radius is not None:
        print "  Creating cylindrical part base (radius={:.4f}).".format(circular_burrito_radius)
        axes_map = {1: ('x', YZPLANE, YAXIS), 2: ('y', XZPLANE, ZAXIS), 3: ('z', XYPLANE, YAXIS)}
        pbc_ax_char, sketch_plane, sketch_up_axis = axes_map[pbc_direction]
        
        pbc_ax_min = overall_bounds_dict[pbc_ax_char + '_min']
        extrusion_depth = overall_bounds_dict[pbc_ax_char + '_max'] - pbc_ax_min
        if extrusion_depth < tol_geom: raise ValueError("Cylindrical extrusion depth ({:.4f}) too small.".format(extrusion_depth))

        # Calculate explicit global center coordinates
        cx = (overall_bounds_dict['x_min'] + overall_bounds_dict['x_max']) / 2.0
        cy = (overall_bounds_dict['y_min'] + overall_bounds_dict['y_max']) / 2.0
        cz = (overall_bounds_dict['z_min'] + overall_bounds_dict['z_max']) / 2.0

        # Determine sketch origin and center based on PBC direction mapping
        origin = (0, 0, 0)
        sketch_center = (0, 0)

        if pbc_direction == 1: # X-axis cylinder (Plane YZ, Up Y)
            # YZ Plane: Normal +X. Up +Y. Right is -Z.
            # Local X (Right) -> Global -Z
            # Local Y (Up)    -> Global Y
            origin = (pbc_ax_min, 0, 0)
            sketch_center = (-cz, cy) 
        elif pbc_direction == 2: # Y-axis cylinder (Plane XZ, Up Z)
            # XZ Plane: Normal +Y. Up +Z. Right is -X.
            # Local X (Right) -> Global -X
            # Local Y (Up)    -> Global Z
            origin = (0, pbc_ax_min, 0)
            sketch_center = (-cx, cz)
        elif pbc_direction == 3: # Z-axis cylinder (Plane XY, Up Y)
            # XY Plane: Normal +Z. Up +Y. Right is +X.
            # Local X (Right) -> Global X
            # Local Y (Up)    -> Global Y
            origin = (0, 0, pbc_ax_min)
            sketch_center = (cx, cy)

        temp_plane = part_obj.DatumPlaneByPrincipalPlane(principalPlane=sketch_plane, offset=pbc_ax_min)
        temp_axis = part_obj.DatumAxisByPrincipalAxis(principalAxis=sketch_up_axis)
        
        # Note: Using RIGHT orientation assumes standard Abaqus plane conventions
        transform = part_obj.MakeSketchTransform(sketchPlane=part_obj.datums[temp_plane.id], sketchUpEdge=part_obj.datums[temp_axis.id],
                                                 sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=origin)
                                                 
        sketch_name = '{}_profile'.format(part_name_str)
        if sketch_name in model_obj.sketches: del model_obj.sketches[sketch_name]
        s = model_obj.ConstrainedSketch(name=sketch_name, sheetSize=max(extrusion_depth, circular_burrito_radius*2.5), transform=transform)
        s.CircleByCenterPerimeter(center=sketch_center, point1=(sketch_center[0] + circular_burrito_radius, sketch_center[1]))
        part_obj.BaseSolidExtrude(sketch=s, depth=extrusion_depth)
        del model_obj.sketches[sketch_name]
    else:
        z_min_part = overall_bounds_dict['z_min']
        extrusion_depth = overall_bounds_dict['z_max'] - z_min_part
        if extrusion_depth < tol_geom: raise ValueError("Extrusion depth ({:.4f}) too small.".format(extrusion_depth))

        max_dim_overall = max(abs(overall_bounds_dict['x_max']-overall_bounds_dict['x_min']),
                               abs(overall_bounds_dict['y_max']-overall_bounds_dict['y_min']),
                               extrusion_depth, 1.0) * 1.2
        sheet_size_to_use = max_dim_overall if max_dim_overall > tol_geom else 200.0

        temp_sketch_plane_datum = part_obj.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=z_min_part)
        temp_sketch_axis_datum = part_obj.DatumAxisByPrincipalAxis(principalAxis=YAXIS)

        transform_to_use = part_obj.MakeSketchTransform(
            sketchPlane=part_obj.datums[temp_sketch_plane_datum.id],
            sketchUpEdge=part_obj.datums[temp_sketch_axis_datum.id],
            sketchPlaneSide=SIDE1,
            sketchOrientation=RIGHT,
            origin=(0.0, 0.0, z_min_part)
        )
        sketch_name = '{}_profile'.format(part_name_str)
        if sketch_name in model_obj.sketches:
            del model_obj.sketches[sketch_name]

        s_profile_transformed = model_obj.ConstrainedSketch(name=sketch_name,
                                                            sheetSize=sheet_size_to_use,
                                                            transform=transform_to_use)

        part_obj.projectReferencesOntoSketch(sketch=s_profile_transformed, filter=COPLANAR_EDGES)
        s_profile_transformed.rectangle(point1=(overall_bounds_dict['x_min'], overall_bounds_dict['y_min']),
                                        point2=(overall_bounds_dict['x_max'], overall_bounds_dict['y_max']))
        part_obj.BaseSolidExtrude(sketch=s_profile_transformed, depth=extrusion_depth)
        del model_obj.sketches[sketch_name]
        print "  Base part created. Actual Z range of solid: [{:.4f},{:.4f}]".format(z_min_part, z_min_part + extrusion_depth)

    if part_obj.cells:
        apply_partitions_to_part(part_obj, physical_partition_coords_dict, overall_bounds_dict, tol_geom)
    else:
        print "  WARNING: Base part extrusion resulted in no cells. Skipping partitioning."
        has_internal_partitions = False
        for axis in ['x', 'y', 'z']:
            if len(physical_partition_coords_dict[axis]) > 2:
                has_internal_partitions = True; break
        if has_internal_partitions:
             raise RuntimeError("Base part extrusion resulted in zero cells, but internal partitions were expected for part '{}'.".format(part_obj.name))

    return part_obj

# --- HELPER FUNCTIONS FOR CELL CLASSIFICATION (CUT/NOTCH) ---
def is_centroid_in_cut_region(cx, cy, cz, overall_bounds, cut_thick_list, cut_shift_list, tol_c):
    is_within_all_defined_cuts = True
    at_least_one_cut_dim_defined = False
    for i_ax, ax_char in enumerate(['x', 'y', 'z']):
        if cut_thick_list[i_ax] is not None:
            at_least_one_cut_dim_defined = True
            thk = abs(cut_thick_list[i_ax])
            if thk < tol_c: thk = tol_c
            shf = cut_shift_list[i_ax] if cut_shift_list[i_ax] is not None else 0.0
            c_coord_val = [cx, cy, cz][i_ax]

            center_ax = (overall_bounds[ax_char + '_min'] + overall_bounds[ax_char + '_max']) / 2.0
            p1_cut = center_ax + shf - thk / 2.0
            p2_cut = center_ax + shf + thk / 2.0

            if not (p1_cut - tol_c <= c_coord_val <= p2_cut + tol_c):
                is_within_all_defined_cuts = False
                break
    return at_least_one_cut_dim_defined and is_within_all_defined_cuts

def is_centroid_in_notch_region(cx, cy, cz, overall_bounds, n_thick_val, n_normal_idx, n_plane_param, n_shift_val, n_depth_user_val, tol_c):
    if n_thick_val is None: return False

    eff_n_thick = abs(n_thick_val)
    if eff_n_thick < tol_c: eff_n_thick = tol_c

    eff_n_shift = n_shift_val if n_shift_val is not None else 0.0

    axes_chars = ['x', 'y', 'z']
    normal_ax_char = axes_chars[n_normal_idx]
    c_normal = [cx, cy, cz][n_normal_idx]

    plane_ax_idx = abs(n_plane_param) - 1
    plane_ax_char = axes_chars[plane_ax_idx]
    c_plane = [cx, cy, cz][plane_ax_idx]
    is_plane_pos_side = n_plane_param > 0

    if normal_ax_char == plane_ax_char: return False

    part_min_n, part_max_n = overall_bounds[normal_ax_char+'_min'], overall_bounds[normal_ax_char+'_max']
    center_n = (part_min_n + part_max_n) / 2.0
    t_p1 = center_n + eff_n_shift - eff_n_thick / 2.0
    t_p2 = center_n + eff_n_shift + eff_n_thick / 2.0
    if not (t_p1 - tol_c <= c_normal <= t_p2 + tol_c):
        return False

    part_min_p, part_max_p = overall_bounds[plane_ax_char+'_min'], overall_bounds[plane_ax_char+'_max']
    surface_coord = part_max_p if is_plane_pos_side else part_min_p

    actual_depth = abs(part_max_p - part_min_p) / 2.0 if n_depth_user_val is None else abs(n_depth_user_val)
    if actual_depth > abs(part_max_p - part_min_p) - tol_c : actual_depth = abs(part_max_p - part_min_p)

    depth_limit_coord = surface_coord - actual_depth if is_plane_pos_side else surface_coord + actual_depth

    if is_plane_pos_side:
        if not (depth_limit_coord - tol_c <= c_plane <= surface_coord + tol_c): return False
    else:
        if not (surface_coord - tol_c <= c_plane <= depth_limit_coord + tol_c): return False

    return True

def create_part_geometry_sets(part_obj, md_iface_bounds, physical_partition_coords_dict, tol_set, bridge_cellset_flag, lammps_original_box_bounds_for_sets, overall_part_bounds_for_empty_check,in_cut_width, in_cut_shift, in_notch_width, in_notch_normal, in_notch_plane, in_notch_shift, in_notch_depth, padding_flag, boun_set_bounds):
    print "\nCreating Part Cell Sets (EMPTY, MD, FE, PAD if applicable) for '{}'...".format(part_obj.name)
    md_xlo, md_xhi = md_iface_bounds['xlo'], md_iface_bounds['xhi']
    md_ylo, md_yhi = md_iface_bounds['ylo'], md_iface_bounds['yhi']
    md_zlo, md_zhi = md_iface_bounds['zlo'], md_iface_bounds['zhi']

    empty_cells, md_cells, fe_cells, pad_cells = [], [], [], []
    
    if not part_obj.cells:
        print "  ERROR: Part '{}' has no cells for set creation.".format(part_obj.name); return

    is_cut_globally_active = any(ct is not None for ct in in_cut_width)
    is_notch_globally_active = in_notch_width is not None
    user_min_map = {'x': domain_min_coords[0], 'y': domain_min_coords[1], 'z': domain_min_coords[2]}
    user_max_map = {'x': domain_max_coords[0], 'y': domain_max_coords[1], 'z': domain_max_coords[2]}
    
    if circular_burrito_radius is not None:
        axes = ['x', 'y', 'z']
        pbc_axis = axes[pbc_direction - 1]
        user_min_map[pbc_axis], user_max_map[pbc_axis] = None, None

    # --- ROBUST CLASSIFICATION LOOP ---
    for cell_to_classify in part_obj.cells:
        point_inside = cell_to_classify.pointOn[0]
        cx, cy, cz = point_inside[0], point_inside[1], point_inside[2]

        # 1. Check for EMPTY cells (cuts/notches) first.
        is_empty = False
        if is_cut_globally_active:
            if is_centroid_in_cut_region(cx, cy, cz, overall_part_bounds_for_empty_check, in_cut_width, in_cut_shift, tol_set):
                is_empty = True
        if not is_empty and is_notch_globally_active:
            if is_centroid_in_notch_region(cx, cy, cz, overall_part_bounds_for_empty_check, in_notch_width, in_notch_normal-1, in_notch_plane, in_notch_shift, in_notch_depth, tol_set):
                is_empty = True
        if is_empty:
            empty_cells.append(cell_to_classify)
            continue

        # 2. Check for MD cells.
        is_in_pure_md = (md_xlo - tol_set <= cx <= md_xhi + tol_set and
                         md_ylo - tol_set <= cy <= md_yhi + tol_set and
                         md_zlo - tol_set <= cz <= md_zhi + tol_set)
        if is_in_pure_md:
            md_cells.append(cell_to_classify)
            continue

        # 3. Classify remaining cells as PAD or FE. The bridge region is initially part of FE.
        is_pad = False
        if padding_flag:
            l_orig_xlo, l_orig_xhi = lammps_original_box_bounds_for_sets['xlo'], lammps_original_box_bounds_for_sets['xhi']
            l_orig_ylo, l_orig_yhi = lammps_original_box_bounds_for_sets['ylo'], lammps_original_box_bounds_for_sets['yhi']
            l_orig_zlo, l_orig_zhi = lammps_original_box_bounds_for_sets['zlo'], lammps_original_box_bounds_for_sets['zhi']
            if circular_burrito_radius is not None:
                coords = [cx, cy, cz]
                l_mins = [l_orig_xlo, l_orig_ylo, l_orig_zlo]
                l_maxs = [l_orig_xhi, l_orig_yhi, l_orig_zhi]
                pbc_idx = pbc_direction - 1
                if (coords[pbc_idx] < l_mins[pbc_idx] - tol_set or coords[pbc_idx] > l_maxs[pbc_idx] + tol_set):
                    is_pad = True
            else:
                coords = [cx, cy, cz]
                for ax_idx, ax_char in enumerate(['x', 'y', 'z']):
                    if user_min_map[ax_char] is None and coords[ax_idx] < boun_set_bounds[ax_char+'_min'] - tol_set: is_pad = True; break
                    if user_max_map[ax_char] is None and coords[ax_idx] > boun_set_bounds[ax_char+'_max'] + tol_set: is_pad = True; break
        
        if is_pad:
            pad_cells.append(cell_to_classify)
        else:
            fe_cells.append(cell_to_classify)

    # Create the sets from the classified lists.
    if 'EMPTY' in part_obj.sets: del part_obj.sets['EMPTY']
    if 'MD' in part_obj.sets: del part_obj.sets['MD']
    if 'FE' in part_obj.sets: del part_obj.sets['FE']
    if 'PAD' in part_obj.sets: del part_obj.sets['PAD']
    if 'BD' in part_obj.sets: del part_obj.sets['BD'] # Always remove old BD set

    if empty_cells: part_obj.Set(cells=part_module.CellArray(cells=empty_cells), name='EMPTY'); print "  Created Part set 'EMPTY' ({} cells)".format(len(empty_cells))
    else: print "  No cells classified as EMPTY."
    if md_cells: part_obj.Set(cells=part_module.CellArray(cells=md_cells), name='MD'); print "  Created Part set 'MD' ({} cells)".format(len(md_cells))
    else: print "  No cells classified as MD."
    if fe_cells: part_obj.Set(cells=part_module.CellArray(cells=fe_cells), name='FE'); print "  Created Part set 'FE' ({} cells)".format(len(fe_cells))
    else: print "  No cells classified as FE."
    if pad_cells: part_obj.Set(cells=part_module.CellArray(cells=pad_cells), name='PAD'); print "  Created Part set 'PAD' ({} cells)".format(len(pad_cells))
    else: print "  No cells classified as PAD."

def create_bd_element_set(part_obj, lammps_original_box_bounds, tol_set):
    """Creates a BD element set by finding all elements fully inside the LAMMPS box."""
    print "\nCreating 'BD' element set based on 'fully covered' elements..."
    if 'FE' not in part_obj.sets or not part_obj.sets['FE'].elements:
        print "  'FE' set has no elements. Cannot create 'BD' set."
        return
    
    # Use a relaxed tolerance for this check to handle floating point issues at boundaries
    relaxed_tol = 1e-3
        
    l_xlo, l_xhi = lammps_original_box_bounds['xlo'], lammps_original_box_bounds['xhi']
    l_ylo, l_yhi = lammps_original_box_bounds['ylo'], lammps_original_box_bounds['yhi']
    l_zlo, l_zhi = lammps_original_box_bounds['zlo'], lammps_original_box_bounds['zhi']
    
    bd_element_labels = []
    fe_element_labels = []
    
    all_fe_elements = part_obj.sets['FE'].elements
    
    for elem in all_fe_elements:
        is_fully_covered = True
        nodes = elem.getNodes()
        if not nodes:
            is_fully_covered = False
        else:
            for node in nodes:
                nx, ny, nz = node.coordinates
                if not (l_xlo - relaxed_tol <= nx <= l_xhi + relaxed_tol and
                        l_ylo - relaxed_tol <= ny <= l_yhi + relaxed_tol and
                        l_zlo - relaxed_tol <= nz <= l_zhi + relaxed_tol):
                    is_fully_covered = False
                    break
        
        if is_fully_covered:
            bd_element_labels.append(elem.label)
        else:
            fe_element_labels.append(elem.label)
            
    # Create new BD element set
    if bd_element_labels:
        bd_elements = part_obj.elements.sequenceFromLabels(labels=bd_element_labels)
        part_obj.Set(elements=bd_elements, name='BD')
        print "  Created 'BD' element set with {} fully covered elements.".format(len(bd_elements))
    else:
        print "  No fully covered elements found to create 'BD' set."

    # Overwrite the old FE cell set with a new, smaller FE element set
    if 'FE' in part_obj.sets:
        del part_obj.sets['FE']
        
    if fe_element_labels:
        fe_elements = part_obj.elements.sequenceFromLabels(labels=fe_element_labels)
        part_obj.Set(elements=fe_elements, name='FE')
        print "  Recreated 'FE' element set with remaining {} elements.".format(len(fe_elements))
    else:
        print "  No remaining elements for the 'FE' set."


def create_part_outer_boundary_sets(part_obj, overall_bounds_part, boun_set_bounds, tol_boun):
    if circular_burrito_radius is not None:
        print "\nCreating Part Outer Boundary Face Sets (Cylindrical) for '{}'...".format(part_obj.name)
        if not part_obj.faces: print "  ERROR: Part '{}' has no faces.".format(part_obj.name); return

        for set_name in list(part_obj.sets.keys()):
            if 'BOUN_' in set_name: del part_obj.sets[set_name]

        axes = ['x', 'y', 'z']
        pbc_axis_char = axes[pbc_direction - 1]

        pbc_ax_min_boun, pbc_ax_max_boun = boun_set_bounds[pbc_axis_char + '_min'], boun_set_bounds[pbc_axis_char + '_max']
        min_faces = part_obj.faces.getByBoundingBox(**{pbc_axis_char+'Min': pbc_ax_min_boun - tol_boun, pbc_axis_char+'Max': pbc_ax_min_boun + tol_boun})
        max_faces = part_obj.faces.getByBoundingBox(**{pbc_axis_char+'Min': pbc_ax_max_boun - tol_boun, pbc_axis_char+'Max': pbc_ax_max_boun + tol_boun})
        min_set_name, max_set_name = 'BOUN_' + pbc_axis_char + 'lo', 'BOUN_' + pbc_axis_char + 'hi'

        if min_faces:
            part_obj.Set(faces=min_faces, name=min_set_name)
            print "  Created Part set '{}' ({} faces)".format(min_set_name, len(min_faces))
        else:
            print "  WARNING: No faces for Part set '{}'. Creating empty.".format(min_set_name)
            part_obj.Set(faces=part_module.FaceArray(faces=[]), name=min_set_name)
        
        if max_faces:
            part_obj.Set(faces=max_faces, name=max_set_name)
            print "  Created Part set '{}' ({} faces)".format(max_set_name, len(max_faces))
        else:
            print "  WARNING: No faces for Part set '{}'. Creating empty.".format(max_set_name)
            part_obj.Set(faces=part_module.FaceArray(faces=[]), name=max_set_name)
    else:
        print "\nCreating Part Outer Boundary Face Sets (BOUN_...) for '{}'...".format(part_obj.name)
        if not part_obj.faces: print "  ERROR: Part '{}' has no faces.".format(part_obj.name); return

        x_min_boun, x_max_boun = boun_set_bounds['x_min'], boun_set_bounds['x_max']
        y_min_boun, y_max_boun = boun_set_bounds['y_min'], boun_set_bounds['y_max']
        z_min_boun, z_max_boun = boun_set_bounds['z_min'], boun_set_bounds['z_max']

        x_min_search, x_max_search = overall_bounds_part['x_min'], overall_bounds_part['x_max']
        y_min_search, y_max_search = overall_bounds_part['y_min'], overall_bounds_part['y_max']
        z_min_search, z_max_search = overall_bounds_part['z_min'], overall_bounds_part['z_max']

        boundary_defs = {'BOUN_xlo': ('xMin', x_min_boun), 'BOUN_xhi': ('xMax', x_max_boun),
                         'BOUN_ylo': ('yMin', y_min_boun), 'BOUN_yhi': ('yMax', y_max_boun),
                         'BOUN_zlo': ('zMin', z_min_boun), 'BOUN_zhi': ('zMax', z_max_boun)}

        for set_name, (key, val) in boundary_defs.items():
            if set_name in part_obj.sets: del part_obj.sets[set_name]
            faces_found = None
            if key == 'xMin' or key == 'xMax':
                faces_found = part_obj.faces.getByBoundingBox(xMin=val-tol_boun, xMax=val+tol_boun,
                                                              yMin=y_min_search-tol_boun, yMax=y_max_search+tol_boun,
                                                              zMin=z_min_search-tol_boun, zMax=z_max_search+tol_boun)
            elif key == 'yMin' or key == 'yMax':
                faces_found = part_obj.faces.getByBoundingBox(yMin=val-tol_boun, yMax=val+tol_boun,
                                                              xMin=x_min_search-tol_boun, xMax=x_max_search+tol_boun,
                                                              zMin=z_min_search-tol_boun, zMax=z_max_search+tol_boun)
            elif key == 'zMin' or key == 'zMax':
                faces_found = part_obj.faces.getByBoundingBox(zMin=val-tol_boun, zMax=val+tol_boun,
                                                              xMin=x_min_search-tol_boun, xMax=x_max_search+tol_boun,
                                                              yMin=y_min_search-tol_boun, yMax=y_max_search+tol_boun)

            if faces_found and len(faces_found) > 0:
                part_obj.Set(faces=faces_found, name=set_name)
                print "  Created Part set '{}' ({} faces)".format(set_name, len(faces_found))
            else:
                print "  WARNING: No faces for Part set '{}'. Creating empty.".format(set_name)
                part_obj.Set(faces=part_module.FaceArray(faces=[]), name=set_name)


def assign_mesh_controls_and_generate(part_obj, mesh_size_global, bridge_cellset_flag, mesh_md_region_flag, mesh_cutout_region_flag):
    print "\nAssigning mesh controls and generating mesh for '{}'...".format(part_obj.name)
    if not part_obj.cells: raise ValueError("Part '{}' has no cells to mesh.".format(part_obj.name))

    meshed_region_cells_list = []
    # Note: bridge cells are initially in 'FE' set before reclassification
    if 'FE' in part_obj.sets and part_obj.sets['FE'].cells:
        meshed_region_cells_list.extend(part_obj.sets['FE'].cells)
        print "  FE region ({} cells) marked for meshing.".format(len(part_obj.sets['FE'].cells))

    if 'PAD' in part_obj.sets and part_obj.sets['PAD'].cells:
        meshed_region_cells_list.extend(part_obj.sets['PAD'].cells)
        print "  PAD region ({} cells) marked for meshing.".format(len(part_obj.sets['PAD'].cells))

    if mesh_md_region_flag:
        if 'MD' in part_obj.sets and part_obj.sets['MD'].cells:
            meshed_region_cells_list.extend(part_obj.sets['MD'].cells)
            print "  MD region ({} cells) marked for meshing (mesh_md_region=True).".format(len(part_obj.sets['MD'].cells))

    if mesh_cutout_region_flag:
        if 'EMPTY' in part_obj.sets and part_obj.sets['EMPTY'].cells:
            meshed_region_cells_list.extend(part_obj.sets['EMPTY'].cells)
            print "  EMPTY region ({} cells) marked for meshing (mesh_cutout_region=True).".format(len(part_obj.sets['EMPTY'].cells))

    if not meshed_region_cells_list:
        print "  WARNING: No cells found in designated regions to mesh. Skipping mesh generation."
        return

    unique_cells_to_mesh = []
    seen_indices = set()
    for cell in meshed_region_cells_list:
        if cell.index not in seen_indices:
            unique_cells_to_mesh.append(cell)
            seen_indices.add(cell.index)

    if not unique_cells_to_mesh:
        print "  WARNING: No unique cells found in designated regions to mesh after filtering. Skipping mesh generation."
        return

    all_meshed_cells_region = part_module.CellArray(cells=unique_cells_to_mesh)

    elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)

    part_obj.setElementType(regions=(all_meshed_cells_region,), elemTypes=(elemType1, elemType2, elemType3))
    print "  Assigned C3D8/C3D6/C3D4 element types to meshed regions."
    try:
        part_obj.setMeshControls(regions=all_meshed_cells_region, technique=SWEEP)
        print "  Assigned SWEEP to meshed region cells."
    except AbaqusException as e_sweep: print "  WARNING: Failed to apply SWEEP mesh control to some/all meshed cells: {}. Using default.".format(e_sweep)

    part_obj.seedPart(size=mesh_size_global, deviationFactor=0.1, minSizeFactor=0.1)
    print "  Seeded entire part globally with size={:.3f}".format(mesh_size_global)

    print "  Generating mesh for designated regions ({} cells)...".format(len(all_meshed_cells_region))
    part_obj.generateMesh(regions=(all_meshed_cells_region,))
    if not part_obj.elements: print "  WARNING: Mesh generation resulted in 0 elements for '{}'.".format(part_obj.name)
    else: print "  Mesh generated for designated regions: {} nodes, {} elements.".format(len(part_obj.nodes), len(part_obj.elements))


def create_iface_nodes(part_obj):
    print "\nCreating 'iface' node set for '{}'...".format(part_obj.name)
    if 'iface' in part_obj.sets: del part_obj.sets['iface']

    if not part_obj.elements:
        print "  Part '{}' has no meshed elements at all. Cannot create 'iface' set.".format(part_obj.name)
        part_obj.Set(name='iface', nodes=part_module.NodeArray([])); return

    if 'MD' not in part_obj.sets or not part_obj.sets['MD'].cells:
        print "  'MD' cell set not found or empty for part '{}'. Cannot find 'iface' nodes.".format(part_obj.name)
        part_obj.Set(name='iface', nodes=part_module.NodeArray([])); return

    fe_bd_node_labels = set()
    sets_to_check = ['FE']
    if bridge_cellset and 'BD' in part_obj.sets: # Check if BD set was actually created
        sets_to_check.append('BD')

    for set_name in sets_to_check:
        if set_name in part_obj.sets:
            # Handle both cell sets and element sets
            if part_obj.sets[set_name].cells:
                cell_obj_array = part_obj.sets[set_name].cells
                for i in range(len(cell_obj_array)):
                    cell = cell_obj_array[i]
                    elements_in_cell = cell.getElements()
                    for elem in elements_in_cell:
                        for node_obj in elem.getNodes():
                            fe_bd_node_labels.add(node_obj.label)
            elif part_obj.sets[set_name].elements:
                elements_in_set = part_obj.sets[set_name].elements
                for elem in elements_in_set:
                    for node_obj in elem.getNodes():
                        fe_bd_node_labels.add(node_obj.label)


    if not fe_bd_node_labels:
        print "  No nodes found in meshed FE or BD regions. 'iface' will be empty."
        part_obj.Set(name='iface', nodes=part_module.NodeArray([])); return

    md_boundary_node_labels = set()
    md_cell_array = part_obj.sets['MD'].cells
    for md_cell_index in range(len(md_cell_array)):
        current_md_cell = md_cell_array[md_cell_index]

        md_cell_vertices_indices = current_md_cell.getVertices()
        for vert_idx in md_cell_vertices_indices:
            vertex_obj = part_obj.vertices[vert_idx]
            nodes_on_vert = vertex_obj.getNodes()
            if nodes_on_vert: md_boundary_node_labels.update(n.label for n in nodes_on_vert)

        md_cell_edges_indices = current_md_cell.getEdges()
        for edge_idx in md_cell_edges_indices:
            edge_obj = part_obj.edges[edge_idx]
            nodes_on_edge = edge_obj.getNodes()
            if nodes_on_edge: md_boundary_node_labels.update(n.label for n in nodes_on_edge)

        md_cell_faces_indices = current_md_cell.getFaces()
        for face_idx in md_cell_faces_indices:
            face_obj = part_obj.faces[face_idx]
            nodes_on_face = face_obj.getNodes()
            if nodes_on_face: md_boundary_node_labels.update(n.label for n in nodes_on_face)

    if not md_boundary_node_labels:
        print "  No geometric boundary nodes found for MD cells. 'iface' will be empty."
        part_obj.Set(name='iface', nodes=part_module.NodeArray([])); return

    iface_labels = list(fe_bd_node_labels.intersection(md_boundary_node_labels))

    if iface_labels:
        print "  Identified {} candidate interface node labels from intersection.".format(len(iface_labels))
        try:
            iface_nodes_seq = part_obj.nodes.sequenceFromLabels(labels=iface_labels)
            if not iface_nodes_seq:
                print "  Resulting node sequence for 'iface' is empty. Creating an empty set."
                part_obj.Set(name='iface', nodes=part_module.NodeArray([]))
                return

            part_obj.Set(name='iface', nodes=iface_nodes_seq)
            print "  Successfully created 'iface' node set with {} nodes.".format(len(iface_nodes_seq))
        except AbaqusException as e_iface:
            print "  Error creating 'iface' set from labels: {}. Creating an empty set.".format(e_iface)
            part_obj.Set(name='iface', nodes=part_module.NodeArray([]))
    else:
        print "  No common interface nodes found between (FE/BD element nodes) and (MD cell boundary nodes). Empty 'iface' set."
        part_obj.Set(name='iface', nodes=part_module.NodeArray([]))


def del_empty_models():
    print "\nCleaning up Abaqus session..."
    models_to_delete = [m for m in mdb.models.keys() if not mdb.models[m].parts]
    for model_key in models_to_delete: del mdb.models[model_key]; print "  Deleted empty model '{}'.".format(model_key)
    if not models_to_delete: print "  No empty models found."

# ==============================================================================
# Main mesh generation execution logic
# ==============================================================================
def run_meshing_process():
    print "="*70; print " LAMMPS Box Meshing Script ".center(70, "="); print "="*70
    myModel = None; myPart = None

    script_dir = os.getcwd()
    try:
        data_file, _ = find_lammps_data_file(script_dir, lammps_file)
        data_file = os.path.abspath(data_file)

        output_file_abs = os.path.abspath(output_file)
        output_dir = os.path.dirname(output_file_abs)
        base_name_for_objects = os.path.basename(output_file_abs)

        if output_dir and not os.path.exists(output_dir):
            print "Creating output directory: {}".format(output_dir)
            os.makedirs(output_dir)
        
        if output_dir:
            print "Changing Abaqus working directory to: {}".format(output_dir)
            os.chdir(output_dir)

        model_name, part_name, job_name, abq_step_name = base_name_for_objects, base_name_for_objects, base_name_for_objects, base_name_for_objects
        
        if model_name not in mdb.models: myModel = mdb.Model(name=model_name)
        else:
            myModel = mdb.models[model_name]
            if part_name in myModel.parts: del myModel.parts[part_name]
            if job_name in mdb.jobs: del mdb.jobs[job_name]

        lammps_box_orig = parse_lammps_data(data_file)
        overall_part_bounds, physical_partition_coords, md_interface_bounds, boun_set_bounds = calculate_dimensions(
            lammps_box_orig, domain_min_coords, domain_max_coords, bridge_thickness, tol, full_elem_bridge,
            cut_width, cut_shift, notch_width, notch_normal, notch_plane, notch_shift, notch_depth,
            have_padding, default_mesh_size
        )
        myPart = create_partitioned_part(myModel, part_name, overall_part_bounds, physical_partition_coords, tol)

        if myPart:
            create_part_geometry_sets(myPart, md_interface_bounds, physical_partition_coords, tol, bridge_cellset, lammps_box_orig,
                                      overall_part_bounds,
                                      cut_width, cut_shift, notch_width, notch_normal, notch_plane, notch_shift, notch_depth,
                                      have_padding, boun_set_bounds)
            create_part_outer_boundary_sets(myPart, overall_part_bounds, boun_set_bounds, tol)

            assign_mesh_controls_and_generate(myPart, default_mesh_size, bridge_cellset, mesh_md_region, mesh_cutout_region)
            
            if bridge_cellset and myPart.elements:
                create_bd_element_set(myPart, lammps_box_orig, tol)

            if have_iface:
                create_iface_nodes(myPart)
            else:
                print "\nSkipping 'iface' node set creation as per 'have_iface = False'."

            if myPart.elements :
                myAssembly = myModel.rootAssembly
                instance_name = part_name + '-1'
                if instance_name in myAssembly.instances: del myAssembly.instances[instance_name]
                myAssembly.regenerate()
                inst = myAssembly.Instance(name=instance_name, part=myPart, dependent=ON)
                myAssembly.regenerate(); print "\nAssembly instance '{}' created.".format(instance_name)

                print "\nCreating Step: '{}'...".format(abq_step_name)
                if abq_step_name in myModel.steps: del myModel.steps[abq_step_name]
                myModel.StaticStep(name=abq_step_name, previous='Initial'); print "  Step '{}' created.".format(abq_step_name)

                print "\nCreating Job: '{}'...".format(job_name)
                mdb.Job(name=job_name, model=model_name, type=ANALYSIS, nodalOutputPrecision=FULL, numCpus=1)
                print "  Job '{}' created.".format(job_name)
                mdb.jobs[job_name].writeInput(consistencyChecking=OFF)
                print "  Input file '{}.inp' written to {}.".format(job_name, os.getcwd())

                if myModel:
                    cae_file_name = base_name_for_objects + ".cae"
                    mdb.saveAs(pathName=cae_file_name)
                    print "  CAE file saved as '{}' in {}.".format(cae_file_name, os.getcwd())

                print "\nUpdating viewport..."
                if session.viewports:
                    vp_keys = session.viewports.keys()
                    if vp_keys:
                        vp = session.viewports[vp_keys[0]]
                        vp.setValues(displayedObject=myAssembly); vp.assemblyDisplay.setValues(mesh=ON, loads=ON, bcs=ON)
                        if 'Iso' in session.views: vp.view.setValues(session.views['Iso'])
                        vp.view.fitView(); print "  Viewport updated."
                    else: print "  No active viewport found to update."
            else: print "\nSkipping Assembly/Step/BCs/Job: No elements were generated or part is invalid."
        else:
            print "\nERROR: Part creation failed. Skipping subsequent steps."

        del_empty_models()

    except IOError as e:
        print "IOError during mesh creation: {}".format(e)
        traceback.print_exc()
    except ValueError as e_val:
        print "ValueError during setup: {}".format(e_val)
        traceback.print_exc()
    except AbaqusException as e_abq:
        print "AbaqusException: {}".format(e_abq)
        traceback.print_exc()
    except Exception as e_unexp:
        print "Unexpected error during mesh creation: {}".format(e_unexp)
        traceback.print_exc()
    print "\n"


# ==============================================================================
# ==============================================================================
# Main execution
# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
    run_meshing_process()

#    # Clean up folder
#    print "\nCleaning up temporary files..."
#    for f in os.listdir('.'):
#        if f.endswith(('.pyc', '.jnl')) or 'abaqus.rpy' in f:
#            try: os.remove(f)
#            except OSError: pass