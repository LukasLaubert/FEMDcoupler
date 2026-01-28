import sys
import os
import json
import re
from collections import OrderedDict
import glob

# ==============================================================================
# Logger Class
# ==============================================================================
class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()
        self.stdout.flush()

# --- Data Structures for LAMMPS File ---
class AtomStyleInfo:
    def __init__(self, style_name, id_idx, mol_id_idx, type_idx, q_idx, x_idx, y_idx, z_idx):
        self.style_name = style_name
        self.id_idx = id_idx
        self.mol_id_idx = mol_id_idx
        self.type_idx = type_idx
        self.q_idx = q_idx
        self.x_idx = x_idx
        self.y_idx = y_idx
        self.z_idx = z_idx

class LammpsData:
    def __init__(self):
        self.initial_comment = ""
        self.header_counts = OrderedDict()
        self.box_bounds_lines = []
        self.atom_style_info = None
        self.section_headers = {}
        self.atoms = []
        self.bonds = []
        self.velocities = []
        self.masses = []
        self.other_sections = OrderedDict()

# --- LAMMPS File Parsing ---
def _determine_atom_style_info(atom_section_header_line):
    parts = atom_section_header_line.split('#', 1)
    style = parts[1].strip().lower() if len(parts) > 1 else "atomic"
    
    mol_id_idx = 2 if 'molecular' in style or 'full' in style else None
    type_idx = 3 if mol_id_idx is not None else 2
    current_idx = type_idx + 1

    q_idx = current_idx if 'charge' in style or 'full' in style else None
    if q_idx is not None: current_idx += 1

    return AtomStyleInfo(style, 1, mol_id_idx, type_idx, q_idx, current_idx, current_idx + 1, current_idx + 2)

def _parse_atom_line(line_str, style_info):
    parts = line_str.split()
    try:
        atom = {
            'id': int(parts[style_info.id_idx-1]),
            'type': int(parts[style_info.type_idx-1]),
            'x': float(parts[style_info.x_idx-1]),
            'y': float(parts[style_info.y_idx-1]),
            'z': float(parts[style_info.z_idx-1]),
            'original_parts': parts
        }
        if style_info.mol_id_idx is not None: atom['mol_id'] = int(parts[style_info.mol_id_idx-1])
        if style_info.q_idx is not None: atom['q'] = float(parts[style_info.q_idx-1])
        return atom
    except (ValueError, IndexError):
        return None

def _parse_lammps_data(filepath):
    data = LammpsData()
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    data.initial_comment = lines[0].strip()
    
    header_keywords = ["atoms", "bonds", "angles", "dihedrals", "impropers", "atom types", "bond types", "angle types", "dihedral types", "improper types"]
    known_sections = ["Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", "PairIJ Coeffs"]
    modifiable_sections = ["Atoms", "Bonds", "Velocities", "Masses"]
    
    current_section_name = None
    line_iter = iter(lines[1:])
    for line in line_iter:
        ln = line.strip()
        if not ln: continue
        words = ln.split()

        try:
            potential_kw = " ".join(words[1:])
            if potential_kw in header_keywords:
                data.header_counts[potential_kw] = int(words[0])
                current_section_name = None
                continue
        except (ValueError, IndexError):
            pass

        if any(s in ln for s in ["xlo xhi", "ylo yhi", "zlo zhi"]):
            data.box_bounds_lines.append(ln)
            current_section_name = None
            continue

        section_header_candidate = ln.split('#')[0].strip()
        if section_header_candidate in known_sections:
            current_section_name = section_header_candidate
            data.section_headers[current_section_name] = ln
            if current_section_name not in modifiable_sections:
                data.other_sections[current_section_name] = []
            continue

        if current_section_name:
            if current_section_name == "Atoms":
                if data.atom_style_info is None:
                    data.atom_style_info = _determine_atom_style_info(data.section_headers["Atoms"])
                atom_obj = _parse_atom_line(ln, data.atom_style_info)
                if atom_obj: data.atoms.append(atom_obj)
            elif current_section_name == "Bonds":
                data.bonds.append([int(p) for p in words[:4]])
            elif current_section_name == "Velocities":
                data.velocities.append([int(words[0])] + [float(p) for p in words[1:4]])
            elif current_section_name == "Masses":
                data.masses.append(words[:2])
            else:
                data.other_sections[current_section_name].append(ln)
    return data

# --- FE Mesh Data Extraction ---
def _parse_inp_file(filepath):
    """Robustly parse .inp file for nodes, elements, and element sets, including generate."""
    nodes, elements, elsets = {}, {}, {}
    in_node_section, in_element_section, in_elset_section = False, False, False
    current_elset_name = None
    elset_generate_mode = False

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('**'):
                continue

            if stripped.startswith('*'):
                in_node_section, in_element_section, in_elset_section = False, False, False
                lowered_line = stripped.lower()
                
                if lowered_line.startswith('*node'):
                    in_node_section = True
                elif lowered_line.startswith('*element'):
                    in_element_section = True
                elif lowered_line.startswith('*elset'):
                    in_elset_section = True
                    # Check for generate parameter
                    elset_generate_mode = 'generate' in lowered_line
                    
                    match = re.search(r'elset=([^,]+)', stripped, re.IGNORECASE)
                    if match:
                        current_elset_name = match.group(1).strip()
                        if current_elset_name not in elsets:
                            elsets[current_elset_name] = set()
                    else:
                        current_elset_name = None
                continue
            
            if in_node_section:
                parts = [p.strip() for p in stripped.split(',')]
                try:
                    nodes[int(parts[0])] = [float(p) for p in parts[1:4]]
                except ValueError:
                    pass
            elif in_element_section:
                parts = [p.strip() for p in stripped.split(',')]
                try:
                    elem_id = int(parts[0])
                    elem_nodes = [int(p) for p in parts[1:] if p]
                    elem_type = 'C3D8' if len(elem_nodes) == 8 else 'C3D4' if len(elem_nodes) == 4 else 'C3D6' if len(elem_nodes) == 6 else None
                    elements[elem_id] = {'type': elem_type, 'nodes': elem_nodes}
                except ValueError:
                    pass
            elif in_elset_section and current_elset_name:
                parts = [p.strip() for p in stripped.split(',') if p.strip()]
                if elset_generate_mode:
                    # Syntax: start, end, increment
                    if len(parts) == 3:
                        try:
                            start, end, inc = int(parts[0]), int(parts[1]), int(parts[2])
                            elsets[current_elset_name].update(range(start, end + 1, inc))
                        except ValueError:
                            pass
                else:
                    # Standard list of IDs
                    elsets[current_elset_name].update(int(p) for p in parts if p.isdigit())
                
    return nodes, elements, elsets

def _extract_fe_mesh_data(all_nodes, all_elements, all_sets, cellset_names):
    EDGE_CONNECTIVITY = {
        'C3D8': [(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4), (0,4), (1,5), (2,6), (3,7)],
        'C3D6': [(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (0,3), (1,4), (2,5)],
        'C3D4': [(0,1), (1,2), (2,0), (0,3), (1,3), (2,3)]
    }
    
    relevant_elems = set()
    for name in cellset_names:
        if name in all_sets:
            relevant_elems.update(all_sets[name])

    fe_nodes, fe_edges = set(), set()
    for elem_id in relevant_elems:
        if elem_id not in all_elements: continue
        element = all_elements[elem_id]
        fe_nodes.update(element['nodes'])
        connectivity = EDGE_CONNECTIVITY.get(element['type'])
        if connectivity:
            for i, j in connectivity:
                try:
                    n1, n2 = element['nodes'][i], element['nodes'][j]
                    fe_edges.add(tuple(sorted((n1, n2))))
                except IndexError:
                    pass
    
    return sorted(list(fe_nodes)), sorted(list(fe_edges))

# --- LAMMPS File Writing ---
def _create_fe_atom_line(id, type_id, coords, style_info, template_parts):
    num_cols = len(template_parts)
    parts = ['0'] * num_cols
    
    parts[style_info.id_idx - 1] = str(id)
    parts[style_info.type_idx - 1] = str(type_id)
    if style_info.mol_id_idx is not None: parts[style_info.mol_id_idx - 1] = "0"
    if style_info.q_idx is not None: parts[style_info.q_idx - 1] = "0.0"
    parts[style_info.x_idx - 1] = "{:.15g}".format(coords[0])
    parts[style_info.y_idx - 1] = "{:.15g}".format(coords[1])
    parts[style_info.z_idx - 1] = "{:.15g}".format(coords[2])
    
    return " ".join(parts)

def _write_modified_lammps_file(output_filepath, data, new_atoms, new_bonds, new_velocities, new_atom_type, new_bond_type):
    with open(output_filepath, "w") as f:
        f.write(data.initial_comment + "\n\n")
        
        data.header_counts["atoms"] = data.header_counts.get("atoms", 0) + len(new_atoms)
        data.header_counts["bonds"] = data.header_counts.get("bonds", 0) + len(new_bonds)
        data.header_counts["atom types"] = max(data.header_counts.get("atom types", 0), new_atom_type)
        data.header_counts["bond types"] = max(data.header_counts.get("bond types", 0), new_bond_type)
        for key, val in data.header_counts.items(): f.write("{} {}\n".format(val, key))
        
        f.write("\n" + "\n".join(data.box_bounds_lines) + "\n\n")

        section_order = ["Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"]
        
        for section in section_order:
            if section in data.section_headers or (section in ["Bonds", "Velocities"] and (data.bonds or data.velocities or new_bonds or new_velocities)):
                if section not in data.section_headers: data.section_headers[section] = section
                header_line = data.section_headers[section]
                
                has_content = False
                if section == "Masses" or section == "Atoms": has_content = True
                elif section == "Bonds": has_content = bool(data.bonds or new_bonds)
                elif section == "Velocities": has_content = bool(data.velocities or new_velocities)
                elif section in data.other_sections: has_content = bool(data.other_sections[section])

                if has_content:
                    f.write(header_line + "\n\n")
                    if section == "Masses":
                        for mass_entry in data.masses: f.write(" ".join(mass_entry) + "\n")
                        f.write("{} 1.0e-50\n".format(new_atom_type))
                    elif section == "Atoms":
                        for atom in sorted(data.atoms, key=lambda x: x['id']): f.write(" ".join(atom['original_parts']) + "\n")
                        if new_atoms: f.write("\n".join(new_atoms) + "\n")
                    elif section == "Velocities":
                        for v in sorted(data.velocities, key=lambda x: x[0]): f.write("{} {:.8e} {:.8e} {:.8e}\n".format(*v))
                        if new_velocities: f.write("\n".join(new_velocities) + "\n")
                    elif section == "Bonds":
                        for bond in sorted(data.bonds, key=lambda x: x[0]): f.write(" ".join(map(str, bond)) + "\n")
                        if new_bonds: f.write("\n".join(new_bonds) + "\n")
                    elif section in data.other_sections:
                        f.write("\n".join(data.other_sections[section]) + "\n")
                    f.write("\n")

# --- Main Function ---
def run_femd_merger():
    print("="*70 + "\n" + " FEMD Grid/Particle Merger Initialized ".center(70) + "\n" + "="*70)
    try:
        if len(sys.argv) > 1:
            with open(sys.argv[1], 'r') as f: params = json.load(f)
            
            datafile_path = os.path.splitext(params['output_file'])[0] + '.data'
            inpfile_path = os.path.splitext(params['output_file'])[0] + '.inp'
            cellset_names = ["FE", "BD"] if params['bridge_cellset'] else ["FE"]
        else:
            print("INFO: No config file specified, searching for default files in current directory.")
            data_files = glob.glob('*.data')
            if not data_files: raise IOError("No .data file found in current directory.")
            datafile_path = data_files[0]
            inp_files = glob.glob('*.inp')
            if not inp_files: raise IOError("No .inp file found in current directory.")
            inpfile_path = inp_files[0]
            cellset_names = ["FE", "BD"] # Assume both if no config is given

        print("Using MD data file: {}".format(datafile_path))
        print("Using FE inp file:  {}".format(inpfile_path))
        
        all_nodes, all_elements, all_sets = _parse_inp_file(inpfile_path)
        fe_nodes, fe_edges = _extract_fe_mesh_data(all_nodes, all_elements, all_sets, cellset_names)
        print("Extracted {} nodes and {} edges from FE grid.".format(len(fe_nodes), len(fe_edges)))

        lammps_data = _parse_lammps_data(datafile_path)
        print("Parsed {} atoms and {} bonds from LAMMPS data.".format(len(lammps_data.atoms), len(lammps_data.bonds)))
        
        if not lammps_data.atoms: raise ValueError("LAMMPS data file contains no atoms to use as a format template.")
        template_parts = lammps_data.atoms[0]['original_parts']

        max_atom_id = max(a['id'] for a in lammps_data.atoms) if lammps_data.atoms else 0
        max_bond_id = max(b[0] for b in lammps_data.bonds) if lammps_data.bonds else 0
        max_atom_type = max(int(m[0]) for m in lammps_data.masses) if lammps_data.masses else 0
        max_bond_type = max(b[1] for b in lammps_data.bonds) if lammps_data.bonds else 0
        
        fe_atom_type = max_atom_type + 1
        fe_bond_type = max_bond_type + 1

        new_atoms, new_velocities, new_bonds = [], [], []
        fe_node_id_to_new_atom_id = {}
        
        has_original_velocities = bool(lammps_data.velocities)

        for i, fe_node_id in enumerate(fe_nodes):
            new_atom_id = max_atom_id + i + 1
            fe_node_id_to_new_atom_id[fe_node_id] = new_atom_id
            coords = all_nodes[fe_node_id]
            new_atoms.append(_create_fe_atom_line(new_atom_id, fe_atom_type, coords, lammps_data.atom_style_info, template_parts))
            if has_original_velocities:
                new_velocities.append("{} 0.0 0.0 0.0".format(new_atom_id))

        for i, edge in enumerate(fe_edges):
            new_bond_id = max_bond_id + i + 1
            atom1_id = fe_node_id_to_new_atom_id[edge[0]]
            atom2_id = fe_node_id_to_new_atom_id[edge[1]]
            new_bonds.append("{} {} {} {}".format(new_bond_id, fe_bond_type, atom1_id, atom2_id))

        output_filename = os.path.splitext(os.path.basename(datafile_path))[0] + "_FEMD_merged.data"
        output_path = os.path.join(os.path.dirname(datafile_path) or '.', output_filename)
        
        _write_modified_lammps_file(output_path, lammps_data, new_atoms, new_bonds, new_velocities, fe_atom_type, fe_bond_type)
        
        print("\nSuccessfully merged FE data.")
        print("  - Added {} atoms (type {})".format(len(new_atoms), fe_atom_type))
        print("  - Added {} bonds (type {})".format(len(new_bonds), fe_bond_type))
        print("Output written to: {}".format(output_path))

    except Exception as e:
        print("\nERROR: Process failed.")
        print("Reason: {}".format(e))
    finally:
        print("="*70)

if __name__ == "__main__":
    try:
        config_file_init = None
        if len(sys.argv) > 1 and '.json' in sys.argv[1]:
            config_file_init = sys.argv[1]
        else:
            json_files = glob.glob('*.json')
            if json_files:
                config_file_init = json_files[0]
        
        if config_file_init:
            with open(config_file_init, 'r') as f:
                params_init = json.load(f)
                if params_init.get('generate_logs', False):
                    output_file_init = params_init.get('output_file')
                    if output_file_init:
                        sys.stdout = Tee(output_file_init + "_genFE.log", "a")
    except:
        pass
        
    run_femd_merger()