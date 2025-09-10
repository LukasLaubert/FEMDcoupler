#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
create_periodicity_mesh.py

Processes .inp files to adjust node coordinates for periodicity.
Groups nodes that are approximately aligned along certain directions and
makes their coordinates in the other directions exactly identical.
"""

import os
import re
import sys
from collections import Counter
import json
import glob

def parse_parameters():
    """Parse parameters from a .json file specified by argument or found in the path."""
    config_filename = None
    if len(sys.argv) > 1 and '.json' in sys.argv[1]:
        config_filename = sys.argv[1]
    else:
        json_files = glob.glob('*.json')
        if not json_files:
            print("Error: No .json configuration file found in the current directory.")
            sys.exit(1)
        config_filename = json_files[0]
        print(f"Info: No JSON file given as argument. Using first found: '{config_filename}'")
        
    try:
        with open(config_filename, 'r') as f:
            params = json.load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file '{config_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing JSON file '{config_filename}': {e}")
        sys.exit(1)
    
    return params

def find_inp_file(params):
    """Find .inp files, prioritizing output_file from params, then first available."""
    output_file = params.get('output_file')
    preferred_file = None
    if output_file:
        base_name = os.path.splitext(output_file)[0]
        preferred_file = f"{base_name}.inp"
        if os.path.exists(preferred_file):
            return preferred_file, f"Using file from parameters: {preferred_file}"

    # Fallback to searching the directory
    inp_files = [f for f in os.listdir('.') if f.endswith('.inp')]
    if not inp_files:
        print("Error: No .inp files found in current directory.")
        sys.exit(1)

    selected_file = inp_files[0]
    if preferred_file:
        message = f"File '{preferred_file}' not found. Using first available: {selected_file}"
    else:
        message = f"Using first available .inp file: {selected_file}"

    return selected_file, message

def parse_inp_file(file_path):
    """Parse .inp file to extract node coordinates, structure, and Nsets."""
    nodes_dict = {}
    file_structure = []
    node_formats = {}
    nsets = {}
    
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        in_node_section = False
        in_nset_section = False
        current_nset_name = None
        
        for line in lines:
            file_structure.append(line)
            stripped_line = line.strip()
            
            if stripped_line.startswith('*Node'):
                in_node_section = True
                in_nset_section = False
                continue
            
            if stripped_line.lower().startswith('*nset'):
                in_node_section = False
                in_nset_section = True
                match = re.search(r'nset=([^,]+)', stripped_line, re.IGNORECASE)
                if match:
                    current_nset_name = match.group(1).strip()
                    if current_nset_name.startswith('BOUN_'):
                        nsets[current_nset_name] = set()
                else:
                    current_nset_name = None
                continue

            if stripped_line.startswith('*'):
                in_node_section = False
                in_nset_section = False
                current_nset_name = None
                continue
            
            if in_node_section and stripped_line:
                parts = stripped_line.split(',')
                if len(parts) >= 4:
                    try:
                        node_id = int(parts[0].strip())
                        x, y, z = float(parts[1].strip()), float(parts[2].strip()), float(parts[3].strip())
                        nodes_dict[node_id] = [x, y, z]
                        
                        node_line_match = re.match(r'^(\s*)(\d+),(\s*)([^,]+),(\s*)([^,]+),(\s*)([^,]+)', line)
                        if node_line_match:
                            node_formats[node_id] = {
                                'leading_spaces': node_line_match.group(1), 'after_node_id': node_line_match.group(3),
                                'after_x': node_line_match.group(5), 'after_y': node_line_match.group(7),
                                'x_original': parts[1].strip(), 'y_original': parts[2].strip(),
                                'z_original': parts[3].strip(), 'full_line': line.rstrip('\n')
                            }
                    except ValueError:
                        pass
            
            if in_nset_section and current_nset_name and current_nset_name.startswith('BOUN_') and stripped_line:
                node_ids = [int(n.strip()) for n in stripped_line.split(',') if n.strip().isdigit()]
                nsets[current_nset_name].update(node_ids)
                        
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing .inp file: {e}")
        sys.exit(1)
    
    return nodes_dict, file_structure, node_formats, nsets

def pre_correct_boun_nodes(nodes_dict, nsets):
    """Adjusts nodes in BOUN sets to be perfectly aligned on their respective planes."""
    print("\nPre-correcting BOUN node sets for perfect alignment:")
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    for nset_name, node_ids in nsets.items():
        if not node_ids or not nset_name.startswith('BOUN_'):
            continue

        axis_char = nset_name[5]
        if axis_char not in axis_map:
            continue
        
        axis_index = axis_map[axis_char]
        
        # Voting: find the most common coordinate
        coords = [nodes_dict[nid][axis_index] for nid in node_ids if nid in nodes_dict]
        if not coords:
            continue
        
        voted_coord = Counter(coords).most_common(1)[0][0]
        
        # Adjust nodes that deviate
        adjustment_count = 0
        for nid in node_ids:
            if nid in nodes_dict and abs(nodes_dict[nid][axis_index] - voted_coord) > 1e-12:
                nodes_dict[nid][axis_index] = voted_coord
                adjustment_count += 1
        
        if adjustment_count > 0:
            print(f"  Adjusted {adjustment_count} nodes in {nset_name} to {axis_char}={voted_coord}")

def get_processing_directions(params):
    """Determine which directions need processing based on parameters."""
    if params.get('circular_burrito_radius') is not None:
        # For circular systems, only process the specified pbc_direction
        return [params['pbc_direction']]
    else:
        # For rectangular systems, process directions with None in domain bounds
        directions = []
        for i in range(3):
            if (params['domain_min_coords'][i] is None or 
                params['domain_max_coords'][i] is None):
                directions.append(i + 1)  # 1-indexed
        return directions

def calculate_typical_nodal_distance(nodes_dict, direction):
    """Calculate typical nodal distance in the processing direction."""
    processing_index = direction - 1
    processing_coords = [coords[processing_index] for coords in nodes_dict.values()]
    sorted_coords = sorted(processing_coords)
    
    differences = []
    for i in range(1, len(sorted_coords)):
        diff = sorted_coords[i] - sorted_coords[i-1]
        if diff > 1e-12:
            differences.append(diff)
    
    if not differences:
        return 1.0
    
    differences.sort()
    median_index = len(differences) // 2
    if len(differences) % 2 == 1:
        return differences[median_index]
    else:
        return (differences[median_index - 1] + differences[median_index]) / 2

def get_boundary_node_count(nodes_dict, direction, nsets):
    """Get number of unique boundary node positions for validation using BOUN_...lo Nset."""
    # Determine non-processing indices
    if direction == 1:
        non_processing_indices = [1, 2]  # y, z
    elif direction == 2:
        non_processing_indices = [0, 2]  # x, z
    else:
        non_processing_indices = [0, 1]  # x, y
        
    axis_char = ['x', 'y', 'z'][direction-1]
    nset_name = f'BOUN_{axis_char}lo'
    
    if nset_name not in nsets:
        print(f"Error: Required node set '{nset_name}' not found in .inp file.")
        sys.exit(1)

    boundary_node_ids = nsets[nset_name]
    
    # Find unique boundary positions from the specified Nset
    boundary_positions = set()
    for node_id in boundary_node_ids:
        if node_id in nodes_dict:
            coords = nodes_dict[node_id]
            boundary_pos = tuple(coords[i] for i in non_processing_indices)
            boundary_positions.add(boundary_pos)
    
    return len(boundary_positions), boundary_positions

def _create_groups_with_tolerance(coord_combinations, tolerance):
    """Create groups of nodes within tolerance."""
    groups = {}
    used_coords = set()
    
    for base_coord, base_nodes in coord_combinations.items():
        if base_coord in used_coords:
            continue
            
        # Find all coordinates within tolerance
        group_coords = [base_coord]
        group_nodes = base_nodes[:]
        
        for other_coord, other_nodes in coord_combinations.items():
            if other_coord != base_coord and other_coord not in used_coords:
                # Check if all coordinates are within tolerance
                if all(abs(a - b) <= tolerance for a, b in zip(base_coord, other_coord)):
                    group_coords.append(other_coord)
                    group_nodes.extend(other_nodes)
                    used_coords.add(other_coord)
        
        # Use the base coordinate as key
        groups[base_coord] = group_nodes
        used_coords.add(base_coord)
    
    return groups

def group_nodes_by_direction(nodes_dict, direction, tol, nsets):
    """Group nodes with adaptive tolerance search."""
    target_group_count, _ = get_boundary_node_count(nodes_dict, direction, nsets)
    
    # Determine coordinate indices
    if direction == 1:
        coord_indices = [1, 2]
    elif direction == 2:
        coord_indices = [0, 2]
    else:
        coord_indices = [0, 1]
    
    # Collect coordinate combinations
    coord_combinations = {}
    for node_id, coords in nodes_dict.items():
        group_coords = tuple(coords[i] for i in coord_indices)
        if group_coords not in coord_combinations:
            coord_combinations[group_coords] = []
        coord_combinations[group_coords].append(node_id)
    
    # Adaptive tolerance search
    tolerance_factor = 1.0
    adjustment_factor = 10.0
    current_tolerance = tol * tolerance_factor
    
    best_groups, best_tolerance, best_diff = None, None, float('inf')
    max_iterations, iteration = 30, 0
    history = []
    
    while iteration < max_iterations:
        iteration += 1
        
        groups = _create_groups_with_tolerance(coord_combinations, current_tolerance)
        group_count = len(groups)
        diff = abs(group_count - target_group_count)
        
        print(f"  Iter {iteration}: tol={current_tolerance:.2e}, groups={group_count}, target={target_group_count}, diff={diff}")
        
        if group_count == target_group_count:
            best_groups, best_tolerance = groups, current_tolerance
            break
        
        if diff < best_diff or (diff == best_diff and current_tolerance > best_tolerance):
            best_diff, best_groups, best_tolerance = diff, groups, current_tolerance
        
        if diff <= 1 and iteration > 10:
            best_groups, best_tolerance = groups, current_tolerance
            break
        
        history.append({'group_count': group_count})
        
        if len(history) >= 3:
            recent = history[-3:]
            if ((recent[-2]['group_count'] > target_group_count and recent[-1]['group_count'] < target_group_count) or
                (recent[-2]['group_count'] < target_group_count and recent[-1]['group_count'] > target_group_count)):
                adjustment_factor *= 0.5
                print(f"  Oscillation detected, reduced adjustment factor to {adjustment_factor:.2f}")
        
        if group_count > target_group_count:
            tolerance_factor *= adjustment_factor
        else:
            tolerance_factor /= adjustment_factor
        
        current_tolerance = tol * tolerance_factor
    
    return best_groups, best_tolerance

def adjust_node_coordinates(nodes_dict, groups, direction, tolerance):
    """Adjust node coordinates for periodicity."""
    updated_nodes = nodes_dict.copy()
    adjustment_count = 0
    
    # Determine non-processing indices
    if direction == 1:
        non_processing_indices = [1, 2]
    elif direction == 2:
        non_processing_indices = [0, 2]
    else:
        non_processing_indices = [0, 1]
    
    for base_coord, node_ids in groups.items():
        if len(node_ids) <= 1:
            continue
        
        # Find most frequent values for non-processing coordinates
        coord_values = {i: [] for i in non_processing_indices}
        for node_id in node_ids:
            for idx in non_processing_indices:
                coord_values[idx].append(updated_nodes[node_id][idx])
        
        selected_values = {}
        for idx in non_processing_indices:
            counter = Counter(coord_values[idx])
            selected_values[idx] = counter.most_common(1)[0][0]
        
        # Apply adjustments
        for node_id in node_ids:
            adjusted = False
            for idx in non_processing_indices:
                if abs(updated_nodes[node_id][idx] - selected_values[idx]) > 1e-12:
                    updated_nodes[node_id][idx] = selected_values[idx]
                    adjusted = True
            if adjusted:
                adjustment_count += 1
    
    return updated_nodes, adjustment_count

def validate_periodicity(nodes_dict, direction, tolerance, nsets):
    """Validate that periodicity is achieved."""
    # Determine non-processing indices
    if direction == 1:
        non_processing_indices = [1, 2]
    elif direction == 2:
        non_processing_indices = [0, 2]
    else:
        non_processing_indices = [0, 1]
    
    # Group nodes by non-processing coordinates
    position_groups = {}
    for node_id, coords in nodes_dict.items():
        pos_key = tuple(coords[i] for i in non_processing_indices)
        if pos_key not in position_groups:
            position_groups[pos_key] = []
        position_groups[pos_key].append(node_id)
    
    # Check if all coordinates are identical within groups
    all_identical = True
    for pos_key, node_ids in position_groups.items():
        if len(node_ids) <= 1:
            continue
        
        base_coords = nodes_dict[node_ids[0]]
        for node_id in node_ids[1:]:
            for idx in non_processing_indices:
                if abs(nodes_dict[node_id][idx] - base_coords[idx]) > 1e-12:
                    all_identical = False
                    break
            if not all_identical:
                break
        if not all_identical:
            break
    
    # Get boundary node count for comparison
    boundary_count, _ = get_boundary_node_count(nodes_dict, direction, nsets)
    
    is_valid = (len(position_groups) == boundary_count) and all_identical
    
    status = "PASSED" if is_valid else "FAILED"
    message = f"  Validation: {status} (positions: {len(position_groups)}, boundary: {boundary_count})"
    
    return is_valid, message

def write_output_file(file_structure, nodes_dict, node_formats, output_file):
    """Write the output file with updated node coordinates."""
    try:
        with open(output_file, 'w') as f:
            in_node_section = False
            
            for line in file_structure:
                if line.strip().startswith('*Node'):
                    in_node_section = True
                    f.write(line)
                    continue
                
                if in_node_section and line.strip().startswith('*'):
                    in_node_section = False
                    f.write(line)
                    continue
                
                if in_node_section and line.strip():
                    parts = line.strip().split(',')
                    if len(parts) >= 4:
                        try:
                            node_id = int(parts[0].strip())
                            if node_id in nodes_dict and node_id in node_formats:
                                # Get updated coordinates
                                x, y, z = nodes_dict[node_id]
                                
                                # Format coordinates to match original precision and format
                                fmt = node_formats[node_id]
                                
                                # Use the same formatting logic as the original
                                def format_coord(value, original):
                                    if value == int(value):
                                        return f"{int(value)}."
                                    else:
                                        formatted = f"{value:.10f}".rstrip('0').rstrip('.')
                                        return formatted if formatted else "0"
                                
                                x_str = format_coord(x, fmt['x_original'])
                                y_str = format_coord(y, fmt['y_original'])
                                z_str = format_coord(z, fmt['z_original'])
                                
                                # Build the line with proper spacing and commas
                                new_line = (fmt['leading_spaces'] + f"{node_id}," + fmt['after_node_id'] + x_str + "," +
                                          fmt['after_x'] + y_str + "," + fmt['after_y'] + z_str + '\n')
                                
                                f.write(new_line)
                            else:
                                f.write(line)
                        except ValueError:
                            f.write(line)
                    else:
                        f.write(line)
                else:
                    f.write(line)
                    
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

def main():
    """Main function to execute the periodicity mesh creation process."""
    print("="*60)
    print("PERIODICITY MESH CREATION")
    print("="*60)
    
    # Parse parameters
    params = parse_parameters()
    
    # Check if dPBC_correction is enabled
    if not params.get('dPBC_correction', False):
        print("\nInfo: dPBC_correction is False.")
        print("No changes to ensure dPBC compatibility will be applied.")
        print("Exiting with code 0.")
        sys.exit(0)
    
    # Find and parse input file
    inp_file, message = find_inp_file(params)
    print(f"\nInput file: {message}")
    
    nodes_dict, file_structure, node_formats, nsets = parse_inp_file(inp_file)
    print(f"Nodes found: {len(nodes_dict)}")
    print(f"BOUN Nsets found: {len([k for k in nsets if k.startswith('BOUN_')])}")

    # Pre-correct boundary nodes for perfect alignment
    pre_correct_boun_nodes(nodes_dict, nsets)
    
    # Determine processing directions
    processing_directions = get_processing_directions(params)
    dir_names = [['X', 'Y', 'Z'][d-1] for d in processing_directions]
    print(f"\nProcessing directions: {', '.join(dir_names)}")
    
    # Process each direction
    updated_nodes = nodes_dict.copy()
    total_adjustments = 0
    
    for direction in processing_directions:
        print(f"\n{'='*40}")
        print(f"Direction {direction} ({['X', 'Y', 'Z'][direction-1]})")
        print(f"{'='*40}")
        
        # Group nodes with adaptive tolerance
        print("Adaptive tolerance search:")
        groups, adaptive_tolerance = group_nodes_by_direction(updated_nodes, direction, params.get('tol', 1e-6), nsets)
        print(f"Final: tol={adaptive_tolerance:.2e}, groups={len(groups)}")
        
        # Adjust coordinates
        updated_nodes, adjustments = adjust_node_coordinates(updated_nodes, groups, direction, adaptive_tolerance)
        print(f"Nodes adjusted: {adjustments}")
        total_adjustments += adjustments
        
        # Validate
        is_valid, validation_msg = validate_periodicity(updated_nodes, direction, adaptive_tolerance, nsets)
        print(validation_msg)
    
    # Write output file
    output_file = inp_file # Overwrite the input file
    
    print(f"\n{'='*60}")
    print(f"Writing output to overwrite: {output_file}")
    write_output_file(file_structure, updated_nodes, node_formats, output_file)
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total nodes: {len(nodes_dict)}")
    print(f"Total adjustments: {total_adjustments}")
    print(f"Output file: {output_file}")
    print("="*60)

if __name__ == "__main__":
    print("\nCleaning up temporary Abaqus files...")
    params = parse_parameters()
    output_dir = os.path.dirname(params.get('output_file', ''))

    # Clean current directory and output directory (if different)
    for path in set(['.', output_dir or '.']):
        if os.path.isdir(path):
            for f in os.listdir(path):
                if f.endswith(('.pyc', '.jnl', '.exception')) or 'abaqus.rpy' in f:
                    try: os.remove(os.path.join(path, f))
                    except OSError: pass
                        
    main()