# ==============================================================================
# Parameter Class
# ==============================================================================
class ParamConfig:
    # --- LEGEND / ABBREVIATIONS ---
    # FE: Finite Element       -> The surrounding continuum part of the domain
    # MD: Molecular Dynamics   -> The particle-based core of the domain
    # BD: Bridging Domain      -> The transition/overlap region between FE and MD
    
    # --- PATH & FILE HANDLING ----
    lammps_file = "../equilSystems/epoxy_anisoNPT_20M.data"  # specify a file name with (optional) a path; if "None" (without ""), uses first file with ".data" extension in current path
    output_file = "../epoxy_anisoNPT_20M_Ycut30d"  # file name of the output files excluding file extension with (optional) output path

    # --- FE GEOMETRY PARAMETERS ---
    # (a) Rectangular FE systam: set 'circular_burrito_radius = None' below to use
    domain_min_coords = [80, 80, None]  # X_min, Y_min, Z_min of FE domain (or None for PBCs in this direction)
    domain_max_coords = [320, 230, None]  # X_max, Y_max, Z_max of FE domain (or None for PBCs in this direction)
    # (b) Circular FE system
    circular_burrito_radius = None      # radius of a circular burrito system from MD center | !!overwrites '(a) rectuangular FE system' if not None!!
    pbc_direction = 2                   # direction of periodic boundary conditions; 1: x, 2: y, 3: z
    

    # --- FE REGION & PARTITION PARAMETERS ---
    bridge_thickness = 10           # FE extends this much into MD box, creating the BD
    default_mesh_size = 10          # default mesh size, which is automatically adapted to match the domain sizes well
#    mesh_coarsening_factor = 3     # coarsens mesh with values >1 by the factor given up to the outermost border; refines mesh with 0<values<1; disable with None or 0 or 1 # NOT IMPLEMENTED
    tol = 1e-6

    full_elem_bridge = True         # if True, FE elements are always fully covered with MD particles in the BD (mesh size for BD is seperately calculated from rest of FE); if False, MD can also extend into only parts of FE elements (BD mesh size is calculated together with rest of FE)
    bridge_cellset = True           # if True, BD is a seperate cell set; if False, BD is part of the FE cell set (and no separate BD cellset is created)
    have_padding = False            # if True, creates a padding layer of FE elements whereever there is no FE(/BD) layer (i.e. where we have PBCs in the MD system)
    
    mesh_md_region = False          # if True, meshes the MD region as well
    mesh_cutout_region = False      # if True, meshes the "empty" (cut/notch) region as well
    have_iface = False              # if True, creates a nodeset at the boundary between the FE(/BD) and the pure MD partition
    
    dPBC_correction = True          # if True, postprocesses the mesh outside of Abaqus to guarantee that periodic Dirichlet boundary conditions can be imposed seamlessly; overwrites output .inp file from Abaqus

    # --- MD PARAMETERS ---
    dpd_thickness = 12                  # all particles reaching from MD border up to this thickness in each nonperiodic boundary direction are converted to seperate types; set to 0 or None for disabling DPD conversion
    anchor_placement_probability = 0.8  # probability of putting an anchor atom on a particle in the bridge_thickness region; set to 0 for having no anchors
    ANCHOR_ATOM_MASS_VAL = 1.0e-50      # mass of the anchor atoms (must be non-zero)
    truncate_chains = True              # if True, truncates border chains in directions of specified domain_min_coords and/or domain_max_coords (where we have FE) (this is essential wherever we do not have PBCs in the MD region)
    remove_molecules_threshold = 3      # threshold of atoms per molecule: molecules with <= this number of atoms (after cut) are deleted entirely. None disables this.
    remove_bridge_interactions = False  # if True, removes all bonded, angle, dihedral, and improper interactions between atoms within bridging region if ALL involved atoms are inside brdige (excludes anchor bonds)

    # --- FRACTURE PARAMETERS ---
    # - GENERAL -
    md_remove_shape = 'ellipse'   # cutout shape in MD: 'ellipse', 'cube' | holds for cutout and notch

    # - CUTOUT -
    # Defines a void/shape removed from the inside of the model.
    # Width defines the diameter of the ellipsoid (or side length) in that direction.
    cut_width = [10, 20, None]    # width of cut through the whole model in the specified plane(s) [x, y, z] | 0 will set the width to tol, negative value will be used similar to positive values | None means "no cut" in this plane, resulting in an infinite extension in this direction
    cut_shift = [-25, 20, None]   # shift of the cut from the center in the direction of the specified plane [x, y, z] | None will be treated as 0 | can only be defined if a cut_width is defined in the respective direction

    # - NOTCH -
    # Defines a cut starting from the outer surface/edge of the part.
    notch_width = 12        # width of a centered notch along the specified plane (positive and negative values lead to same result; None means "no notch", function is skipped; 0 will result in width of tol)
    notch_normal = 1        # normal direction of the plane where the notch is cut through; the cut is along that plane; 1: x, 2: y; 3: z
    notch_plane = -2        # part's outer surface plane at which the notch starts to be cut through the whole part (cannot be equal to +/-notch_normal); 1: pos_x, -1 neg_x, 2: pos_y, -2: neg_y, 3: pos_z, -3: neg_z
    notch_shift = 0         # shift perpendicular to the normal direction of the cut notch; "None" is interpreted as 0
    notch_depth = None      # depth of the notch measured from outer notch_plane FE surface (going inside the part no matter if positive or negative); "None" specifies that the notch ranges until the middle of the MD region



# ==============================================================================
# Let's go
# ==============================================================================
import os
import json

def create_json_config():
    """Serializes the ParamConfig class to a JSON file and prints the filename."""    
    params_dict = {
        key: value for key, value in vars(ParamConfig).items()
        if not key.startswith('__')
    }
    
    json_filename = "{}.json".format(ParamConfig.output_file)

    with open(json_filename, 'w') as f:
        json.dump(params_dict, f, indent=4)
        
    print(json_filename)

if __name__ == "__main__":
    create_json_config()