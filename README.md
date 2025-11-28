# FEMDcoupler

FEMDcoupler is a pre-processing toolkit for generating coupled Finite Element (FE) and Molecular Dynamics (MD) models. It creates simulation-ready systems by integrating an MD core with a surrounding FE domain, supporting features like pre-cracks and periodic boundaries.

The workflow is orchestrated through a series of Python scripts that use Abaqus for FE mesh generation and perform data manipulation on LAMMPS data files.

## Features

- **FE & MD Model Generation**: Prepares 3D models with a central MD region integrated into a surrounding FE continuum based on a single LAMMPS data file. Supports both rectangular and cylindrical FE domains.
- **Meshing**: Generates a 3D hexahedral mesh. The mesh density is configurable and adapts to the geometry specified and given by the LAMMPS data file.
- **Fracture Mechanics**: Optionally introduces pre-cracks into the model.
    - **Notches**: Create notches with configurable thickness, depth, and orientation.
    - **Cutouts**: Insert elliptical or cuboid cutouts, either as a cavity/void or a groove through the entire model.
- **Periodic Boundary Support**: Post-processes the Abaqus-generated mesh to ensure node-for-node periodicity, which allows for applying periodic Dirichlet boundary conditions (dPBCs).
- **FE-MD Interface Handling**:
    - **Anchor Atoms**: Probabilistically places anchor atoms in the bridging region between the FE and MD domains to facilitate coupling.
    - **Interaction Management**: Can automatically remove interactions at the boundary of the bridging region to avoid interactions between particles at non-periodic boundaries.
    - **Chain Truncation**: Optional truncation of molecular chains at the FE-MD boundary to ensure clean coupling interfaces.
    - **DPD Support**: Capabilities to re-type atoms near non-periodic boundaries, facilitating Dissipative Particle Dynamics (DPD) setups.
- **Configuration**: A centralized Python configuration file (`FEMDcoupler_params.py`) allows for customization of the model, from geometry and mesh size to pre-crack parameters. This script generates a JSON configuration file used by the subsequent pipeline.

## Software Requirements

- **Abaqus/CAE** (2023 or older): Required for FE part creation and meshing. The scripts are designed to interface with the Abaqus Python environment.
- **Python 3.x**: A standard Python 3 interpreter is needed to run the data processing and orchestration scripts.

## Installation

1.  Clone or download the repository to your local machine.
2.  Ensure that the `abaqus` command is available in your system's PATH.
3.  On Linux/macOS, make the shell script executable before the first run:
    ```bash
    chmod +x FEMDcoupler_run.sh
    ```

## Usage

The primary workflow is managed by the main run scripts, which execute a sequence of Python scripts in the correct order.

### 1. Configuration

All model parameters are defined within the `ParamConfig` class in the **`FEMDcoupler_params.py`** file. Before running the workflow, modify this file to define your desired model.

Key configurable aspects include:
- FE domain geometry (rectangular or circular).
- Meshing and cell set parameters, as well as bridging region thickness.
- Cut and notch parameters for fracture studies.
- Anchor atom placement probability.
- Input and output file paths.

Refer to the comments in `FEMDcoupler_params.py` for a detailed explanation of each parameter.

### 2. Execution

- **On Windows**:
  ```batch
  .\FEMDcoupler_run.bat
  ```
- **On Linux/macOS**:
  ```bash
  ./FEMDcoupler_run.sh
  ```

Running these scripts will initiate the full pre-processing pipeline, culminating in the generation of the final coupled model files.

### Tips & Error Handling

- **Running Scripts Individually**: All Python scripts can be run independently. For a full process, simply follow the execution order stated under "Workflow Overview" below, which is also executed by the `.bat` or `.sh` files.

- **Abaqus GUI Interaction**:
    - When the `create_part_mesh.py` script is executed, it will launch Abaqus/CAE. You can make arbitrary changes to the model within the GUI at this stage.
    - To ensure your changes are recognized by subsequently executed scripts, you must save the model by overwriting the already generated `.inp` file. To do so, expand `Jobs` in the `Analysis` section on the Model Worktree (click `+` left of it) > right click the only existing job listed there > `Write Input` > `OK` > `Yes`.
    - If no changes are made, saving is not required, as the script handles file generation automatically. Abaqus may be closed without saving.

- **Abaqus Version Incompatibilities**: `create_part_mesh.py` is written in Python 2.7, which is only compatible with Abaqus/CAE 2023 or older. For running it in Abaqus CAE 2024 or newer, please refer to this guide: https://tecnodigitalschool.com/upgrade-from-python2-to-python3-abaqus-2024/

- **Remote/Headless Abaqus**:
    - **X2Go Rendering Issues**: If running Abaqus/CAE via a remote desktop solution like X2Go, you may encounter rendering problems. To resolve this, add the ` -mesa` flag after `abaqus cae` in the `FEMDcoupler_run.sh` for starting Abaqus. If the mesh generation script is then not executed automatically, execute it manually via `File > Run script... > create_part_mesh.py` inside the GUI.
    - **Headless Mode**: To prevent the Abaqus GUI from starting altogether (e.g. if no manual changes to the FE mesh and its sets are required), you can modify the run script. After the `abaqus cae` command, replace `script=` with `noGUI=`. This is significantly faster but removes the ability to interactively modify the FE model.

## Workflow Overview

The `FEMDcoupler` toolkit follows this sequence to build the final model:

1.  **`FEMDcoupler_params.py`**: Reads the user-defined parameters and generates a `.json` configuration file.
2.  **`create_anchors_cut_notch.py`**: Modifies the source LAMMPS data file based on the JSON configuration to introduce cuts, notches, and anchor atoms as well as to truncate atom bonds at the boundaries and to re-group boundary atoms for DPD deployment.
3.  **`create_part_mesh.py`**: Launches Abaqus/CAE to build the FE part, partition it, generate the mesh, assign sets, and write an initial `.inp` file.
4.  **`create_periodicity_mesh.py`**: Post-processes the `.inp` file, adjusting node coordinates to ensure perfect periodicity for dPBCs.
5.  **`create_merge_FEMD_data_inp.py`**: Merges the FE mesh data from the (corrected) `.inp` file with the modified MD atom data for inspection of the coupled configuration with tools like OVITO.

## License

This software is distributed under the MIT License. See the `LICENSE` file for more details.