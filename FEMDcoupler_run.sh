#!/bin/bash

JSON_FILE=$(python3 FEMDcoupler_params.py)
ABAQUS_MODE=$(python3 -c "import json; print('noGUI' if json.load(open('$JSON_FILE'))['skip_abaqus_GUI'] else 'script')")
x-terminal-emulator -e "bash -c 'python3 create_anchors_cut_notch.py \"$JSON_FILE\"'" &
abaqus cae $ABAQUS_MODE=create_part_mesh.py -- "$JSON_FILE"
python3 create_periodicity_mesh.py "$JSON_FILE"
python3 create_merge_FEMD_data_inp.py "$JSON_FILE"

read -p "Press Enter to exit."