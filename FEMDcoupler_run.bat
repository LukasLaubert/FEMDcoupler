@echo off
setlocal

FOR /F "tokens=*" %%F IN ('python FEMDcoupler_params.py') DO (set "JSON_FILE=%%F")
start "MD data file processing: anchor/cut/notch handling" python create_anchors_cut_notch.py %JSON_FILE%
start "" /B /WAIT cmd /c "abaqus cae script=create_part_mesh.py -- %JSON_FILE%" > nul
python create_periodicity_mesh.py %JSON_FILE%
python create_merge_FEMD_data_inp.py %JSON_FILE%

endlocal
pause