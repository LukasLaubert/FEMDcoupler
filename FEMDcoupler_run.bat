@echo off
setlocal

FOR /F "tokens=*" %%F IN ('python FEMDcoupler_params.py') DO (set "JSON_FILE=%%F")
FOR /F "tokens=*" %%G IN ('python -c "import json; print(json.load(open('%JSON_FILE%'))['skip_abaqus_GUI'])"') DO (set "SKIP_GUI=%%G")
if "%SKIP_GUI%"=="True" (set "ABAQUS_MODE=noGUI") else (set "ABAQUS_MODE=script")
start "MD data file processing: anchor/cut/notch handling" python create_anchors_cut_notch.py %JSON_FILE%
start "" /B /WAIT cmd /c "abaqus cae %ABAQUS_MODE%=create_part_mesh.py -- %JSON_FILE%"
python create_periodicity_mesh.py %JSON_FILE%
python create_merge_FEMD_data_inp.py %JSON_FILE%

endlocal
pause