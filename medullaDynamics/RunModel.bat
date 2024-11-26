@echo off
REM Activate the conda environment
CALL "C:\Users\munib\anaconda3\Scripts\activate.bat" ephys-torch

REM Run the Python script with arguments passed from MATLAB
python RunModel.py %1 %2 %3 %4 %5 %6 %7 %8 %9

REM Deactivate the environment (optional, as the script will end)
CALL conda deactivate
