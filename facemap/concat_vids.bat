@REM @echo off
@REM setlocal enabledelayedexpansion

@REM rem Paths to the directories containing the videos
@REM set "CAM0_DIR=Cam0"
@REM set "CAM1_DIR=Cam1"
@REM set "OUTPUT_DIR=Output"

@REM rem Create the output directory if it doesn't exist
@REM if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

@REM rem Loop through all the files in Cam0 directory
@REM for %%f in (%CAM0_DIR%\*.avi) do (
@REM     rem Extract the filename without the extension
@REM     set "filename=%%~nxf"
    
@REM     rem Extract the base name without the 'cam0' part
@REM     set "base_name=%%~nxf"
@REM     set "base_name=!base_name:cam0=!"
    
@REM     rem Create the corresponding Cam1 file path
@REM     set "cam1_file=%CAM1_DIR%\!base_name:cam_0=cam_1!"

@REM     rem Print the cam1_file variable
@REM     echo Matching file: !cam1_file!

@REM     rem Check if the corresponding file exists in Cam1 directory
@REM     if exist "!cam1_file!" (
@REM         rem Concatenate the videos vertically using the provided ffmpeg command
@REM         ffmpeg -i "%%f" -i "!cam1_file!" -filter_complex "[0:v]scale=420:trunc(ow/a/2)*2[v0];[1:v]scale=420:trunc(ow/a/2)*2[v1];[v0][v1]vstack=inputs=2" "%OUTPUT_DIR%\!filename:cam0=concatenated!.avi"
@REM     ) else (
@REM         echo "Matching file for %%f not found in Cam1 directory."
@REM     )
@REM )

@REM echo "Video concatenation complete!"
@REM pause

@echo off
setlocal enabledelayedexpansion

rem Check if the root directory is provided as an argument
if "%~1"=="" (
    echo "Please provide the root directory as an argument."
    echo "Usage: script_name.bat [root_directory]"
    exit /b 1
)

rem Set the root directory from the input argument
set "ROOT_DIR=%~1"

rem Paths to the directories containing the videos relative to the root directory
set "CAM0_DIR=%ROOT_DIR%\Cam0"
set "CAM1_DIR=%ROOT_DIR%\Cam1"
set "OUTPUT_DIR=%ROOT_DIR%\Output"

rem Create the output directory if it doesn't exist
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

rem Loop through all the files in Cam0 directory
for %%f in (%CAM0_DIR%\*.avi) do (
    rem Extract the filename without the extension
    set "filename=%%~nxf"
    
    rem Extract the base name without the 'cam0' part
    set "base_name=%%~nxf"
    set "base_name=!base_name:cam0=!"
    
    rem Create the corresponding Cam1 file path
    set "cam1_file=%CAM1_DIR%\!base_name:cam_0=cam_1!"

    rem Print the cam1_file variable
    echo Matching file: !cam1_file!

    rem Check if the corresponding file exists in Cam1 directory
    if exist "!cam1_file!" (
        rem Concatenate the videos vertically using the provided ffmpeg command
        ffmpeg -i "%%f" -i "!cam1_file!" -filter_complex "[0:v]scale=420:trunc(ow/a/2)*2[v0];[1:v]scale=420:trunc(ow/a/2)*2[v1];[v0][v1]vstack=inputs=2" "%OUTPUT_DIR%\!filename:cam0=concatenated!.avi"
    ) else (
        echo "Matching file for %%f not found in Cam1 directory."
        pause
    )
)

echo "Video concatenation complete!"
pause
