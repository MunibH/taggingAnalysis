@echo off
setlocal enabledelayedexpansion

rem cam0 and cam1 vid filenames need to be matched to concat, but they have differing times sometimes by one or two seconds
rem wasn't as easy to account for in this batch script, just using side cam.

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

rem Function to convert time to total seconds
setlocal enableextensions enabledelayedexpansion
set "seconds="
set "time_string="
set "total_seconds="
for /f "tokens=1-3 delims=_" %%a in ("%time_string%") do (
    set /a "total_seconds = %%a * 3600 + %%b * 60 + %%c"
)
endlocal & set "%seconds%=%total_seconds%"

rem Function to calculate absolute difference between two numbers
rem Output: absolute difference in variable 'diff'
setlocal enabledelayedexpansion
set "number1="
set "number2="
set "diff="
set /a "diff = number1 - number2"
if !diff! lss 0 (
    set /a "diff = -diff"
)
endlocal & set "%diff%=!diff!"

rem Loop through all the files in Cam0 directory
for %%f in ("%CAM0_DIR%\*.avi") do (
    rem Extract the filename without the extension
    set "filename=%%~nxf"
    set "basename=%%~nf"  rem Filename without extension

    rem Extract time information from Cam0 file
    for /f "tokens=3-6 delims=_." %%t in ("%%~nf") do (
        set "time_cam0=%%t_%%u_%%v_%%w"
        call :time_to_seconds !time_cam0! time_cam0_seconds
    )

    rem Initialize variables for finding the closest match
    set "min_diff=999999"
    set "closest_cam1_file="

    rem Loop through all the files in Cam1 directory to find the closest match
    for %%g in ("%CAM1_DIR%\*.avi") do (
        rem Extract time information from Cam1 file
        for /f "tokens=3-6 delims=_." %%t in ("%%~ng") do (
            set "time_cam1=%%t_%%u_%%v_%%w"
            call :time_to_seconds !time_cam1! time_cam1_seconds
        )

        rem Calculate the absolute time difference
        call :absolute_diff !time_cam0_seconds! !time_cam1_seconds! time_diff

        rem Check if this is the closest match
        if !time_diff! lss !min_diff! (
            set "min_diff=!time_diff!"
            set "closest_cam1_file=%%g"
        )
    )

    rem Check if a matching file was found
    if not "!closest_cam1_file!"=="" (
        set "output_file=%OUTPUT_DIR%\!basename:cam0=concatenated!.avi"
        echo "Concatenating %%f and !closest_cam1_file! -> !output_file!"
        ffmpeg -i "%%f" -i "!closest_cam1_file!" -filter_complex "[0:v]scale=420:trunc(ow/a/2)*2[v0];[1:v]scale=420:trunc(ow/a/2)*2[v1];[v0][v1]vstack=inputs=2" "!output_file!"
        
        if !ERRORLEVEL! NEQ 0 (
            echo "Failed to concatenate videos for: %%f and !closest_cam1_file!"
        ) else (
            echo "Successfully concatenated to: !output_file!"
        )
    ) else (
        echo "No close matching file found for %%f in Cam1 directory."
    )
)

goto :eof

rem Function to convert time in "hour_minute_second" format to total seconds
:time_to_seconds
setlocal enabledelayedexpansion
set "time_string=%~1"
set /a "total_seconds=0"
for /f "tokens=1-3 delims=_" %%a in ("%time_string%") do (
    set /a "total_seconds = %%a * 3600 + %%b * 60 + %%c"
)
endlocal & set "%2=%total_seconds%"
goto :eof

rem Function to calculate absolute difference between two numbers
:absolute_diff
setlocal enabledelayedexpansion
set /a "diff = %1 - %2"
if !diff! lss 0 (
    set /a "diff = -diff"
)
endlocal & set "%3=!diff!"
goto :eof

echo "Video concatenation complete!"
pause
