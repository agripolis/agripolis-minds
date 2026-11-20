@echo off
rem ---------------------------------------------------------------------------
rem Run AgriPoliS for every surrogate model in RedzoneBestModels and for every
rem red zone fraction (FL-0 .. FL-50 = RedZone_fraction 0 .. 0.5).
rem
rem For each model folder RedzoneBestModels\<WS>\<Best>\FarmLin_<timestamp> it
rem   1. points MODEL_DIR in inputfiles-RR\surrogate\model.txt to that folder
rem   2. for each FL-xx writes a scenario file (copy of scenario.txt) with
rem      Scenario = FL-xx_<timestamp> and the matching RedZone_fraction,
rem      so results go to outputfiles_FL-xx_<timestamp>
rem   3. runs agpCC.exe, console output goes to Logs_models\log_FL-xx_<timestamp>.txt
rem The original model.txt is restored at the end.
rem
rem Usage:  run_redzone_models.bat [extra agpCC args, e.g. --RUNS=3 --TEILER=10]
rem         set DRYRUN=1 first to only print what would be run
rem ---------------------------------------------------------------------------
setlocal EnableExtensions EnableDelayedExpansion

set "ROOT=%~dp0"
set "ROOT=%ROOT:~0,-1%"
set "EXE=%ROOT%\src\agpCC\x64\Release\agpCC.exe"
set "INDIR=%ROOT%\inputfiles-RR"
set "MODELS=%ROOT%\RedzoneBestModels"
set "BASE_SCEN_FILE=%INDIR%\scenario.txt"
set "MODEL_TXT=%INDIR%\surrogate\model.txt"
set "EXTRA_ARGS=%*"
set "LOGDIR=%ROOT%\Logs_models"

rem "base-name fraction" pairs
set CASES="FL-0 0" "FL-10 0.1" "FL-20 0.2" "FL-30 0.3" "FL-40 0.4" "FL-50 0.5"

if not exist "%EXE%" (
    echo ERROR: %EXE% not found. Build the Release x64 configuration first.
    exit /b 1
)

if not exist "%LOGDIR%" mkdir "%LOGDIR%"
copy /y "%MODEL_TXT%" "%MODEL_TXT%.bak" >nul

set /a N=0
set /a FAILED=0
for /d %%W in ("%MODELS%\*") do (
    for /d %%B in ("%%W\*") do (
        for /d %%M in ("%%B\FarmLin_*") do (
            if exist "%%M\saved_model.pb" (
                set "MODEL=%%~nxM"
                set "TSTAMP=!MODEL:FarmLin_=!"

                rem model.txt: same as original, only MODEL_DIR replaced
                (
                    for /f "usebackq delims=" %%L in ("%MODEL_TXT%.bak") do (
                        set "LINE=%%L"
                        if /i "!LINE:~0,9!"=="MODEL_DIR" (
                            echo MODEL_DIR	=	%%~fM
                        ) else (
                            echo(%%L
                        )
                    )
                ) > "%MODEL_TXT%"

                for %%C in (%CASES%) do (
                    for /f "tokens=1,2" %%a in ("%%~C") do (
                        set /a N+=1
                        set "SCEN=%%a_!TSTAMP!"
                        set "SCEN_FILE=scenario_!SCEN!.txt"

                        rem scenario file: base scenario with new name and fraction
                        (
                            for /f "usebackq delims=" %%L in ("%BASE_SCEN_FILE%") do (
                                set "LINE=%%L"
                                if /i "!LINE:~0,9!"=="Scenario:" (
                                    echo Scenario:	!SCEN!
                                ) else if /i "!LINE:~0,16!"=="RedZone_fraction" (
                                    echo RedZone_fraction	%%b
                                ) else (
                                    echo(%%L
                                )
                            )
                        ) > "%INDIR%\!SCEN_FILE!"

                        echo.
                        echo ===== [!N!] !SCEN!  RedZone_fraction=%%b =====
                        echo model: %%~fM
                        if "%DRYRUN%"=="1" (
                            findstr /b /c:"MODEL_DIR" "%MODEL_TXT%"
                            findstr /b /c:"Scenario:" /c:"RedZone" "%INDIR%\!SCEN_FILE!"
                            echo "%EXE%" "%INDIR%" "!SCEN_FILE!" %EXTRA_ARGS%
                        ) else (
                            "%EXE%" "%INDIR%" "!SCEN_FILE!" %EXTRA_ARGS% > "%LOGDIR%\log_!SCEN!.txt" 2>&1
                        )
                        if errorlevel 1 (
                            echo FAILED, see Logs_models\log_!SCEN!.txt
                            set /a FAILED+=1
                        ) else (
                            echo done, output in outputfiles_!SCEN!
                        )
                        del "%INDIR%\!SCEN_FILE!"
                    )
                )
            )
        )
    )
)

move /y "%MODEL_TXT%.bak" "%MODEL_TXT%" >nul

echo.
echo Finished %N% runs, %FAILED% failed.
endlocal
