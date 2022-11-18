setlocal

SET BUILD=mayabuild
SET MAYA_VERSION=2020
SET COMPILER=Visual Studio 16 2019

SET PFX=%~dp0
cd %PFX%
rmdir %BUILD% /s /q
mkdir %BUILD%
cd %BUILD%

cmake ^
    -DMAYA_VERSION=%MAYA_VERSION% ^
    -G "%COMPILER%" ..\

cmake --build . --config Release

pause
