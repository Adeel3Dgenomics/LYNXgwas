@echo off
echo === Running LYNXgwas ===

REM Build classpath with lib dependencies
set "CP=bin"
for %%f in (lib\*.jar) do call set "CP=%%CP%%;%%f"

java -Xmx8g -cp "%CP%" Main %*
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Run failed. See above for errors.
    exit /b 1
)
echo.
echo Opening home page...
start http://localhost:8765/
