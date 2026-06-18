@echo off
echo === Running LYNXgwas ===
java -Xmx4g -cp bin Main --config config.properties %*
if %ERRORLEVEL% neq 0 (
    echo [FAIL] Run failed. See above for errors.
    exit /b 1
)
echo.
echo Opening viewer...
start output\index.html
