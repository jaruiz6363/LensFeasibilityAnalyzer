@echo off
REM Build script for Windows
REM Requires .NET 8.0 SDK: https://dotnet.microsoft.com/download

echo Building Lens Feasibility Analyzer...

cd FeasibilityAnalyzer
dotnet build -c Release

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Build successful!
    echo Executable: FeasibilityAnalyzer\bin\Release\net8.0\LensFeasibility.exe
    echo.
    echo Run examples:
    echo   FeasibilityAnalyzer\bin\Release\net8.0\LensFeasibility.exe --example
) else (
    echo Build failed!
)

pause
