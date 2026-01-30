#!/bin/bash
# Build script for Linux/Mac
# Requires .NET 8.0 SDK: https://dotnet.microsoft.com/download

echo "Building Lens Feasibility Analyzer..."

cd FeasibilityAnalyzer
dotnet build -c Release

if [ $? -eq 0 ]; then
    echo ""
    echo "Build successful!"
    echo "Executable: FeasibilityAnalyzer/bin/Release/net8.0/LensFeasibility"
    echo ""
    echo "Run examples:"
    echo "  ./FeasibilityAnalyzer/bin/Release/net8.0/LensFeasibility --example"
else
    echo "Build failed!"
fi
