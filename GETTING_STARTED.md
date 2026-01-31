# Getting Started with Lens Feasibility Analyzer

This guide helps users who are new to C# or .NET run the Lens Feasibility Analyzer.

## Step 1: Install .NET SDK

The Lens Feasibility Analyzer requires .NET 8.0 or later.

### Windows

1. Go to https://dotnet.microsoft.com/download/dotnet/8.0
2. Click "Download .NET SDK x64" (under the SDK column, not Runtime)
3. Run the downloaded installer
4. Follow the installation prompts (default options are fine)

### macOS

**Option A: Download installer**
1. Go to https://dotnet.microsoft.com/download/dotnet/8.0
2. Download the macOS installer (.pkg file)
3. Run the installer

**Option B: Using Homebrew**
```bash
brew install dotnet-sdk
```

### Linux (Ubuntu/Debian)

```bash
# Add Microsoft package repository
wget https://packages.microsoft.com/config/ubuntu/22.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb
sudo dpkg -i packages-microsoft-prod.deb
rm packages-microsoft-prod.deb

# Install .NET SDK
sudo apt-get update
sudo apt-get install -y dotnet-sdk-8.0
```

### Verify Installation

Open a terminal (Command Prompt on Windows, Terminal on macOS/Linux) and run:

```bash
dotnet --version
```

You should see a version number like `8.0.xxx`. If you see an error, restart your terminal or computer.

## Step 2: Download the Source Code

### Option A: Using Git

If you have Git installed:

```bash
git clone https://github.com/jaruiz6363/LensFeasibilityAnalyzer.git
cd LensFeasibilityAnalyzer
```

### Option B: Download ZIP

1. Go to the GitHub repository page
2. Click the green "Code" button
3. Click "Download ZIP"
4. Extract the ZIP file to a folder of your choice

## Step 3: Build the Program

Open a terminal and navigate to the project folder:

```bash
cd path/to/LensFeasibilityAnalyzer/FeasibilityAnalyzer
```

Build the program:

```bash
dotnet build -c Release
```

You should see "Build succeeded" with 0 errors.

## Step 4: Run the Program

### Run Built-in Examples

The easiest way to see the program in action:

```bash
dotnet run -c Release -- --example
```

This runs analysis on four example lens systems.

### Interactive Mode

To enter your own lens parameters interactively:

```bash
dotnet run -c Release
```

The program will prompt you for each parameter with sensible defaults.

### Analyze from JSON File

Create a JSON file with your lens parameters (see `README.md` for format), then:

```bash
dotnet run -c Release -- --file your_input.json
```

To save the report to a file:

```bash
dotnet run -c Release -- --file your_input.json output_report.txt
```

## Troubleshooting

### "dotnet: command not found"

- Make sure you installed the .NET SDK (not just the Runtime)
- Restart your terminal after installation
- On Windows, you may need to restart your computer

### "The project file could not be found"

Make sure you're in the correct directory. You should be in the `FeasibilityAnalyzer` folder that contains `FeasibilityAnalyzer.csproj`.

### Build errors about .NET version

If you see errors about unsupported .NET version, make sure you have .NET 8.0 or later:

```bash
dotnet --list-sdks
```

If you have an older version, download .NET 8.0 from the link above.

### Permission denied (Linux/macOS)

If you get permission errors, you may need to use `sudo` for installation steps, but NOT for running `dotnet build` or `dotnet run`.

## Quick Reference

| Command | Description |
|---------|-------------|
| `dotnet build -c Release` | Build the program |
| `dotnet run -c Release` | Run in interactive mode |
| `dotnet run -c Release -- --example` | Run built-in examples |
| `dotnet run -c Release -- --file input.json` | Analyze from JSON file |
| `dotnet run -c Release -- --help` | Show help message |

## Creating a Standalone Executable

If you want to create an executable that doesn't require .NET to be installed:

```bash
# Windows
dotnet publish -c Release -r win-x64 --self-contained true

# macOS (Intel)
dotnet publish -c Release -r osx-x64 --self-contained true

# macOS (Apple Silicon)
dotnet publish -c Release -r osx-arm64 --self-contained true

# Linux
dotnet publish -c Release -r linux-x64 --self-contained true
```

The executable will be in `bin/Release/net8.0/{runtime}/publish/`.

## Need More Help?

- See `README.md` for detailed documentation on input parameters and output format
- Check the `Research` folder for reference papers on the underlying theory
- Open an issue on GitHub if you encounter problems
