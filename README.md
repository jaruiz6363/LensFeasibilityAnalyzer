# Lens Feasibility Analyzer

A standalone C# tool for pre-design feasibility analysis of optical lens systems using eikonal-based methods. This tool helps answer "Can I build this?" before investing time in detailed lens optimization.

## Features

- **Eikonal Analysis**: Uses Walther's eikonal methods to determine fundamental performance limits
- **Sine Condition Analysis**: Checks compatibility between distortion requirements and coma correction
- **Reversibility Analysis**: Computes unavoidable aberrations for symmetric/reversible lenses
- **Working Distance Range**: Analyzes performance across the specified WD range
- **Chromatic Analysis**: Estimates color correction requirements and recommends glass pairs
- **Tolerance Sensitivity**: Predicts manufacturing tolerance contributions to wavefront error
- **Complexity Estimation**: Estimates element count and design difficulty

## Building

### Requirements
- .NET 8.0 SDK or later

### Build Commands
```bash
cd FeasibilityAnalyzer
dotnet build -c Release
```

The executable will be in `bin/Release/net8.0/LensFeasibility.exe` (Windows) or `LensFeasibility` (Linux/Mac).

## Usage

### Interactive Mode (default)
```bash
./LensFeasibility
```
Prompts for all parameters with sensible defaults.

### Run Built-in Examples
```bash
./LensFeasibility --example
```
Runs analysis on several example lens types:
1. 5× Microscope Objective
2. F-Theta Scan Lens
3. 1:1 Symmetric Relay
4. Challenging High-NA Wide-Field

### Analyze from JSON File
```bash
./LensFeasibility --file input.json
./LensFeasibility --file input.json report.txt
```
Reads parameters from JSON file, optionally saves report to specified output file.

## Input Parameters

### Basic Optical Specs
| Parameter | Description | Example |
|-----------|-------------|---------|
| `NAImage` | Image-side numerical aperture | 0.25 |
| `FieldImage` | Image semi-field height (mm) | 2.5 |
| `Magnification` | Paraxial magnification (negative for real image) | -5.0 |

### Working Distance
| Parameter | Description | Example |
|-----------|-------------|---------|
| `WDNominal` | Nominal working distance (mm) | 10.0 |
| `WDMin` | Minimum working distance (mm) | 8.0 |
| `WDMax` | Maximum working distance (mm) | 15.0 |

### Distortion
| Parameter | Description | Example |
|-----------|-------------|---------|
| `DistortionType` | None, FTheta, Orthographic, Barrel, Pincushion | "None" |
| `DistortionMax` | Maximum allowed distortion (%) | 1.0 |

### Constraints
| Parameter | Description | Example |
|-----------|-------------|---------|
| `IsReversible` | Must lens work in both directions? | false |
| `TelecentricObject` | Object-space telecentric? | false |
| `TelecentricImage` | Image-space telecentric? | true |
| `FlatField` | Require flat field? | true |

### Wavelength
| Parameter | Description | Example |
|-----------|-------------|---------|
| `WavelengthPrimary` | Primary wavelength (nm) | 550.0 |
| `WavelengthMin` | Minimum wavelength (nm) | 486.1 |
| `WavelengthMax` | Maximum wavelength (nm) | 656.3 |

### Tolerance & Performance
| Parameter | Description | Example |
|-----------|-------------|---------|
| `ToleranceGrade` | 1=Commercial, 2=Precision, 3=High-Precision | 2 |
| `TargetRMSWavefront` | Target RMS wavefront error (waves) | 0.07 |

## Sample JSON Input

```json
{
  "NAImage": 0.25,
  "FieldImage": 2.5,
  "Magnification": -5.0,
  "WDNominal": 10.0,
  "WDMin": 8.0,
  "WDMax": 15.0,
  "DistortionType": "None",
  "DistortionMax": 1.0,
  "IsReversible": false,
  "TelecentricObject": false,
  "TelecentricImage": true,
  "FlatField": true,
  "WavelengthPrimary": 550.0,
  "WavelengthMin": 486.1,
  "WavelengthMax": 656.3,
  "ToleranceGrade": 2,
  "TargetRMSWavefront": 0.07
}
```

## Output Report Sections

1. **Overall Assessment**: Feasibility score (0-100) and limiting factors
2. **Derived Quantities**: Computed object-side NA, field angles, etc.
3. **Sine Condition Analysis**: Distortion vs coma compatibility
4. **Reversibility Analysis**: Unavoidable aberrations from symmetry
5. **Working Distance Analysis**: Performance across WD range
6. **Telecentricity Analysis**: Constraints and conflicts
7. **Petzval Analysis**: Field curvature correction requirements
8. **Chromatic Analysis**: Color correction level and glass recommendations
9. **Tolerance Analysis**: Sensitivity estimates and yield prediction
10. **Aberration Floors**: Theoretical minimum achievable aberrations
11. **Complexity Estimate**: Expected element count and cost class
12. **Recommendations**: Design guidance and relaxation suggestions

## Theory Background

This tool is based on eikonal methods from:

- A. Walther, "Systematic Approach to the Teaching of Lens Theory" (1967)
- A. Walther, "Lenses, Wave Optics, and Eikonal Functions" (1969)
- A. Walther, "The perils of symmetry" (1985)

The eikonal function describes the optical path length as a function of object and image coordinates. By analyzing the ideal eikonal for perfect imaging and comparing it to constraints (reversibility, distortion type, etc.), we can determine fundamental limits before detailed design.

## Tolerance Grades

| Grade | Name | Radius | Thickness | Decenter | Tilt |
|-------|------|--------|-----------|----------|------|
| 1 | Commercial | ±0.5% | ±0.10mm | ±0.05mm | ±3' |
| 2 | Precision | ±0.2% | ±0.05mm | ±0.02mm | ±1' |
| 3 | High Precision | ±0.1% | ±0.02mm | ±0.01mm | ±0.5' |

## License

This tool is provided for educational and professional use in optical design.

## Version History

- 1.0.0: Initial release with full analysis suite
