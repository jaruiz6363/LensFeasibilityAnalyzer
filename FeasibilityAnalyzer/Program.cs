using System;
using System.IO;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace FeasibilityAnalyzer
{
    /// <summary>
    /// Main entry point for the Lens Feasibility Analyzer
    /// </summary>
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine();
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║              LENS FEASIBILITY ANALYZER v1.1                              ║");
            Console.WriteLine("║              Eikonal-Based Pre-Design Assessment Tool                    ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();

            if (args.Length == 0)
            {
                ShowUsage();
                RunInteractiveMode();
            }
            else if (args[0] == "--help" || args[0] == "-h")
            {
                ShowUsage();
            }
            else if (args[0] == "--example" || args[0] == "-e")
            {
                RunExamples();
            }
            else if (args[0] == "--file" || args[0] == "-f")
            {
                if (args.Length < 2)
                {
                    Console.WriteLine("Error: Please specify input file path");
                    return;
                }
                RunFromFile(args[1], args.Length > 2 ? args[2] : null);
            }
            else if (args[0] == "--interactive" || args[0] == "-i")
            {
                RunInteractiveMode();
            }
            else
            {
                Console.WriteLine($"Unknown option: {args[0]}");
                ShowUsage();
            }
        }

        static void ShowUsage()
        {
            Console.WriteLine("Usage:");
            Console.WriteLine("  LensFeasibility                     Run interactive mode");
            Console.WriteLine("  LensFeasibility -e, --example       Run built-in examples");
            Console.WriteLine("  LensFeasibility -f, --file <path>   Analyze from JSON input file");
            Console.WriteLine("  LensFeasibility -i, --interactive   Run interactive parameter entry");
            Console.WriteLine("  LensFeasibility -h, --help          Show this help");
            Console.WriteLine();
            Console.WriteLine("Output options (after -f):");
            Console.WriteLine("  LensFeasibility -f input.json output.txt   Save report to file");
            Console.WriteLine();
        }

        static void RunExamples()
        {
            var engine = new FeasibilityEngine();

            // Example 1: Microscope objective
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine("EXAMPLE 1: 5× Microscope Objective");
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine();

            var input1 = FeasibilityInput.CreateMicroscopeObjective5x();
            var output1 = engine.Analyze(input1);
            Console.WriteLine(ReportGenerator.GenerateReport(output1));

            Console.WriteLine();
            Console.WriteLine("Press Enter to continue to next example...");
            Console.ReadLine();

            // Example 2: F-theta scan lens
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine("EXAMPLE 2: F-Theta Scan Lens");
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine();

            var input2 = FeasibilityInput.CreateFThetaScanLens();
            var output2 = engine.Analyze(input2);
            Console.WriteLine(ReportGenerator.GenerateReport(output2));

            Console.WriteLine();
            Console.WriteLine("Press Enter to continue to next example...");
            Console.ReadLine();

            // Example 3: Symmetric relay lens
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine("EXAMPLE 3: 1:1 Symmetric Relay Lens");
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine();

            var input3 = FeasibilityInput.CreateRelayLens();
            var output3 = engine.Analyze(input3);
            Console.WriteLine(ReportGenerator.GenerateReport(output3));

            // Example 4: Challenging high-NA wide-field
            Console.WriteLine();
            Console.WriteLine("Press Enter to continue to next example...");
            Console.ReadLine();

            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine("EXAMPLE 4: Challenging High-NA Wide-Field Lens");
            Console.WriteLine("═══════════════════════════════════════════════════════════════════════════");
            Console.WriteLine();

            var input4 = new FeasibilityInput
            {
                NAImage = 0.5,
                FieldImage = 15.0,
                Magnification = -10.0,
                WDNominal = 5.0,
                WDMin = 3.0,
                WDMax = 10.0,
                DistortionType = DistortionType.None,
                DistortionMax = 0.5,
                IsReversible = false,
                TelecentricObject = false,
                TelecentricImage = true,
                FlatField = true,
                WavelengthPrimary = 550.0,
                WavelengthMin = 400.0,
                WavelengthMax = 700.0,
                ToleranceGrade = 3,
                TargetRMSWavefront = 0.05
            };
            var output4 = engine.Analyze(input4);
            Console.WriteLine(ReportGenerator.GenerateReport(output4));

            Console.WriteLine();
            Console.WriteLine("Examples complete.");
        }

        static void RunFromFile(string inputPath, string outputPath)
        {
            try
            {
                if (!File.Exists(inputPath))
                {
                    Console.WriteLine($"Error: Input file not found: {inputPath}");
                    return;
                }

                string json = File.ReadAllText(inputPath);
                var options = new JsonSerializerOptions
                {
                    PropertyNameCaseInsensitive = true,
                    Converters = { new JsonStringEnumConverter() }
                };

                var input = JsonSerializer.Deserialize<FeasibilityInput>(json, options);
                
                if (input == null)
                {
                    Console.WriteLine("Error: Failed to parse input file");
                    return;
                }

                var validation = input.Validate();
                if (!validation.IsValid)
                {
                    Console.WriteLine("Input validation errors:");
                    foreach (var error in validation.Errors)
                        Console.WriteLine($"  - {error}");
                    return;
                }

                var engine = new FeasibilityEngine();
                var output = engine.Analyze(input);
                string report = ReportGenerator.GenerateReport(output);

                if (!string.IsNullOrEmpty(outputPath))
                {
                    File.WriteAllText(outputPath, report);
                    Console.WriteLine($"Report saved to: {outputPath}");
                    Console.WriteLine();
                    Console.WriteLine(ReportGenerator.GenerateSummary(output));
                }
                else
                {
                    Console.WriteLine(report);
                }
            }
            catch (JsonException ex)
            {
                Console.WriteLine($"JSON parsing error: {ex.Message}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error: {ex.Message}");
            }
        }

        static void RunInteractiveMode()
        {
            Console.WriteLine("Interactive Parameter Entry");
            Console.WriteLine("(Press Enter to use default values shown in brackets)");
            Console.WriteLine();

            var input = new FeasibilityInput();

            // Basic parameters
            Console.WriteLine("─── Basic Optical Parameters ───");
            input.NAImage = ReadDouble("Image NA", 0.25, 0.01, 0.99);
            input.FieldImage = ReadDouble("Image semi-field (mm)", 5.0, 0.1, 100);
            input.Magnification = ReadDouble("Magnification (negative for real image)", -5.0, -100, 100);

            // Working distance
            Console.WriteLine();
            Console.WriteLine("─── Working Distance ───");
            input.WDNominal = ReadDouble("Nominal WD (mm)", 20.0, 0.1, 10000);
            input.WDMin = ReadDouble("Minimum WD (mm)", input.WDNominal * 0.8, 0.1, input.WDNominal);
            input.WDMax = ReadDouble("Maximum WD (mm)", input.WDNominal * 1.5, input.WDNominal, 10000);

            // Distortion
            Console.WriteLine();
            Console.WriteLine("─── Distortion ───");
            Console.WriteLine("  Distortion types: 0=None, 1=FTheta, 2=Orthographic, 3=Barrel, 4=Pincushion");
            int distType = (int)ReadDouble("Distortion type", 0, 0, 4);
            input.DistortionType = (DistortionType)distType;
            if (input.DistortionType != DistortionType.None)
            {
                input.DistortionMax = ReadDouble("Max distortion (%)", 1.0, 0, 50);
            }

            // Constraints
            Console.WriteLine();
            Console.WriteLine("─── Constraints ───");
            input.IsReversible = ReadBool("Reversible lens", false);
            input.TelecentricObject = ReadBool("Object telecentric", false);
            input.TelecentricImage = ReadBool("Image telecentric", false);
            input.FlatField = ReadBool("Flat field required", true);

            // Wavelength
            Console.WriteLine();
            Console.WriteLine("─── Wavelength ───");
            input.WavelengthPrimary = ReadDouble("Primary wavelength (nm)", 550.0, 200, 2000);
            input.WavelengthMin = ReadDouble("Min wavelength (nm)", 486.1, 200, input.WavelengthPrimary);
            input.WavelengthMax = ReadDouble("Max wavelength (nm)", 656.3, input.WavelengthPrimary, 2000);

            // Tolerance
            Console.WriteLine();
            Console.WriteLine("─── Tolerance & Performance ───");
            Console.WriteLine("  Tolerance grades: 1=Commercial, 2=Precision, 3=High-Precision");
            input.ToleranceGrade = (int)ReadDouble("Tolerance grade", 2, 1, 3);
            input.TargetRMSWavefront = ReadDouble("Target RMS wavefront (waves)", 0.07, 0.01, 1.0);

            // Validate
            var validation = input.Validate();
            if (!validation.IsValid)
            {
                Console.WriteLine();
                Console.WriteLine("Input validation errors:");
                foreach (var error in validation.Errors)
                    Console.WriteLine($"  - {error}");
                return;
            }

            // Run analysis
            Console.WriteLine();
            Console.WriteLine("Running analysis...");
            Console.WriteLine();

            var engine = new FeasibilityEngine();
            var output = engine.Analyze(input);

            // Show report
            string report = ReportGenerator.GenerateReport(output);
            Console.WriteLine(report);

            // Offer to save
            Console.WriteLine();
            Console.Write("Save report to file? (Enter filename or press Enter to skip): ");
            string filename = Console.ReadLine()?.Trim();
            if (!string.IsNullOrEmpty(filename))
            {
                try
                {
                    File.WriteAllText(filename, report);
                    Console.WriteLine($"Report saved to: {filename}");

                    // Also save input as JSON for future reference
                    string jsonFile = Path.ChangeExtension(filename, ".json");
                    var jsonOptions = new JsonSerializerOptions 
                    { 
                        WriteIndented = true,
                        Converters = { new JsonStringEnumConverter() }
                    };
                    string inputJson = JsonSerializer.Serialize(input, jsonOptions);
                    File.WriteAllText(jsonFile, inputJson);
                    Console.WriteLine($"Input saved to: {jsonFile}");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"Error saving file: {ex.Message}");
                }
            }
        }

        static double ReadDouble(string prompt, double defaultValue, double min, double max)
        {
            while (true)
            {
                Console.Write($"  {prompt} [{defaultValue}]: ");
                string input = Console.ReadLine()?.Trim();

                if (string.IsNullOrEmpty(input))
                    return defaultValue;

                if (double.TryParse(input, out double value))
                {
                    if (value >= min && value <= max)
                        return value;
                    Console.WriteLine($"    Value must be between {min} and {max}");
                }
                else
                {
                    Console.WriteLine("    Invalid number format");
                }
            }
        }

        static bool ReadBool(string prompt, bool defaultValue)
        {
            string defaultStr = defaultValue ? "y" : "n";
            Console.Write($"  {prompt} (y/n) [{defaultStr}]: ");
            string input = Console.ReadLine()?.Trim().ToLower();

            if (string.IsNullOrEmpty(input))
                return defaultValue;

            return input == "y" || input == "yes" || input == "true" || input == "1";
        }

        /// <summary>
        /// Generate a sample JSON input file
        /// </summary>
        public static void GenerateSampleInputFile(string path)
        {
            var input = FeasibilityInput.CreateMicroscopeObjective5x();
            var options = new JsonSerializerOptions 
            { 
                WriteIndented = true,
                Converters = { new JsonStringEnumConverter() }
            };
            string json = JsonSerializer.Serialize(input, options);
            File.WriteAllText(path, json);
            Console.WriteLine($"Sample input file created: {path}");
        }
    }
}
