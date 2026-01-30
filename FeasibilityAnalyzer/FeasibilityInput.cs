using System;

namespace FeasibilityAnalyzer
{
    /// <summary>
    /// Distortion type enumeration
    /// </summary>
    public enum DistortionType
    {
        None,           // Zero distortion (rectilinear)
        FTheta,         // f-theta for laser scanning
        Orthographic,   // f*sin(theta)
        Barrel,         // Negative distortion
        Pincushion,     // Positive distortion
        Mustache        // Mixed barrel-pincushion
    }

    /// <summary>
    /// Input specification for feasibility analysis
    /// </summary>
    public class FeasibilityInput
    {
        // ============ Basic Optical Specs ============
        
        /// <summary>Image-side numerical aperture</summary>
        public double NAImage { get; set; }
        
        /// <summary>Image semi-field height in mm</summary>
        public double FieldImage { get; set; }
        
        /// <summary>Paraxial magnification (negative for real inverted image)</summary>
        public double Magnification { get; set; }

        // ============ Working Distance ============
        
        /// <summary>Nominal working distance in mm</summary>
        public double WDNominal { get; set; }
        
        /// <summary>Minimum working distance in mm</summary>
        public double WDMin { get; set; }
        
        /// <summary>Maximum working distance in mm</summary>
        public double WDMax { get; set; }

        // ============ Distortion ============
        
        /// <summary>Type of distortion required/allowed</summary>
        public DistortionType DistortionType { get; set; } = DistortionType.None;
        
        /// <summary>Maximum allowed distortion in percent</summary>
        public double DistortionMax { get; set; } = 1.0;

        // ============ Symmetry Constraints ============
        
        /// <summary>Must the lens work in both directions?</summary>
        public bool IsReversible { get; set; } = false;
        
        /// <summary>Object-space telecentric?</summary>
        public bool TelecentricObject { get; set; } = false;
        
        /// <summary>Image-space telecentric?</summary>
        public bool TelecentricImage { get; set; } = false;

        // ============ Field Shape ============
        
        /// <summary>Require flat field?</summary>
        public bool FlatField { get; set; } = true;
        
        /// <summary>Maximum field curvature if not flat, in mm</summary>
        public double FieldCurvatureMax { get; set; } = 1.0;

        // ============ Wavelength ============
        
        /// <summary>Primary wavelength in nm</summary>
        public double WavelengthPrimary { get; set; } = 550.0;
        
        /// <summary>Minimum wavelength in nm</summary>
        public double WavelengthMin { get; set; } = 486.1;  // F line
        
        /// <summary>Maximum wavelength in nm</summary>
        public double WavelengthMax { get; set; } = 656.3;  // C line
        
        /// <summary>Spectral bandwidth in nm (convenience property)</summary>
        public double Bandwidth => WavelengthMax - WavelengthMin;

        // ============ Additional Constraints ============
        
        /// <summary>Maximum total track length in mm (0 = no constraint)</summary>
        public double MaxTotalTrack { get; set; } = 0;
        
        /// <summary>Maximum number of elements (0 = no constraint)</summary>
        public int MaxElements { get; set; } = 0;
        
        /// <summary>Allow aspheric surfaces?</summary>
        public bool AllowAspheres { get; set; } = true;
        
        /// <summary>Target RMS wavefront error in waves</summary>
        public double TargetRMSWavefront { get; set; } = 0.07;  // Maréchal criterion

        // ============ Tolerance Requirements ============
        
        /// <summary>Tolerance grade: 1=commercial, 2=precision, 3=high-precision</summary>
        public int ToleranceGrade { get; set; } = 2;

        /// <summary>
        /// Validate the input parameters
        /// </summary>
        public ValidationResult Validate()
        {
            var result = new ValidationResult();

            if (NAImage <= 0 || NAImage >= 1)
                result.AddError("NAImage must be between 0 and 1");

            if (FieldImage <= 0)
                result.AddError("FieldImage must be positive");

            if (Math.Abs(Magnification) < 1e-10)
                result.AddError("Magnification cannot be zero");

            if (WDNominal <= 0)
                result.AddError("WDNominal must be positive");

            if (WDMin <= 0 || WDMin > WDNominal)
                result.AddError("WDMin must be positive and <= WDNominal");

            if (WDMax < WDNominal)
                result.AddError("WDMax must be >= WDNominal");

            if (WavelengthPrimary <= 0)
                result.AddError("WavelengthPrimary must be positive");

            if (WavelengthMin >= WavelengthMax)
                result.AddError("WavelengthMin must be less than WavelengthMax");

            if (DistortionMax < 0)
                result.AddError("DistortionMax cannot be negative");

            if (TargetRMSWavefront <= 0)
                result.AddError("TargetRMSWavefront must be positive");

            return result;
        }

        /// <summary>
        /// Create a copy of this input
        /// </summary>
        public FeasibilityInput Clone()
        {
            return (FeasibilityInput)this.MemberwiseClone();
        }

        /// <summary>
        /// Create example input for a 5x microscope objective
        /// </summary>
        public static FeasibilityInput CreateMicroscopeObjective5x()
        {
            return new FeasibilityInput
            {
                NAImage = 0.25,
                FieldImage = 2.5,
                Magnification = -5.0,
                WDNominal = 10.0,
                WDMin = 8.0,
                WDMax = 15.0,
                DistortionType = DistortionType.None,
                DistortionMax = 1.0,
                IsReversible = false,
                TelecentricObject = false,
                TelecentricImage = true,
                FlatField = true,
                WavelengthPrimary = 550.0,
                WavelengthMin = 486.1,
                WavelengthMax = 656.3,
                ToleranceGrade = 2
            };
        }

        /// <summary>
        /// Create example input for an f-theta scan lens
        /// </summary>
        public static FeasibilityInput CreateFThetaScanLens()
        {
            return new FeasibilityInput
            {
                NAImage = 0.02,
                FieldImage = 50.0,
                Magnification = -1.0,
                WDNominal = 100.0,
                WDMin = 100.0,
                WDMax = 100.0,
                DistortionType = DistortionType.FTheta,
                DistortionMax = 0.1,
                IsReversible = false,
                TelecentricObject = true,
                TelecentricImage = false,
                FlatField = true,
                WavelengthPrimary = 1064.0,
                WavelengthMin = 1060.0,
                WavelengthMax = 1068.0,
                ToleranceGrade = 2
            };
        }

        /// <summary>
        /// Create example input for a relay lens
        /// </summary>
        public static FeasibilityInput CreateRelayLens()
        {
            return new FeasibilityInput
            {
                NAImage = 0.15,
                FieldImage = 10.0,
                Magnification = -1.0,
                WDNominal = 50.0,
                WDMin = 50.0,
                WDMax = 50.0,
                DistortionType = DistortionType.None,
                DistortionMax = 0.5,
                IsReversible = true,
                TelecentricObject = true,
                TelecentricImage = true,
                FlatField = true,
                WavelengthPrimary = 550.0,
                WavelengthMin = 486.1,
                WavelengthMax = 656.3,
                ToleranceGrade = 2
            };
        }
    }

    /// <summary>
    /// Validation result container
    /// </summary>
    public class ValidationResult
    {
        public bool IsValid => Errors.Count == 0;
        public System.Collections.Generic.List<string> Errors { get; } = new System.Collections.Generic.List<string>();
        public System.Collections.Generic.List<string> Warnings { get; } = new System.Collections.Generic.List<string>();

        public void AddError(string message) => Errors.Add(message);
        public void AddWarning(string message) => Warnings.Add(message);
    }
}
