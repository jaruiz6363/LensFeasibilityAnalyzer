using System;
using System.Collections.Generic;

namespace FeasibilityAnalyzer
{
    /// <summary>
    /// Derived quantities computed from input
    /// </summary>
    public class DerivedQuantities
    {
        public double NAObject { get; set; }
        public double FieldObject { get; set; }
        public double ThetaImage { get; set; }      // radians
        public double ThetaObject { get; set; }     // radians
        public double FieldAngleImage { get; set; } // radians
        public double FieldAngleObject { get; set; } // radians
        public double LagrangeInvariant { get; set; }
        public double FNumber { get; set; }
        public double ImageDistance { get; set; }
        public double WavelengthMM { get; set; }
        public double WDRatio { get; set; }
    }

    /// <summary>
    /// Sine condition vs distortion analysis result
    /// </summary>
    public class SineConditionResult
    {
        public bool Compatible { get; set; }
        public double MaxBetaDiscrepancy { get; set; }
        public double InducedComaWaves { get; set; }
        public List<FieldComaPoint> ComaVsField { get; set; } = new List<FieldComaPoint>();
        public string Recommendation { get; set; }
    }

    public class FieldComaPoint
    {
        public double NormalizedField { get; set; }
        public double ComaWaves { get; set; }
    }

    /// <summary>
    /// Reversibility analysis result (Walther's method)
    /// </summary>
    public class ReversibilityResult
    {
        public bool Constrained { get; set; }
        public double UnavoidableRMSWaves { get; set; }
        public double UnavoidablePVWaves { get; set; }
        public string DominantAberration { get; set; }
        public bool IsFeasible { get; set; }
        public string Recommendation { get; set; }
    }

    /// <summary>
    /// Working distance analysis at a single point
    /// </summary>
    public class WDPointResult
    {
        public double WorkingDistance { get; set; }
        public double Magnification { get; set; }
        public double NAObject { get; set; }
        public double FieldAngle { get; set; }
        public double RefocusNeeded { get; set; }
        public double SphericalScale { get; set; }
        public double ComaScale { get; set; }
        public double AstigmatismScale { get; set; }
    }

    /// <summary>
    /// Working distance range analysis result
    /// </summary>
    public class WorkingDistanceResult
    {
        public double WDRatio { get; set; }
        public WDPointResult AtMin { get; set; }
        public WDPointResult AtNominal { get; set; }
        public WDPointResult AtMax { get; set; }
        public double MaxDefocusWaves { get; set; }
        public bool NeedsRefocusCompensation { get; set; }
        public double DifficultyMetric { get; set; }
        public string Recommendation { get; set; }
    }

    /// <summary>
    /// Telecentricity analysis result
    /// </summary>
    public class TelecentricityResult
    {
        public bool ObjectTelecentric { get; set; }
        public bool ImageTelecentric { get; set; }
        public bool DoubleTelecentricFeasible { get; set; }
        public string ConstraintLevel { get; set; }
        public List<string> Conflicts { get; set; } = new List<string>();
        public List<string> Benefits { get; set; } = new List<string>();
        public string Recommendation { get; set; }
    }

    /// <summary>
    /// Petzval / field curvature analysis result
    /// </summary>
    public class PetzvalResult
    {
        public double PetzvalRadius { get; set; }
        public double PetzvalSag { get; set; }
        public double PetzvalDefocusWaves { get; set; }
        public bool FlatFieldRequired { get; set; }
        public string CorrectionDifficulty { get; set; }
        public List<string> CorrectionMethods { get; set; } = new List<string>();
        public string Warning { get; set; }
        public string Recommendation { get; set; }
    }

    /// <summary>
    /// Chromatic aberration analysis result
    /// </summary>
    public class ChromaticResult
    {
        public double AxialColorUncorrected { get; set; }
        public double LateralColorUncorrected { get; set; }
        public double AxialColorCorrected { get; set; }
        public double LateralColorCorrected { get; set; }
        public double SecondarySpectrum { get; set; }
        public int MinDoubletsRequired { get; set; }
        public bool NeedsApochromat { get; set; }
        public bool NeedsAnomalousDispersion { get; set; }
        public string CorrectionLevel { get; set; }
        public string Difficulty { get; set; }
        public List<string> RecommendedGlassPairs { get; set; } = new List<string>();
        public string Recommendation { get; set; }
    }

    /// <summary>
    /// Tolerance sensitivity analysis result
    /// </summary>
    public class ToleranceResult
    {
        public int ToleranceGrade { get; set; }
        public double ToleranceWavefrontRMS { get; set; }
        public double DesignMarginRequired { get; set; }
        
        // Sensitivities (waves per unit error)
        public double RadiusSensitivity { get; set; }      // waves per 0.1% radius error
        public double ThicknessSensitivity { get; set; }   // waves per 0.01mm thickness error
        public double DecenterSensitivity { get; set; }    // waves per 0.01mm decenter
        public double TiltSensitivity { get; set; }        // waves per arcmin tilt
        public double IndexSensitivity { get; set; }       // waves per 0.0001 index error
        public double AbbeSensitivity { get; set; }        // waves per 1% Abbe error
        public double WedgeSensitivity { get; set; }       // waves per arcmin wedge
        public double SurfaceIrregSensitivity { get; set; } // waves per 0.1 wave irregularity
        
        public string DominantSensitivity { get; set; }
        public double EstimatedYield { get; set; }
        public List<string> CriticalTolerances { get; set; } = new List<string>();
        public string Recommendation { get; set; }
        
        // Tolerance values for the selected grade
        public ToleranceValues ToleranceSpec { get; set; }
    }

    /// <summary>
    /// Tolerance specification values
    /// </summary>
    public class ToleranceValues
    {
        public double RadiusPercent { get; set; }
        public double ThicknessMM { get; set; }
        public double DecenterMM { get; set; }
        public double TiltArcmin { get; set; }
        public double IndexError { get; set; }
        public double AbbePercent { get; set; }
        public double WedgeArcmin { get; set; }
        public double SurfaceIrregWaves { get; set; }
        public double ScratchDig { get; set; }
    }

    /// <summary>
    /// Relative illumination analysis result
    /// </summary>
    public class RelativeIlluminationResult
    {
        public List<IlluminationFieldPoint> IlluminationVsField { get; set; } = new List<IlluminationFieldPoint>();
        public double MinRelativeIllumination { get; set; }  // Minimum RI across field (0-1)
        public double IlluminationUniformity { get; set; }   // 1 - (max-min)/(max+min)
        public double NaturalVignettingAtEdge { get; set; }  // cos^4 law contribution
        public double DistortionContribution { get; set; }   // Distortion effect on RI
        public double PupilAberrationContribution { get; set; } // Pupil aberration effect
        public string Severity { get; set; }  // "Negligible", "Minor", "Moderate", "Severe"
        public string Recommendation { get; set; }
    }

    public class IlluminationFieldPoint
    {
        public double NormalizedField { get; set; }
        public double RelativeIllumination { get; set; }  // 0-1, where 1 = on-axis
        public double NaturalVignetting { get; set; }     // cos^4 contribution
        public double DistortionEffect { get; set; }      // Magnification gradient effect
    }

    /// <summary>
    /// Differential distortion analysis result
    /// </summary>
    public class DifferentialDistortionResult
    {
        public List<DifferentialDistortionPoint> DifferentialVsField { get; set; } = new List<DifferentialDistortionPoint>();
        public double MaxDifferentialDistortion { get; set; }  // Maximum dβ/dθ deviation (%)
        public double TangentialRadialDifference { get; set; } // T-R difference at edge (%)
        public double LocalMagnificationVariation { get; set; } // Max local mag variation (%)
        public bool IsUniform { get; set; }  // True if differential distortion < 0.5%
        public string MetrologyImpact { get; set; }  // Impact on measurement applications
        public string Recommendation { get; set; }
    }

    public class DifferentialDistortionPoint
    {
        public double NormalizedField { get; set; }
        public double LocalMagnification { get; set; }      // β(θ) at this field
        public double MagnificationGradient { get; set; }   // dβ/dθ
        public double DifferentialDistortion { get; set; }  // Deviation from ideal (%)
        public double TangentialMagnification { get; set; } // Tangential direction
        public double RadialMagnification { get; set; }     // Radial direction
    }

    /// <summary>
    /// Illumination-distortion-differential coupling conflict analysis
    /// </summary>
    public class IlluminationDistortionConflictResult
    {
        public bool HasConflict { get; set; }
        public string ConflictSeverity { get; set; }  // "None", "Minor", "Moderate", "Severe"
        public double ConflictScore { get; set; }     // 0-100, higher = more conflict

        // Specific conflicts
        public bool IlluminationDistortionConflict { get; set; }
        public string IlluminationDistortionDetail { get; set; }

        public bool DifferentialSineConditionConflict { get; set; }
        public string DifferentialSineConditionDetail { get; set; }

        public bool IlluminationTelecentricityConflict { get; set; }
        public string IlluminationTelecentricityDetail { get; set; }

        public List<string> Conflicts { get; set; } = new List<string>();
        public List<string> Recommendations { get; set; } = new List<string>();
    }

    /// <summary>
    /// Aberration floor estimates
    /// </summary>
    public class AberrationFloors
    {
        public double Spherical { get; set; }
        public double Coma { get; set; }
        public double Astigmatism { get; set; }
        public double FieldCurvature { get; set; }
        public double Distortion { get; set; }
        public double ChromaticAxial { get; set; }
        public double ChromaticLateral { get; set; }
        public double TotalRMS { get; set; }
        public bool DiffractionLimited { get; set; }
        public string LimitingAberrations { get; set; }
    }

    /// <summary>
    /// Complexity estimate
    /// </summary>
    public class ComplexityEstimate
    {
        public int MinElements { get; set; }
        public int MaxElements { get; set; }
        public int EstimatedGroups { get; set; }
        public string Confidence { get; set; }
        public Dictionary<string, int> ElementBreakdown { get; set; } = new Dictionary<string, int>();
        public List<string> DominantDrivers { get; set; } = new List<string>();
        public List<string> SimilarDesigns { get; set; } = new List<string>();
        public int CostClass { get; set; }  // 1-5
    }

    /// <summary>
    /// Complete feasibility analysis output
    /// </summary>
    public class FeasibilityOutput
    {
        public DateTime AnalysisTime { get; set; } = DateTime.Now;
        public string Version { get; set; } = "1.1.0";
        
        // Overall assessment
        public bool IsFeasible { get; set; }
        public double FeasibilityScore { get; set; }  // 0-100
        public List<string> LimitingFactors { get; set; } = new List<string>();
        
        // Input and derived
        public FeasibilityInput Input { get; set; }
        public DerivedQuantities Derived { get; set; }
        
        // Analysis results
        public SineConditionResult SineCondition { get; set; }
        public ReversibilityResult Reversibility { get; set; }
        public WorkingDistanceResult WorkingDistance { get; set; }
        public TelecentricityResult Telecentricity { get; set; }
        public PetzvalResult Petzval { get; set; }
        public ChromaticResult Chromatic { get; set; }
        public ToleranceResult Tolerance { get; set; }
        public AberrationFloors AberrationFloors { get; set; }
        public ComplexityEstimate Complexity { get; set; }

        // New illumination and differential distortion analysis
        public RelativeIlluminationResult RelativeIllumination { get; set; }
        public DifferentialDistortionResult DifferentialDistortion { get; set; }
        public IlluminationDistortionConflictResult IlluminationDistortionConflict { get; set; }
        
        // Recommendations
        public List<string> Recommendations { get; set; } = new List<string>();
        public List<string> RelaxationSuggestions { get; set; } = new List<string>();
    }
}
