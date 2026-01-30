using System;
using System.Text;

namespace FeasibilityAnalyzer
{
    /// <summary>
    /// Generates formatted text reports from analysis results
    /// </summary>
    public static class ReportGenerator
    {
        private const double RAD_TO_DEG = 180.0 / Math.PI;
        
        /// <summary>
        /// Generate complete text report
        /// </summary>
        public static string GenerateReport(FeasibilityOutput output)
        {
            var sb = new StringBuilder();
            
            // Header
            sb.AppendLine(new string('═', 78));
            sb.AppendLine("                    LENS FEASIBILITY ANALYSIS REPORT");
            sb.AppendLine("                    Eikonal-Based Pre-Design Assessment");
            sb.AppendLine(new string('═', 78));
            sb.AppendLine();
            sb.AppendLine($"  Analysis Date:    {output.AnalysisTime:yyyy-MM-dd HH:mm:ss}");
            sb.AppendLine($"  Version:          {output.Version}");
            sb.AppendLine();
            
            // Overall assessment
            WriteSection(sb, "OVERALL ASSESSMENT");
            string status = output.IsFeasible ? "FEASIBLE" : "CHALLENGING";
            sb.AppendLine($"  Status:             {status}");
            sb.AppendLine($"  Feasibility Score:  {output.FeasibilityScore:F0} / 100");
            sb.AppendLine();
            
            if (output.LimitingFactors.Count > 0)
            {
                sb.AppendLine("  Limiting Factors:");
                foreach (var factor in output.LimitingFactors)
                {
                    sb.AppendLine($"    • {factor}");
                }
                sb.AppendLine();
            }
            
            // Input summary
            WriteSection(sb, "INPUT SPECIFICATION");
            WriteSubSection(sb, "Basic Optical Parameters");
            sb.AppendLine($"  Image NA:           {output.Input.NAImage:F3}");
            sb.AppendLine($"  Image semi-field:   {output.Input.FieldImage:F2} mm ({output.Input.FieldImage * 2:F2} mm diameter)");
            sb.AppendLine($"  Magnification:      {output.Input.Magnification:F2}×");
            sb.AppendLine();
            
            WriteSubSection(sb, "Working Distance");
            sb.AppendLine($"  Nominal:            {output.Input.WDNominal:F2} mm");
            sb.AppendLine($"  Range:              {output.Input.WDMin:F2} - {output.Input.WDMax:F2} mm");
            sb.AppendLine($"  Ratio:              {output.Derived.WDRatio:F2}×");
            sb.AppendLine();
            
            WriteSubSection(sb, "Constraints");
            sb.AppendLine($"  Distortion:         {output.Input.DistortionType} (max {output.Input.DistortionMax:F1}%)");
            sb.AppendLine($"  Reversible:         {(output.Input.IsReversible ? "Yes" : "No")}");
            sb.AppendLine($"  Object telecentric: {(output.Input.TelecentricObject ? "Yes" : "No")}");
            sb.AppendLine($"  Image telecentric:  {(output.Input.TelecentricImage ? "Yes" : "No")}");
            sb.AppendLine($"  Flat field:         {(output.Input.FlatField ? "Yes" : "No")}");
            sb.AppendLine();
            
            WriteSubSection(sb, "Wavelength");
            sb.AppendLine($"  Primary:            {output.Input.WavelengthPrimary:F1} nm");
            sb.AppendLine($"  Range:              {output.Input.WavelengthMin:F1} - {output.Input.WavelengthMax:F1} nm");
            sb.AppendLine($"  Bandwidth:          {output.Input.Bandwidth:F1} nm");
            sb.AppendLine();
            
            WriteSubSection(sb, "Tolerance & Performance");
            sb.AppendLine($"  Tolerance grade:    {output.Input.ToleranceGrade} ({GetToleranceGradeName(output.Input.ToleranceGrade)})");
            sb.AppendLine($"  Target RMS WFE:     {output.Input.TargetRMSWavefront:F3} waves");
            sb.AppendLine();
            
            // Derived quantities
            WriteSection(sb, "DERIVED QUANTITIES");
            sb.AppendLine($"  Object NA:          {output.Derived.NAObject:F4}");
            sb.AppendLine($"  Object semi-field:  {output.Derived.FieldObject:F3} mm");
            sb.AppendLine($"  Object field angle: {output.Derived.FieldAngleObject * RAD_TO_DEG:F2}°");
            sb.AppendLine($"  Image field angle:  {output.Derived.FieldAngleImage * RAD_TO_DEG:F2}°");
            sb.AppendLine($"  Lagrange invariant: {output.Derived.LagrangeInvariant:F4} mm");
            sb.AppendLine($"  F-number:           f/{output.Derived.FNumber:F2}");
            sb.AppendLine($"  Est. image distance:{output.Derived.ImageDistance:F1} mm");
            sb.AppendLine();
            
            // Sine condition analysis
            WriteSection(sb, "SINE CONDITION / DISTORTION ANALYSIS");
            sb.AppendLine($"  Compatible:         {(output.SineCondition.Compatible ? "Yes" : "No")}");
            sb.AppendLine($"  Max β discrepancy:  {output.SineCondition.MaxBetaDiscrepancy * 100:F2}%");
            sb.AppendLine($"  Induced coma:       {output.SineCondition.InducedComaWaves:F3} waves");
            sb.AppendLine();
            
            if (output.SineCondition.ComaVsField.Count > 0)
            {
                sb.AppendLine("  Coma vs Field:");
                sb.AppendLine("    Field    Coma (waves)");
                sb.AppendLine("    ─────    ────────────");
                foreach (var pt in output.SineCondition.ComaVsField)
                {
                    sb.AppendLine($"    {pt.NormalizedField:F1}      {pt.ComaWaves:F4}");
                }
                sb.AppendLine();
            }
            
            sb.AppendLine($"  » {output.SineCondition.Recommendation}");
            sb.AppendLine();
            
            // Reversibility analysis
            WriteSection(sb, "REVERSIBILITY ANALYSIS");
            if (output.Reversibility.Constrained)
            {
                sb.AppendLine($"  Constrained:        Yes (lens must work both directions)");
                sb.AppendLine($"  Unavoidable RMS:    {output.Reversibility.UnavoidableRMSWaves:F3} waves");
                sb.AppendLine($"  Unavoidable P-V:    {output.Reversibility.UnavoidablePVWaves:F3} waves");
                sb.AppendLine($"  Dominant aberr:     {output.Reversibility.DominantAberration}");
                sb.AppendLine($"  DL achievable:      {(output.Reversibility.IsFeasible ? "Yes" : "No")}");
            }
            else
            {
                sb.AppendLine($"  Constrained:        No");
            }
            sb.AppendLine();
            sb.AppendLine($"  » {output.Reversibility.Recommendation}");
            sb.AppendLine();
            
            // Working distance analysis
            WriteSection(sb, "WORKING DISTANCE RANGE ANALYSIS");
            sb.AppendLine($"  WD range ratio:     {output.WorkingDistance.WDRatio:F2}×");
            sb.AppendLine($"  Max defocus:        {output.WorkingDistance.MaxDefocusWaves:F2} waves (if not refocused)");
            sb.AppendLine($"  Refocus needed:     {(output.WorkingDistance.NeedsRefocusCompensation ? "Yes" : "No")}");
            sb.AppendLine($"  Difficulty metric:  {output.WorkingDistance.DifficultyMetric:F3}");
            sb.AppendLine();
            
            WriteSubSection(sb, "At Minimum WD");
            WriteWDPoint(sb, output.WorkingDistance.AtMin);
            
            WriteSubSection(sb, "At Nominal WD");
            WriteWDPoint(sb, output.WorkingDistance.AtNominal);
            
            WriteSubSection(sb, "At Maximum WD");
            WriteWDPoint(sb, output.WorkingDistance.AtMax);
            
            sb.AppendLine($"  » {output.WorkingDistance.Recommendation}");
            sb.AppendLine();
            
            // Telecentricity analysis
            WriteSection(sb, "TELECENTRICITY ANALYSIS");
            sb.AppendLine($"  Object telecentric: {(output.Telecentricity.ObjectTelecentric ? "Required" : "Not required")}");
            sb.AppendLine($"  Image telecentric:  {(output.Telecentricity.ImageTelecentric ? "Required" : "Not required")}");
            sb.AppendLine($"  Constraint level:   {output.Telecentricity.ConstraintLevel}");
            
            if (output.Telecentricity.ObjectTelecentric && output.Telecentricity.ImageTelecentric)
            {
                sb.AppendLine($"  Double telecentric: {(output.Telecentricity.DoubleTelecentricFeasible ? "Feasible" : "Not achievable")}");
            }
            
            if (output.Telecentricity.Benefits.Count > 0)
            {
                sb.AppendLine("  Benefits:");
                foreach (var b in output.Telecentricity.Benefits)
                    sb.AppendLine($"    + {b}");
            }
            
            if (output.Telecentricity.Conflicts.Count > 0)
            {
                sb.AppendLine("  Conflicts:");
                foreach (var c in output.Telecentricity.Conflicts)
                    sb.AppendLine($"    - {c}");
            }
            sb.AppendLine();
            sb.AppendLine($"  » {output.Telecentricity.Recommendation}");
            sb.AppendLine();
            
            // Petzval analysis
            WriteSection(sb, "FIELD CURVATURE / PETZVAL ANALYSIS");
            sb.AppendLine($"  Flat field required:{(output.Petzval.FlatFieldRequired ? "Yes" : "No")}");
            sb.AppendLine($"  Est. Petzval radius:{output.Petzval.PetzvalRadius:F1} mm");
            sb.AppendLine($"  Petzval sag:        {output.Petzval.PetzvalSag:F3} mm at edge");
            sb.AppendLine($"  Petzval defocus:    {output.Petzval.PetzvalDefocusWaves:F2} waves (uncorrected)");
            sb.AppendLine($"  Correction level:   {output.Petzval.CorrectionDifficulty}");
            
            if (!string.IsNullOrEmpty(output.Petzval.Warning))
            {
                sb.AppendLine($"  ⚠ {output.Petzval.Warning}");
            }
            
            if (output.Petzval.CorrectionMethods.Count > 0)
            {
                sb.AppendLine("  Correction methods:");
                foreach (var m in output.Petzval.CorrectionMethods)
                    sb.AppendLine($"    • {m}");
            }
            sb.AppendLine();
            sb.AppendLine($"  » {output.Petzval.Recommendation}");
            sb.AppendLine();
            
            // Chromatic analysis
            WriteSection(sb, "CHROMATIC ABERRATION ANALYSIS");
            sb.AppendLine($"  Bandwidth:          {output.Input.Bandwidth:F0} nm");
            sb.AppendLine($"  Correction level:   {output.Chromatic.CorrectionLevel}");
            sb.AppendLine($"  Difficulty:         {output.Chromatic.Difficulty}");
            sb.AppendLine();
            sb.AppendLine($"  Axial color (uncorrected):  {output.Chromatic.AxialColorUncorrected:F2} waves");
            sb.AppendLine($"  Axial color (corrected):    {output.Chromatic.AxialColorCorrected:F3} waves");
            sb.AppendLine($"  Lateral color (uncorrected):{output.Chromatic.LateralColorUncorrected:F2} waves");
            sb.AppendLine($"  Lateral color (corrected):  {output.Chromatic.LateralColorCorrected:F3} waves");
            sb.AppendLine($"  Secondary spectrum:         {output.Chromatic.SecondarySpectrum:F4} waves");
            sb.AppendLine();
            sb.AppendLine($"  Min doublets required:      {output.Chromatic.MinDoubletsRequired}");
            sb.AppendLine($"  Needs apochromat:           {(output.Chromatic.NeedsApochromat ? "Yes" : "No")}");
            sb.AppendLine($"  Needs anomalous glass:      {(output.Chromatic.NeedsAnomalousDispersion ? "Yes" : "No")}");
            
            if (output.Chromatic.RecommendedGlassPairs.Count > 0)
            {
                sb.AppendLine("  Recommended glass pairs:");
                foreach (var g in output.Chromatic.RecommendedGlassPairs)
                    sb.AppendLine($"    • {g}");
            }
            sb.AppendLine();
            sb.AppendLine($"  » {output.Chromatic.Recommendation}");
            sb.AppendLine();
            
            // Tolerance analysis
            WriteSection(sb, "TOLERANCE SENSITIVITY ANALYSIS");
            sb.AppendLine($"  Tolerance grade:    {output.Tolerance.ToleranceGrade} ({GetToleranceGradeName(output.Tolerance.ToleranceGrade)})");
            sb.AppendLine($"  Est. tolerance RMS: {output.Tolerance.ToleranceWavefrontRMS:F3} waves");
            sb.AppendLine($"  Design margin need: {output.Tolerance.DesignMarginRequired:F3} waves");
            sb.AppendLine($"  Estimated yield:    {output.Tolerance.EstimatedYield:P0}");
            sb.AppendLine($"  Dominant sensitivity: {output.Tolerance.DominantSensitivity}");
            sb.AppendLine();
            
            WriteSubSection(sb, "Tolerance Specification (Grade " + output.Tolerance.ToleranceGrade + ")");
            var spec = output.Tolerance.ToleranceSpec;
            sb.AppendLine($"    Radius:           ±{spec.RadiusPercent:F2}%");
            sb.AppendLine($"    Thickness:        ±{spec.ThicknessMM:F3} mm");
            sb.AppendLine($"    Decenter:         ±{spec.DecenterMM:F3} mm");
            sb.AppendLine($"    Tilt:             ±{spec.TiltArcmin:F1} arcmin");
            sb.AppendLine($"    Index (Δn):       ±{spec.IndexError:F5}");
            sb.AppendLine($"    Abbe (ΔV/V):      ±{spec.AbbePercent:F1}%");
            sb.AppendLine($"    Wedge:            ±{spec.WedgeArcmin:F1} arcmin");
            sb.AppendLine($"    Surface irreg:    ±{spec.SurfaceIrregWaves:F2} waves");
            sb.AppendLine();
            
            WriteSubSection(sb, "Sensitivity Breakdown (waves per unit error)");
            sb.AppendLine($"    Radius:           {output.Tolerance.RadiusSensitivity:F3} per 0.1%");
            sb.AppendLine($"    Thickness:        {output.Tolerance.ThicknessSensitivity:F3} per 0.01mm");
            sb.AppendLine($"    Decenter:         {output.Tolerance.DecenterSensitivity:F3} per 0.01mm");
            sb.AppendLine($"    Tilt:             {output.Tolerance.TiltSensitivity:F3} per arcmin");
            sb.AppendLine($"    Index:            {output.Tolerance.IndexSensitivity:F3} per 0.0001");
            sb.AppendLine($"    Surface irreg:    {output.Tolerance.SurfaceIrregSensitivity:F3} per 0.1 wave");
            sb.AppendLine();
            
            if (output.Tolerance.CriticalTolerances.Count > 0)
            {
                sb.AppendLine("  Critical tolerances:");
                foreach (var t in output.Tolerance.CriticalTolerances)
                    sb.AppendLine($"    ⚠ {t}");
                sb.AppendLine();
            }
            
            sb.AppendLine($"  » {output.Tolerance.Recommendation}");
            sb.AppendLine();
            
            // Aberration floors
            WriteSection(sb, "THEORETICAL ABERRATION FLOORS");
            sb.AppendLine("  These are MINIMUM achievable aberrations - actual design may be worse.");
            sb.AppendLine();
            sb.AppendLine($"  Spherical:          {output.AberrationFloors.Spherical:F3} waves");
            sb.AppendLine($"  Coma:               {output.AberrationFloors.Coma:F3} waves");
            sb.AppendLine($"  Astigmatism:        {output.AberrationFloors.Astigmatism:F3} waves");
            sb.AppendLine($"  Field curvature:    {output.AberrationFloors.FieldCurvature:F3} waves");
            sb.AppendLine($"  Distortion:         {output.AberrationFloors.Distortion:F3}%");
            sb.AppendLine($"  Axial chromatic:    {output.AberrationFloors.ChromaticAxial:F3} waves");
            sb.AppendLine($"  Lateral chromatic:  {output.AberrationFloors.ChromaticLateral:F3} waves");
            sb.AppendLine("  ────────────────────────────────");
            sb.AppendLine($"  Total RMS floor:    {output.AberrationFloors.TotalRMS:F3} waves");
            sb.AppendLine();
            sb.AppendLine($"  Diffraction limited:{(output.AberrationFloors.DiffractionLimited ? "Achievable" : "Not achievable")}");
            sb.AppendLine($"  Limiting factors:   {output.AberrationFloors.LimitingAberrations}");
            sb.AppendLine();
            
            // Complexity estimate
            WriteSection(sb, "COMPLEXITY ESTIMATE");
            sb.AppendLine($"  Estimated elements: {output.Complexity.MinElements} - {output.Complexity.MaxElements}");
            sb.AppendLine($"  Estimated groups:   {output.Complexity.EstimatedGroups}");
            sb.AppendLine($"  Confidence:         {output.Complexity.Confidence}");
            sb.AppendLine($"  Cost class:         {output.Complexity.CostClass}/5 ({GetCostClassName(output.Complexity.CostClass)})");
            sb.AppendLine();
            
            if (output.Complexity.ElementBreakdown.Count > 0)
            {
                sb.AppendLine("  Element breakdown:");
                foreach (var kv in output.Complexity.ElementBreakdown)
                    sb.AppendLine($"    {kv.Key}: +{kv.Value}");
                sb.AppendLine();
            }
            
            if (output.Complexity.DominantDrivers.Count > 0)
            {
                sb.AppendLine($"  Dominant drivers:   {string.Join(", ", output.Complexity.DominantDrivers)}");
            }
            
            if (output.Complexity.SimilarDesigns.Count > 0)
            {
                sb.AppendLine($"  Similar designs:    {string.Join(", ", output.Complexity.SimilarDesigns)}");
            }
            sb.AppendLine();
            
            // Recommendations
            WriteSection(sb, "RECOMMENDATIONS");
            for (int i = 0; i < output.Recommendations.Count; i++)
            {
                sb.AppendLine($"  {i + 1}. {output.Recommendations[i]}");
            }
            sb.AppendLine();
            
            if (output.RelaxationSuggestions.Count > 0)
            {
                WriteSubSection(sb, "If Design Proves Difficult, Consider");
                foreach (var s in output.RelaxationSuggestions)
                {
                    sb.AppendLine($"    • {s}");
                }
                sb.AppendLine();
            }
            
            // Footer
            sb.AppendLine(new string('═', 78));
            sb.AppendLine("                          END OF FEASIBILITY REPORT");
            sb.AppendLine(new string('═', 78));
            
            return sb.ToString();
        }
        
        private static void WriteSection(StringBuilder sb, string title)
        {
            sb.AppendLine(new string('─', 78));
            sb.AppendLine($"  {title}");
            sb.AppendLine(new string('─', 78));
        }
        
        private static void WriteSubSection(StringBuilder sb, string title)
        {
            sb.AppendLine($"  [{title}]");
        }
        
        private static void WriteWDPoint(StringBuilder sb, WDPointResult pt)
        {
            sb.AppendLine($"    WD = {pt.WorkingDistance:F2} mm:");
            sb.AppendLine($"      Magnification:    {pt.Magnification:F3}×");
            sb.AppendLine($"      Object NA:        {pt.NAObject:F4}");
            sb.AppendLine($"      Field angle:      {pt.FieldAngle * RAD_TO_DEG:F2}°");
            sb.AppendLine($"      Refocus needed:   {pt.RefocusNeeded:F2} mm");
            sb.AppendLine($"      Aberr. scaling:   SA×{pt.SphericalScale:F2}, Coma×{pt.ComaScale:F2}, Astig×{pt.AstigmatismScale:F2}");
            sb.AppendLine();
        }
        
        private static string GetToleranceGradeName(int grade)
        {
            return grade switch
            {
                1 => "Commercial",
                2 => "Precision",
                3 => "High Precision",
                _ => "Unknown"
            };
        }
        
        private static string GetCostClassName(int costClass)
        {
            return costClass switch
            {
                1 => "Low",
                2 => "Moderate",
                3 => "Medium-High",
                4 => "High",
                5 => "Very High",
                _ => "Unknown"
            };
        }
        
        /// <summary>
        /// Generate a brief summary (single paragraph)
        /// </summary>
        public static string GenerateSummary(FeasibilityOutput output)
        {
            var sb = new StringBuilder();
            
            string status = output.IsFeasible ? "appears feasible" : "may be challenging";
            sb.Append($"This {output.Input.NAImage:F2} NA, {output.Input.FieldImage * 2:F1}mm field, ");
            sb.Append($"{Math.Abs(output.Input.Magnification):F1}× magnification lens {status} ");
            sb.Append($"(score: {output.FeasibilityScore:F0}/100). ");
            
            sb.Append($"Estimated {output.Complexity.MinElements}-{output.Complexity.MaxElements} elements. ");
            
            if (output.AberrationFloors.DiffractionLimited)
            {
                sb.Append("Diffraction-limited performance is theoretically achievable. ");
            }
            else
            {
                sb.Append($"Aberration floor: {output.AberrationFloors.TotalRMS:F3}λ RMS ");
                sb.Append($"(limited by {output.AberrationFloors.LimitingAberrations}). ");
            }
            
            if (output.LimitingFactors.Count > 0)
            {
                sb.Append($"Key challenges: {string.Join("; ", output.LimitingFactors)}.");
            }
            
            return sb.ToString();
        }
    }
}
