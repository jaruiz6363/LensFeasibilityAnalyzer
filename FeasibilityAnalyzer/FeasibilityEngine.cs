using System;
using System.Collections.Generic;

namespace FeasibilityAnalyzer
{
    /// <summary>
    /// Main feasibility analysis engine
    /// Coordinates all analysis modules and produces comprehensive report
    /// </summary>
    public class FeasibilityEngine
    {
        private const double DEG_TO_RAD = Math.PI / 180.0;
        private const double RAD_TO_DEG = 180.0 / Math.PI;
        
        /// <summary>
        /// Perform complete feasibility analysis
        /// </summary>
        public FeasibilityOutput Analyze(FeasibilityInput input)
        {
            // Validate input
            var validation = input.Validate();
            if (!validation.IsValid)
            {
                throw new ArgumentException(
                    "Invalid input: " + string.Join("; ", validation.Errors));
            }
            
            var output = new FeasibilityOutput { Input = input };
            
            // Step 1: Compute derived quantities
            output.Derived = ComputeDerivedQuantities(input);
            
            // Step 2: Run all analysis modules
            output.SineCondition = AnalyzeSineCondition(input, output.Derived);
            output.Reversibility = AnalyzeReversibility(input, output.Derived);
            output.WorkingDistance = AnalyzeWorkingDistance(input, output.Derived);
            output.Telecentricity = AnalyzeTelecentricity(input, output.Derived);
            output.Petzval = AnalyzePetzval(input, output.Derived);
            output.Chromatic = AnalyzeChromatic(input, output.Derived);
            output.Tolerance = AnalyzeTolerance(input, output.Derived);
            
            // Step 3: Analyze relative illumination and differential distortion
            output.RelativeIllumination = AnalyzeRelativeIllumination(input, output.Derived);
            output.DifferentialDistortion = AnalyzeDifferentialDistortion(input, output.Derived);

            // Step 4: Check for conflicts between illumination, distortion, and differential distortion
            output.IlluminationDistortionConflict = CheckIlluminationDistortionConflict(
                input, output.Derived, output.RelativeIllumination, output.DifferentialDistortion,
                output.SineCondition, output.Telecentricity);

            // Step 5: Compute aberration floors
            output.AberrationFloors = ComputeAberrationFloors(
                input, output.Derived, output.SineCondition,
                output.Reversibility, output.Chromatic, output.Petzval);

            // Step 7: Estimate complexity
            output.Complexity = EstimateComplexity(input, output);

            // Step 8: Compute overall feasibility score
            ComputeFeasibilityScore(output);

            // Step 9: Generate recommendations
            GenerateRecommendations(output);
            
            return output;
        }
        
        /// <summary>
        /// Compute derived quantities from input
        /// </summary>
        private DerivedQuantities ComputeDerivedQuantities(FeasibilityInput input)
        {
            var d = new DerivedQuantities();
            
            double absMag = Math.Abs(input.Magnification);
            
            // Object-side NA
            d.NAObject = input.NAImage / absMag;
            
            // Object semi-field
            d.FieldObject = input.FieldImage / absMag;
            
            // Ray angles
            d.ThetaImage = Math.Asin(Math.Min(input.NAImage, 0.9999));
            d.ThetaObject = Math.Asin(Math.Min(d.NAObject, 0.9999));
            
            // Image distance estimate
            d.ImageDistance = absMag * input.WDNominal;
            
            // Field angles
            d.FieldAngleObject = Math.Atan(d.FieldObject / input.WDNominal);
            d.FieldAngleImage = Math.Atan(input.FieldImage / d.ImageDistance);
            
            // Lagrange invariant
            d.LagrangeInvariant = input.NAImage * input.FieldImage;
            
            // F-number
            d.FNumber = 1.0 / (2.0 * input.NAImage);
            
            // Wavelength in mm
            d.WavelengthMM = input.WavelengthPrimary * 1e-6;
            
            // WD ratio
            d.WDRatio = input.WDMax / input.WDMin;
            
            return d;
        }
        
        /// <summary>
        /// Analyze sine condition vs distortion compatibility
        /// </summary>
        private SineConditionResult AnalyzeSineCondition(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new SineConditionResult();
            
            double[] fieldFracs = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
            double maxDiscrepancy = 0;
            double maxComa = 0;
            double beta0 = input.Magnification;
            
            foreach (double frac in fieldFracs)
            {
                double theta = frac * derived.FieldAngleObject;
                if (theta < 1e-10) continue;
                
                // Sine condition magnification
                double betaSine = beta0 / Math.Cos(theta);
                
                // Distortion-required magnification
                double betaDist = ComputeDistortionMagnification(
                    theta, input.DistortionType, beta0, 
                    input.DistortionMax / 100.0, input.WDNominal);
                
                double deltaBeta = (betaDist - betaSine) / Math.Abs(betaSine);
                
                // Coma estimate (waves)
                double comaWaves = Math.Abs(deltaBeta) * 
                    Math.Pow(derived.NAObject, 2) * 
                    Math.Sin(theta) / derived.WavelengthMM * 0.5;
                
                result.ComaVsField.Add(new FieldComaPoint 
                { 
                    NormalizedField = frac, 
                    ComaWaves = comaWaves 
                });
                
                maxDiscrepancy = Math.Max(maxDiscrepancy, Math.Abs(deltaBeta));
                maxComa = Math.Max(maxComa, comaWaves);
            }
            
            result.MaxBetaDiscrepancy = maxDiscrepancy;
            result.InducedComaWaves = maxComa;
            result.Compatible = maxDiscrepancy < 0.01;
            
            if (maxDiscrepancy > 0.05)
            {
                result.Recommendation = 
                    "Significant conflict between distortion and coma correction. " +
                    $"Induced coma: {maxComa:F2} waves. Consider relaxing distortion " +
                    "requirement or accepting residual coma at edge of field.";
            }
            else if (maxDiscrepancy > 0.01)
            {
                result.Recommendation = 
                    "Minor sine condition conflict. Optimization can likely balance " +
                    $"distortion and coma. Expected residual: {maxComa:F3} waves.";
            }
            else
            {
                result.Recommendation = 
                    "Distortion and coma correction are compatible. No fundamental " +
                    "conflict from sine condition.";
            }
            
            return result;
        }
        
        /// <summary>
        /// Compute magnification required by distortion type
        /// </summary>
        private double ComputeDistortionMagnification(
            double theta, DistortionType type, double beta0, 
            double distMax, double WD)
        {
            switch (type)
            {
                case DistortionType.None:
                    return beta0;
                    
                case DistortionType.FTheta:
                    if (theta < 1e-10) return beta0;
                    return beta0 * theta / Math.Tan(theta);
                    
                case DistortionType.Orthographic:
                    return beta0 * Math.Cos(theta);
                    
                case DistortionType.Barrel:
                    return beta0 * (1.0 - distMax * Math.Pow(Math.Tan(theta), 2));
                    
                case DistortionType.Pincushion:
                    return beta0 * (1.0 + distMax * Math.Pow(Math.Tan(theta), 2));
                    
                default:
                    return beta0;
            }
        }
        
        /// <summary>
        /// Analyze reversibility limits (Walther's method)
        /// </summary>
        private ReversibilityResult AnalyzeReversibility(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new ReversibilityResult();
            result.Constrained = input.IsReversible;
            
            if (!input.IsReversible)
            {
                result.UnavoidableRMSWaves = 0;
                result.UnavoidablePVWaves = 0;
                result.IsFeasible = true;
                result.Recommendation = "No reversibility constraint.";
                return result;
            }
            
            // For reversible lens, S(u,v,w) must equal S(w,v,u)
            // The antisymmetric part is unavoidable aberration
            
            double t = derived.ImageDistance;
            double G = input.Magnification;
            
            double[] samples = { 0.0, 0.25, 0.5, 0.75, 1.0 };
            double maxAberr = 0;
            double sumSq = 0;
            int n = 0;
            
            foreach (double u in samples)
            {
                foreach (double w in samples)
                {
                    double vMax = Math.Sqrt(u * w);
                    double[] vSamples = vMax > 0.01 ? 
                        new[] { -vMax, 0, vMax } : new[] { 0.0 };
                    
                    foreach (double v in vSamples)
                    {
                        double S_uvw = ComputeIdealEikonal(u, v, w, t, G);
                        double S_wvu = ComputeIdealEikonal(w, v, u, t, G);
                        
                        double W_anti = 0.5 * (S_uvw - S_wvu);
                        double W_waves = W_anti / derived.WavelengthMM;
                        
                        maxAberr = Math.Max(maxAberr, Math.Abs(W_waves));
                        sumSq += W_waves * W_waves;
                        n++;
                    }
                }
            }
            
            result.UnavoidablePVWaves = maxAberr;
            result.UnavoidableRMSWaves = Math.Sqrt(sumSq / n);
            result.IsFeasible = result.UnavoidableRMSWaves < 0.07;
            
            // Classify dominant aberration
            if (Math.Abs(input.Magnification) > 1.5)
                result.DominantAberration = "Coma and field-dependent aberrations";
            else if (Math.Abs(input.Magnification) < 0.67)
                result.DominantAberration = "Spherical-like aberration";
            else
                result.DominantAberration = "Mixed aberrations";
            
            if (result.UnavoidableRMSWaves > 0.25)
            {
                result.Recommendation = 
                    $"Reversibility causes {result.UnavoidableRMSWaves:F2}λ RMS " +
                    "unavoidable aberration. Consider asymmetric design.";
            }
            else if (result.UnavoidableRMSWaves > 0.07)
            {
                result.Recommendation = 
                    $"Reversibility limits performance to {result.UnavoidableRMSWaves:F2}λ RMS. " +
                    "Diffraction limit not achievable at full field.";
            }
            else
            {
                result.Recommendation = 
                    $"Reversibility constraint acceptable ({result.UnavoidableRMSWaves:F3}λ RMS).";
            }
            
            return result;
        }
        
        /// <summary>
        /// Compute ideal point eikonal
        /// </summary>
        private double ComputeIdealEikonal(double u, double v, double w, double t, double G)
        {
            double arg = 1.0 + G * G * u - 2.0 * G * v + w;
            if (arg < 0) arg = 1e-10;
            
            double A = t * Math.Sqrt(arg);
            double A_chief = t * Math.Sqrt(1.0 + G * G * u);
            double F = t - A_chief;
            
            return F + A;
        }
        
        /// <summary>
        /// Analyze working distance range
        /// </summary>
        private WorkingDistanceResult AnalyzeWorkingDistance(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new WorkingDistanceResult();
            result.WDRatio = derived.WDRatio;
            
            // Analyze at three WD points
            result.AtMin = AnalyzeAtWD(input, derived, input.WDMin);
            result.AtNominal = AnalyzeAtWD(input, derived, input.WDNominal);
            result.AtMax = AnalyzeAtWD(input, derived, input.WDMax);
            
            // Maximum refocus needed
            double maxRefocus = Math.Max(
                Math.Abs(result.AtMin.RefocusNeeded),
                Math.Abs(result.AtMax.RefocusNeeded));
            
            // Defocus aberration if not refocused
            result.MaxDefocusWaves = maxRefocus * 
                Math.Pow(input.NAImage, 2) / (2.0 * derived.WavelengthMM);
            
            result.NeedsRefocusCompensation = result.MaxDefocusWaves > 0.25;
            
            // Difficulty metric
            result.DifficultyMetric = Math.Log(derived.WDRatio) * 
                Math.Pow(derived.NAObject, 2) * derived.FieldAngleObject;
            
            // Recommendation
            if (derived.WDRatio > 3 && input.NAImage > 0.3)
            {
                result.Recommendation = 
                    "Large WD range with high NA is very challenging. " +
                    "Consider variable focus element or reduced WD range.";
            }
            else if (derived.WDRatio > 2)
            {
                result.Recommendation = 
                    "Moderate WD range. Design must balance aberrations across range. " +
                    "Optimize at nominal WD, verify at extremes.";
            }
            else if (result.NeedsRefocusCompensation)
            {
                result.Recommendation = 
                    "WD range requires refocusing. Consider motorized focus or " +
                    "extended depth of focus design.";
            }
            else
            {
                result.Recommendation = 
                    "WD range is manageable with standard design approach.";
            }
            
            return result;
        }
        
        private WDPointResult AnalyzeAtWD(
            FeasibilityInput input, DerivedQuantities derived, double WD)
        {
            var point = new WDPointResult();
            point.WorkingDistance = WD;
            
            // Magnification if image distance fixed
            point.Magnification = input.Magnification * input.WDNominal / WD;
            point.NAObject = input.NAImage / Math.Abs(point.Magnification);
            point.FieldAngle = Math.Atan(derived.FieldObject / WD);
            
            // Refocus needed to maintain magnification
            point.RefocusNeeded = (WD - input.WDNominal) * 
                Math.Pow(input.Magnification, 2);
            
            // Aberration scaling
            double naRatio = point.NAObject / derived.NAObject;
            double angleRatio = point.FieldAngle / derived.FieldAngleObject;
            
            point.SphericalScale = Math.Pow(naRatio, 4);
            point.ComaScale = Math.Pow(naRatio, 3) * angleRatio;
            point.AstigmatismScale = Math.Pow(naRatio, 2) * Math.Pow(angleRatio, 2);
            
            return point;
        }
        
        /// <summary>
        /// Analyze telecentricity requirements
        /// </summary>
        private TelecentricityResult AnalyzeTelecentricity(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new TelecentricityResult();
            result.ObjectTelecentric = input.TelecentricObject;
            result.ImageTelecentric = input.TelecentricImage;
            
            // Analyze conflicts and benefits
            if (input.TelecentricObject)
            {
                result.Benefits.Add("Constant object-side NA across field");
                result.Benefits.Add("Reduced perspective distortion");
                
                if (input.DistortionType == DistortionType.FTheta)
                    result.Conflicts.Add("Object telecentricity conflicts with f-theta");
            }
            
            if (input.TelecentricImage)
            {
                result.Benefits.Add("Constant image-side NA across field");
                result.Benefits.Add("Optimal for flat sensors");
                result.Benefits.Add("Reduced vignetting at sensor");
            }
            
            // Double telecentricity analysis
            if (input.TelecentricObject && input.TelecentricImage)
            {
                if (Math.Abs(Math.Abs(input.Magnification) - 1.0) < 0.01)
                {
                    result.DoubleTelecentricFeasible = true;
                    result.ConstraintLevel = "Very High";
                    result.Recommendation = 
                        "Double telecentric 1:1 relay is feasible but highly constrained. " +
                        "Consider symmetric design.";
                }
                else
                {
                    result.DoubleTelecentricFeasible = false;
                    result.ConstraintLevel = "Not Achievable";
                    result.Conflicts.Add(
                        "True double telecentricity requires |m|=1. " +
                        "Consider quasi-telecentric design.");
                    result.Recommendation = 
                        "Double telecentricity not achievable with m≠1. " +
                        "Use quasi-telecentric (small chief ray angle) or " +
                        "relax one telecentricity requirement.";
                }
            }
            else if (input.TelecentricObject || input.TelecentricImage)
            {
                result.ConstraintLevel = "Moderate";
                result.Recommendation = 
                    "Single-side telecentricity is achievable. Adds complexity but " +
                    "well-understood design approaches exist.";
            }
            else
            {
                result.ConstraintLevel = "Low";
                result.Recommendation = "No telecentricity constraint.";
            }
            
            return result;
        }
        
        /// <summary>
        /// Analyze Petzval / field curvature
        /// </summary>
        private PetzvalResult AnalyzePetzval(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new PetzvalResult();
            result.FlatFieldRequired = input.FlatField;
            
            // Estimate Petzval radius for a simple lens
            // Petzval radius ≈ -n*f where n≈1.5 for typical glass
            double f_eff = derived.ImageDistance / Math.Abs(input.Magnification);
            result.PetzvalRadius = -1.5 * f_eff;
            
            // Petzval sag at edge of field
            result.PetzvalSag = Math.Pow(input.FieldImage, 2) / 
                (2.0 * Math.Abs(result.PetzvalRadius));
            
            // Defocus in waves
            result.PetzvalDefocusWaves = result.PetzvalSag * 
                Math.Pow(input.NAImage, 2) / (2.0 * derived.WavelengthMM);
            
            // Correction methods
            result.CorrectionMethods.Add("Field flattener near image plane");
            result.CorrectionMethods.Add("Thick meniscus elements");
            result.CorrectionMethods.Add("Negative power element near stop");
            result.CorrectionMethods.Add("Astigmatism balance (trades T/S quality)");
            
            // Difficulty assessment
            if (result.PetzvalDefocusWaves > 5)
            {
                result.CorrectionDifficulty = "Very High";
                result.Warning = "Large Petzval sum makes flat field very difficult.";
            }
            else if (result.PetzvalDefocusWaves > 2)
            {
                result.CorrectionDifficulty = "High";
            }
            else if (result.PetzvalDefocusWaves > 0.5)
            {
                result.CorrectionDifficulty = "Moderate";
            }
            else
            {
                result.CorrectionDifficulty = "Low";
            }
            
            if (input.FlatField)
            {
                if (input.NAImage > 0.4 && input.FieldImage > 10)
                {
                    result.Warning = 
                        "High NA + large field + flat field is a demanding combination.";
                }
                
                result.Recommendation = 
                    $"Flat field correction needed ({result.CorrectionDifficulty} difficulty). " +
                    $"Petzval defocus: {result.PetzvalDefocusWaves:F1} waves uncorrected.";
            }
            else
            {
                result.Recommendation = 
                    "Curved field allowed - no Petzval correction required.";
            }
            
            return result;
        }
        
        /// <summary>
        /// Analyze chromatic aberration
        /// </summary>
        private ChromaticResult AnalyzeChromatic(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new ChromaticResult();
            
            // Estimate uncorrected chromatic aberration
            // Axial color ≈ f / V where V is Abbe number (≈60 for typical glass)
            double f_eff = derived.ImageDistance / Math.Abs(input.Magnification);
            double V_typical = 60.0;
            
            double axialColorMM = f_eff / V_typical * 
                (input.Bandwidth / input.WavelengthPrimary);
            result.AxialColorUncorrected = axialColorMM * 
                Math.Pow(input.NAImage, 2) / (2.0 * derived.WavelengthMM);
            
            // Lateral color ≈ H / V where H is Lagrange invariant
            double latColorMM = derived.LagrangeInvariant / V_typical *
                (input.Bandwidth / input.WavelengthPrimary);
            result.LateralColorUncorrected = latColorMM / derived.WavelengthMM * 0.1;
            
            // Determine correction level needed
            double bandwidth = input.Bandwidth;
            bool needsApochromat = false;
            bool needsAnomalous = false;
            
            if (bandwidth > 300)
            {
                result.CorrectionLevel = "Superachromatic";
                needsApochromat = true;
                needsAnomalous = true;
                result.MinDoubletsRequired = 3;
            }
            else if (bandwidth > 150)
            {
                result.CorrectionLevel = "Apochromatic";
                needsApochromat = true;
                result.MinDoubletsRequired = 2;
            }
            else if (bandwidth > 50)
            {
                result.CorrectionLevel = "Achromatic";
                result.MinDoubletsRequired = 1;
            }
            else
            {
                result.CorrectionLevel = "Near-monochromatic";
                result.MinDoubletsRequired = 0;
            }
            
            result.NeedsApochromat = needsApochromat;
            result.NeedsAnomalousDispersion = needsAnomalous;
            
            // Corrected values (estimates)
            if (needsApochromat)
            {
                result.AxialColorCorrected = result.AxialColorUncorrected * 0.02;
                result.LateralColorCorrected = result.LateralColorUncorrected * 0.05;
                result.SecondarySpectrum = result.AxialColorUncorrected * 0.01;
            }
            else if (result.MinDoubletsRequired > 0)
            {
                result.AxialColorCorrected = result.AxialColorUncorrected * 0.05;
                result.LateralColorCorrected = result.LateralColorUncorrected * 0.1;
                result.SecondarySpectrum = result.AxialColorUncorrected * 0.02;
            }
            else
            {
                result.AxialColorCorrected = result.AxialColorUncorrected;
                result.LateralColorCorrected = result.LateralColorUncorrected;
                result.SecondarySpectrum = 0;
            }
            
            // Recommended glass pairs
            if (needsAnomalous)
            {
                result.RecommendedGlassPairs.Add("S-FPL51 + S-LAH79 (Ohara)");
                result.RecommendedGlassPairs.Add("N-FK51A + N-KZFS11 (Schott)");
                result.RecommendedGlassPairs.Add("FCD1 + TAC8 (Hoya)");
            }
            else if (needsApochromat)
            {
                result.RecommendedGlassPairs.Add("S-FPL53 + S-LAL18 (Ohara)");
                result.RecommendedGlassPairs.Add("N-FK58 + N-LAF21 (Schott)");
            }
            else
            {
                result.RecommendedGlassPairs.Add("N-BK7 + N-SF5 (standard crown-flint)");
                result.RecommendedGlassPairs.Add("N-BAK4 + N-SF10 (improved)");
            }
            
            // Difficulty
            if (needsAnomalous)
            {
                result.Difficulty = "Very High";
                result.Recommendation = 
                    "Broad spectrum requires anomalous dispersion glasses. " +
                    "These are expensive and temperature sensitive. " +
                    "Consider narrowing wavelength range if possible.";
            }
            else if (needsApochromat)
            {
                result.Difficulty = "High";
                result.Recommendation = 
                    "Apochromatic correction required. Use ED glass pairs. " +
                    $"Secondary spectrum: {result.SecondarySpectrum:F3} waves.";
            }
            else if (result.MinDoubletsRequired > 0)
            {
                result.Difficulty = "Moderate";
                result.Recommendation = 
                    "Standard achromatic correction sufficient. " +
                    "Crown-flint doublet will provide good correction.";
            }
            else
            {
                result.Difficulty = "Low";
                result.Recommendation = 
                    "Narrow bandwidth - minimal chromatic correction needed.";
            }
            
            return result;
        }
        
        /// <summary>
        /// Analyze tolerance sensitivity
        /// </summary>
        private ToleranceResult AnalyzeTolerance(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new ToleranceResult();
            result.ToleranceGrade = input.ToleranceGrade;
            
            // Get tolerance values for grade
            result.ToleranceSpec = GetToleranceSpec(input.ToleranceGrade);
            
            // Compute sensitivities
            // These are empirical scaling based on NA and field
            double naFactor = Math.Pow(input.NAImage / 0.1, 2);
            double fieldFactor = Math.Pow(derived.FieldAngleImage / 0.1, 1);
            double magFactor = Math.Pow(Math.Abs(input.Magnification), 0.5);
            
            // Base sensitivities (waves per unit error)
            result.RadiusSensitivity = 0.5 * naFactor;           // per 0.1%
            result.ThicknessSensitivity = 0.3 * naFactor;        // per 0.01mm
            result.DecenterSensitivity = 1.0 * naFactor * fieldFactor;  // per 0.01mm
            result.TiltSensitivity = 0.8 * naFactor * fieldFactor;      // per arcmin
            result.IndexSensitivity = 0.2 * naFactor;            // per 0.0001
            result.AbbeSensitivity = 0.1 * naFactor;             // per 1%
            result.WedgeSensitivity = 0.4 * naFactor;            // per arcmin
            result.SurfaceIrregSensitivity = 1.0;                // per 0.1 wave
            
            // Compute total tolerance contribution (RSS)
            var spec = result.ToleranceSpec;
            double sumSq = 0;
            
            double radiusContrib = result.RadiusSensitivity * spec.RadiusPercent / 0.1;
            double thickContrib = result.ThicknessSensitivity * spec.ThicknessMM / 0.01;
            double decenterContrib = result.DecenterSensitivity * spec.DecenterMM / 0.01;
            double tiltContrib = result.TiltSensitivity * spec.TiltArcmin;
            double indexContrib = result.IndexSensitivity * spec.IndexError / 0.0001;
            double abbeContrib = result.AbbeSensitivity * spec.AbbePercent;
            double wedgeContrib = result.WedgeSensitivity * spec.WedgeArcmin;
            double irregContrib = result.SurfaceIrregSensitivity * spec.SurfaceIrregWaves / 0.1;
            
            sumSq = radiusContrib * radiusContrib +
                    thickContrib * thickContrib +
                    decenterContrib * decenterContrib +
                    tiltContrib * tiltContrib +
                    indexContrib * indexContrib +
                    abbeContrib * abbeContrib +
                    wedgeContrib * wedgeContrib +
                    irregContrib * irregContrib;
            
            result.ToleranceWavefrontRMS = Math.Sqrt(sumSq);
            
            // Design margin required
            double targetNominal = input.TargetRMSWavefront;
            result.DesignMarginRequired = Math.Sqrt(
                Math.Max(0, targetNominal * targetNominal - 
                        result.ToleranceWavefrontRMS * result.ToleranceWavefrontRMS));
            
            // Find dominant sensitivity
            var contribs = new Dictionary<string, double>
            {
                {"Radius", radiusContrib},
                {"Thickness", thickContrib},
                {"Decenter", decenterContrib},
                {"Tilt", tiltContrib},
                {"Index", indexContrib},
                {"Wedge", wedgeContrib},
                {"Surface irregularity", irregContrib}
            };
            
            double maxContrib = 0;
            foreach (var kv in contribs)
            {
                if (kv.Value > maxContrib)
                {
                    maxContrib = kv.Value;
                    result.DominantSensitivity = kv.Key;
                }
                if (kv.Value > 0.3 * result.ToleranceWavefrontRMS)
                {
                    result.CriticalTolerances.Add($"{kv.Key}: {kv.Value:F3} waves");
                }
            }
            
            // Estimated yield (simplified model)
            double toleranceRatio = result.ToleranceWavefrontRMS / targetNominal;
            if (toleranceRatio < 0.5)
                result.EstimatedYield = 0.95;
            else if (toleranceRatio < 0.7)
                result.EstimatedYield = 0.85;
            else if (toleranceRatio < 0.9)
                result.EstimatedYield = 0.70;
            else if (toleranceRatio < 1.0)
                result.EstimatedYield = 0.50;
            else
                result.EstimatedYield = 0.30;
            
            // Recommendation
            if (result.ToleranceWavefrontRMS > targetNominal)
            {
                result.Recommendation = 
                    $"Warning: Tolerance contribution ({result.ToleranceWavefrontRMS:F3}λ) " +
                    $"exceeds target ({targetNominal:F3}λ). " +
                    "Consider tighter tolerances or relaxed performance.";
            }
            else if (result.DesignMarginRequired < 0.02)
            {
                result.Recommendation = 
                    "Very tight tolerance budget. Nominal design must be near-perfect. " +
                    "Consider relaxing target or tightening key tolerances.";
            }
            else
            {
                result.Recommendation = 
                    $"Tolerance budget OK. Design to {result.DesignMarginRequired:F3}λ RMS " +
                    $"nominal to achieve {targetNominal:F3}λ as-built. " +
                    $"Estimated yield: {result.EstimatedYield:P0}";
            }
            
            return result;
        }
        
        /// <summary>
        /// Get tolerance specification for grade
        /// </summary>
        private ToleranceValues GetToleranceSpec(int grade)
        {
            switch (grade)
            {
                case 1: // Commercial
                    return new ToleranceValues
                    {
                        RadiusPercent = 0.5,
                        ThicknessMM = 0.10,
                        DecenterMM = 0.05,
                        TiltArcmin = 3.0,
                        IndexError = 0.001,
                        AbbePercent = 1.0,
                        WedgeArcmin = 3.0,
                        SurfaceIrregWaves = 0.5,
                        ScratchDig = 80
                    };
                    
                case 2: // Precision
                    return new ToleranceValues
                    {
                        RadiusPercent = 0.2,
                        ThicknessMM = 0.05,
                        DecenterMM = 0.02,
                        TiltArcmin = 1.0,
                        IndexError = 0.0005,
                        AbbePercent = 0.5,
                        WedgeArcmin = 1.0,
                        SurfaceIrregWaves = 0.25,
                        ScratchDig = 60
                    };
                    
                case 3: // High Precision
                    return new ToleranceValues
                    {
                        RadiusPercent = 0.1,
                        ThicknessMM = 0.02,
                        DecenterMM = 0.01,
                        TiltArcmin = 0.5,
                        IndexError = 0.0002,
                        AbbePercent = 0.2,
                        WedgeArcmin = 0.5,
                        SurfaceIrregWaves = 0.1,
                        ScratchDig = 40
                    };
                    
                default:
                    return GetToleranceSpec(2);
            }
        }
        
        /// <summary>
        /// Analyze relative illumination across the field
        /// </summary>
        private RelativeIlluminationResult AnalyzeRelativeIllumination(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new RelativeIlluminationResult();

            double[] fieldFracs = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
            double minRI = 1.0;
            double maxRI = 1.0;
            double beta0 = Math.Abs(input.Magnification);

            foreach (double frac in fieldFracs)
            {
                double theta = frac * derived.FieldAngleObject;
                var point = new IlluminationFieldPoint { NormalizedField = frac };

                // Natural vignetting: cos^4(θ) law
                // This accounts for: obliquity of rays, projection onto tilted surface,
                // inverse square law, and solid angle compression
                double cos4 = Math.Pow(Math.Cos(theta), 4);
                point.NaturalVignetting = cos4;

                // Distortion effect on illumination
                // Local magnification gradient affects light concentration
                // RI_distortion ∝ (β(θ)/β₀)² × |d(tanθ')/d(tanθ)|
                double betaTheta = ComputeDistortionMagnification(
                    theta, input.DistortionType, input.Magnification,
                    input.DistortionMax / 100.0, input.WDNominal);
                double betaRatio = Math.Abs(betaTheta / input.Magnification);

                // Compute magnification gradient effect
                double distortionEffect = 1.0;
                if (theta > 1e-6)
                {
                    // For different distortion types, the light concentration varies
                    switch (input.DistortionType)
                    {
                        case DistortionType.FTheta:
                            // F-theta compresses outer field → improves RI
                            distortionEffect = theta / Math.Tan(theta);
                            break;
                        case DistortionType.Orthographic:
                            // Orthographic: y' = f·sin(θ) → moderate effect
                            distortionEffect = Math.Cos(theta);
                            break;
                        case DistortionType.Barrel:
                            // Barrel stretches outer field → worsens RI
                            distortionEffect = 1.0 - 0.5 * input.DistortionMax / 100.0 * Math.Pow(Math.Tan(theta), 2);
                            break;
                        case DistortionType.Pincushion:
                            // Pincushion compresses outer field → improves RI
                            distortionEffect = 1.0 + 0.5 * input.DistortionMax / 100.0 * Math.Pow(Math.Tan(theta), 2);
                            break;
                        default:
                            distortionEffect = 1.0;
                            break;
                    }
                }
                point.DistortionEffect = distortionEffect;

                // Combined relative illumination
                // RI = cos^4(θ) × distortion_effect × telecentricity_factor
                double telecentricFactor = 1.0;
                if (input.TelecentricImage)
                {
                    // Image telecentricity maintains constant NA → improves uniformity
                    telecentricFactor = 1.0 / Math.Max(0.5, cos4);
                    telecentricFactor = Math.Min(telecentricFactor, 1.2); // Cap improvement
                }

                point.RelativeIllumination = cos4 * distortionEffect * telecentricFactor;
                point.RelativeIllumination = Math.Min(1.0, Math.Max(0.0, point.RelativeIllumination));

                result.IlluminationVsField.Add(point);

                minRI = Math.Min(minRI, point.RelativeIllumination);
                maxRI = Math.Max(maxRI, point.RelativeIllumination);
            }

            result.MinRelativeIllumination = minRI;
            result.IlluminationUniformity = 1.0 - (maxRI - minRI) / (maxRI + minRI);

            // Analyze contributions at edge of field
            var edgePoint = result.IlluminationVsField[result.IlluminationVsField.Count - 1];
            result.NaturalVignettingAtEdge = edgePoint.NaturalVignetting;
            result.DistortionContribution = Math.Abs(1.0 - edgePoint.DistortionEffect);

            // Pupil aberration contribution (simplified estimate)
            // Pupil coma and spherical aberration of pupil affect RI
            result.PupilAberrationContribution = 0.02 * Math.Pow(derived.NAObject, 2) *
                Math.Pow(derived.FieldAngleObject, 2);

            // Severity classification
            if (minRI > 0.9)
            {
                result.Severity = "Negligible";
                result.Recommendation = "Relative illumination is excellent (>90% at edge). " +
                    "No special measures needed.";
            }
            else if (minRI > 0.7)
            {
                result.Severity = "Minor";
                result.Recommendation = "Moderate RI roll-off. Consider image-space telecentricity " +
                    "or vignetting control if uniformity is critical.";
            }
            else if (minRI > 0.5)
            {
                result.Severity = "Moderate";
                result.Recommendation = $"Significant RI drop to {minRI:P0} at edge. " +
                    "May need pupil aberration control or post-processing correction.";
            }
            else
            {
                result.Severity = "Severe";
                result.Recommendation = $"Severe RI roll-off ({minRI:P0} at edge). " +
                    "Consider telecentricity, reducing field angle, or accepting non-uniform illumination.";
            }

            return result;
        }

        /// <summary>
        /// Analyze differential distortion across the field
        /// </summary>
        private DifferentialDistortionResult AnalyzeDifferentialDistortion(
            FeasibilityInput input, DerivedQuantities derived)
        {
            var result = new DifferentialDistortionResult();

            double[] fieldFracs = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
            double beta0 = input.Magnification;
            double maxDiff = 0;
            double maxLocalVariation = 0;

            // Compute ideal magnification gradient (for rectilinear projection)
            // For rectilinear: y' = f·tan(θ), so β_ideal = f/cos²(θ) for image height
            // The "ideal" differential is dβ/dθ = 0 for constant magnification

            List<double> magnifications = new List<double>();
            List<double> thetas = new List<double>();

            foreach (double frac in fieldFracs)
            {
                double theta = frac * derived.FieldAngleObject;
                thetas.Add(theta);

                double betaTheta = ComputeDistortionMagnification(
                    theta, input.DistortionType, beta0,
                    input.DistortionMax / 100.0, input.WDNominal);
                magnifications.Add(betaTheta);
            }

            // Compute differential distortion (magnification gradient)
            for (int i = 0; i < fieldFracs.Length; i++)
            {
                var point = new DifferentialDistortionPoint
                {
                    NormalizedField = fieldFracs[i],
                    LocalMagnification = magnifications[i]
                };

                // Compute gradient using central differences
                double dBeta_dTheta;
                if (i == 0)
                {
                    // Forward difference at origin
                    if (thetas[1] - thetas[0] > 1e-10)
                        dBeta_dTheta = (magnifications[1] - magnifications[0]) / (thetas[1] - thetas[0]);
                    else
                        dBeta_dTheta = 0;
                }
                else if (i == fieldFracs.Length - 1)
                {
                    // Backward difference at edge
                    dBeta_dTheta = (magnifications[i] - magnifications[i - 1]) / (thetas[i] - thetas[i - 1]);
                }
                else
                {
                    // Central difference
                    dBeta_dTheta = (magnifications[i + 1] - magnifications[i - 1]) / (thetas[i + 1] - thetas[i - 1]);
                }

                point.MagnificationGradient = dBeta_dTheta;

                // Differential distortion as deviation from constant magnification
                // Expressed as percentage variation
                if (Math.Abs(beta0) > 1e-10)
                {
                    point.DifferentialDistortion = Math.Abs(dBeta_dTheta) / Math.Abs(beta0) * 100.0;
                }

                // Tangential vs radial magnification difference
                // In presence of astigmatism, these differ
                // Simplified model: T-R difference ∝ field² × NA²
                double theta = thetas[i];
                double astigFactor = Math.Pow(theta, 2) * Math.Pow(derived.NAObject, 2);
                point.TangentialMagnification = magnifications[i] * (1.0 + 0.1 * astigFactor);
                point.RadialMagnification = magnifications[i] * (1.0 - 0.1 * astigFactor);

                result.DifferentialVsField.Add(point);

                maxDiff = Math.Max(maxDiff, point.DifferentialDistortion);

                // Local variation from nominal
                double localVar = Math.Abs((magnifications[i] - beta0) / beta0) * 100.0;
                maxLocalVariation = Math.Max(maxLocalVariation, localVar);
            }

            result.MaxDifferentialDistortion = maxDiff;
            result.LocalMagnificationVariation = maxLocalVariation;

            // T-R difference at edge
            var edgePt = result.DifferentialVsField[result.DifferentialVsField.Count - 1];
            if (Math.Abs(edgePt.LocalMagnification) > 1e-10)
            {
                result.TangentialRadialDifference = Math.Abs(
                    (edgePt.TangentialMagnification - edgePt.RadialMagnification) /
                    edgePt.LocalMagnification) * 100.0;
            }

            result.IsUniform = maxDiff < 0.5;

            // Metrology impact assessment
            if (maxDiff < 0.1)
            {
                result.MetrologyImpact = "Excellent - suitable for precision metrology";
                result.Recommendation = "Differential distortion is negligible. " +
                    "Suitable for high-precision measurement applications.";
            }
            else if (maxDiff < 0.5)
            {
                result.MetrologyImpact = "Good - suitable for general metrology with calibration";
                result.Recommendation = "Minor differential distortion. Calibration can correct " +
                    "for measurement applications. No fundamental conflict.";
            }
            else if (maxDiff < 2.0)
            {
                result.MetrologyImpact = "Moderate - requires careful calibration";
                result.Recommendation = $"Noticeable differential distortion ({maxDiff:F1}%). " +
                    "Measurement accuracy depends on calibration quality.";
            }
            else
            {
                result.MetrologyImpact = "Limited - challenging for precision metrology";
                result.Recommendation = $"Significant differential distortion ({maxDiff:F1}%). " +
                    "Consider rectilinear projection if measurement accuracy is critical.";
            }

            return result;
        }

        /// <summary>
        /// Check for conflicts between relative illumination, distortion, and differential distortion
        /// </summary>
        private IlluminationDistortionConflictResult CheckIlluminationDistortionConflict(
            FeasibilityInput input, DerivedQuantities derived,
            RelativeIlluminationResult riResult, DifferentialDistortionResult ddResult,
            SineConditionResult sineResult, TelecentricityResult teleResult)
        {
            var result = new IlluminationDistortionConflictResult();
            double conflictScore = 0;

            // ═══════════════════════════════════════════════════════════════════
            // CONFLICT 1: Illumination vs Distortion Type
            // ═══════════════════════════════════════════════════════════════════
            // Some distortion types inherently worsen or improve relative illumination

            result.IlluminationDistortionConflict = false;
            bool wantsUniformIllumination = riResult.MinRelativeIllumination < 0.7 ||
                                            input.TelecentricImage;

            switch (input.DistortionType)
            {
                case DistortionType.Barrel:
                    // Barrel distortion stretches outer field → worsens RI
                    if (wantsUniformIllumination && input.DistortionMax > 2.0)
                    {
                        result.IlluminationDistortionConflict = true;
                        result.IlluminationDistortionDetail =
                            $"Barrel distortion ({input.DistortionMax:F1}%) stretches outer field, " +
                            "reducing relative illumination at edges. This conflicts with uniform " +
                            "illumination requirements.";
                        result.Conflicts.Add("Barrel distortion worsens edge illumination");
                        result.Recommendations.Add("Consider f-theta or pincushion if uniform RI is needed");
                        conflictScore += 15;
                    }
                    break;

                case DistortionType.FTheta:
                    // F-theta compresses outer field → helps RI but conflicts with sine condition
                    if (sineResult.MaxBetaDiscrepancy > 0.02)
                    {
                        // Note: F-theta actually helps illumination but hurts coma
                        result.IlluminationDistortionDetail =
                            "F-theta distortion improves edge illumination but creates " +
                            $"sine condition conflict (β discrepancy: {sineResult.MaxBetaDiscrepancy * 100:F1}%).";
                    }
                    break;

                case DistortionType.Orthographic:
                    // Orthographic has moderate illumination effect
                    if (wantsUniformIllumination && derived.FieldAngleObject > 0.3)
                    {
                        result.IlluminationDistortionConflict = true;
                        result.IlluminationDistortionDetail =
                            "Orthographic projection at wide field angles causes significant " +
                            "illumination roll-off (follows cos law).";
                        result.Conflicts.Add("Orthographic projection causes RI roll-off at wide angles");
                        conflictScore += 10;
                    }
                    break;
            }

            // ═══════════════════════════════════════════════════════════════════
            // CONFLICT 2: Differential Distortion vs Sine Condition
            // ═══════════════════════════════════════════════════════════════════
            // The sine condition requires specific magnification variation for coma correction
            // Distortion types impose different magnification variations
            // The mismatch creates both coma AND differential distortion issues

            result.DifferentialSineConditionConflict = false;

            // Compare required magnification gradients
            // Sine condition: β(θ) = β₀/cos(θ) → dβ/dθ = β₀·sin(θ)/cos²(θ)
            // Distortion: varies by type

            double theta_edge = derived.FieldAngleObject;
            double sineBetaGradient = 0;
            if (theta_edge > 1e-6)
            {
                sineBetaGradient = Math.Abs(input.Magnification) *
                    Math.Sin(theta_edge) / Math.Pow(Math.Cos(theta_edge), 2);
            }

            // Get actual distortion-imposed gradient at edge
            double actualGradient = 0;
            if (ddResult.DifferentialVsField.Count > 0)
            {
                actualGradient = Math.Abs(ddResult.DifferentialVsField[
                    ddResult.DifferentialVsField.Count - 1].MagnificationGradient);
            }

            double gradientMismatch = 0;
            if (sineBetaGradient > 1e-6)
            {
                gradientMismatch = Math.Abs(actualGradient - sineBetaGradient) / sineBetaGradient;
            }

            if (gradientMismatch > 0.1 && sineResult.InducedComaWaves > 0.05)
            {
                result.DifferentialSineConditionConflict = true;
                result.DifferentialSineConditionDetail =
                    $"The distortion type requires a magnification gradient that differs " +
                    $"from the sine condition by {gradientMismatch * 100:F0}%. This causes both " +
                    $"coma ({sineResult.InducedComaWaves:F2} waves) and non-uniform differential " +
                    $"distortion ({ddResult.MaxDifferentialDistortion:F1}%).";
                result.Conflicts.Add("Distortion gradient conflicts with sine condition");
                result.Recommendations.Add("Accept coma at edge OR relax distortion requirement");
                conflictScore += 20 * gradientMismatch;
            }

            // ═══════════════════════════════════════════════════════════════════
            // CONFLICT 3: Illumination vs Telecentricity
            // ═══════════════════════════════════════════════════════════════════
            // Double telecentricity can't exist with m≠1, but single telecentricity
            // may conflict with illumination requirements

            result.IlluminationTelecentricityConflict = false;

            if (input.TelecentricObject && !input.TelecentricImage)
            {
                // Object telecentric without image telecentric
                // May have illumination issues due to pupil position
                if (riResult.MinRelativeIllumination < 0.6)
                {
                    result.IlluminationTelecentricityConflict = true;
                    result.IlluminationTelecentricityDetail =
                        "Object-side telecentricity without image-side telecentricity " +
                        "can lead to illumination non-uniformity due to stop position.";
                    result.Conflicts.Add("Object telecentricity without image telecentricity affects RI");
                    result.Recommendations.Add("Add image-side telecentricity or accept RI variation");
                    conflictScore += 10;
                }
            }

            if (input.TelecentricObject && input.TelecentricImage &&
                Math.Abs(Math.Abs(input.Magnification) - 1.0) > 0.01)
            {
                // Already flagged in telecentricity analysis, but affects illumination too
                result.IlluminationTelecentricityConflict = true;
                result.IlluminationTelecentricityDetail =
                    "True double telecentricity requires |m|=1. With m≠1, " +
                    "illumination uniformity cannot be optimized simultaneously.";
                conflictScore += 15;
            }

            // ═══════════════════════════════════════════════════════════════════
            // ADDITIONAL COUPLING: Differential distortion affects RI measurement
            // ═══════════════════════════════════════════════════════════════════

            if (ddResult.MaxDifferentialDistortion > 1.0 && riResult.Severity != "Negligible")
            {
                result.Conflicts.Add(
                    $"Differential distortion ({ddResult.MaxDifferentialDistortion:F1}%) " +
                    "complicates relative illumination characterization");
                result.Recommendations.Add(
                    "Characterize RI and distortion together for accurate calibration");
                conflictScore += 5;
            }

            // Overall assessment
            result.ConflictScore = Math.Min(100, conflictScore);
            result.HasConflict = result.ConflictScore > 10;

            if (result.ConflictScore < 10)
            {
                result.ConflictSeverity = "None";
            }
            else if (result.ConflictScore < 25)
            {
                result.ConflictSeverity = "Minor";
            }
            else if (result.ConflictScore < 50)
            {
                result.ConflictSeverity = "Moderate";
            }
            else
            {
                result.ConflictSeverity = "Severe";
            }

            // Add general recommendations if conflicts exist
            if (result.HasConflict && result.Recommendations.Count == 0)
            {
                result.Recommendations.Add(
                    "Review distortion, illumination, and coma requirements together - " +
                    "they are fundamentally coupled through the sine condition.");
            }

            return result;
        }

        /// <summary>
        /// Compute aberration floors
        /// </summary>
        private AberrationFloors ComputeAberrationFloors(
            FeasibilityInput input, DerivedQuantities derived,
            SineConditionResult sine, ReversibilityResult rev,
            ChromaticResult chrom, PetzvalResult petzval)
        {
            var floors = new AberrationFloors();
            
            // Spherical: generally correctable
            floors.Spherical = 0.0;
            
            // Coma: limited by sine condition
            floors.Coma = sine.InducedComaWaves;
            
            // Astigmatism: limited by Petzval tradeoff
            floors.Astigmatism = input.FlatField ? 0.02 : 0.0;
            
            // Field curvature
            floors.FieldCurvature = input.FlatField ? 
                Math.Min(0.05, petzval.PetzvalDefocusWaves * 0.05) : 0.0;
            
            // Distortion
            floors.Distortion = 0.0;
            
            // Chromatic
            floors.ChromaticAxial = chrom.AxialColorCorrected;
            floors.ChromaticLateral = chrom.LateralColorCorrected;
            
            // Reversibility contribution
            double revContrib = rev.Constrained ? rev.UnavoidableRMSWaves : 0.0;
            
            // Total RMS (RSS)
            floors.TotalRMS = Math.Sqrt(
                floors.Spherical * floors.Spherical +
                floors.Coma * floors.Coma +
                floors.Astigmatism * floors.Astigmatism +
                floors.FieldCurvature * floors.FieldCurvature +
                floors.ChromaticAxial * floors.ChromaticAxial +
                floors.ChromaticLateral * floors.ChromaticLateral +
                revContrib * revContrib);
            
            floors.DiffractionLimited = floors.TotalRMS < 0.07;
            
            // Identify limiting aberrations
            var limiting = new List<string>();
            if (floors.Coma > 0.02) limiting.Add("coma");
            if (floors.ChromaticAxial > 0.02) limiting.Add("axial color");
            if (floors.ChromaticLateral > 0.02) limiting.Add("lateral color");
            if (revContrib > 0.02) limiting.Add("reversibility");
            if (floors.FieldCurvature > 0.02) limiting.Add("field curvature");
            
            floors.LimitingAberrations = limiting.Count > 0 ?
                string.Join(", ", limiting) : "none (diffraction limited achievable)";
            
            return floors;
        }
        
        /// <summary>
        /// Estimate design complexity
        /// </summary>
        private ComplexityEstimate EstimateComplexity(
            FeasibilityInput input, FeasibilityOutput output)
        {
            var est = new ComplexityEstimate();
            
            // Base elements from NA and field
            int baseElements = 3;
            
            if (input.NAImage > 0.5) 
            {
                baseElements += 4;
                est.ElementBreakdown["High NA"] = 4;
            }
            else if (input.NAImage > 0.3)
            {
                baseElements += 2;
                est.ElementBreakdown["Moderate NA"] = 2;
            }
            else if (input.NAImage > 0.15)
            {
                baseElements += 1;
                est.ElementBreakdown["Low-moderate NA"] = 1;
            }
            
            double fieldAngleDeg = output.Derived.FieldAngleImage * RAD_TO_DEG;
            if (fieldAngleDeg > 30)
            {
                baseElements += 3;
                est.ElementBreakdown["Wide field"] = 3;
            }
            else if (fieldAngleDeg > 15)
            {
                baseElements += 2;
                est.ElementBreakdown["Moderate field"] = 2;
            }
            else if (fieldAngleDeg > 7)
            {
                baseElements += 1;
                est.ElementBreakdown["Field correction"] = 1;
            }
            
            // Constraint contributions
            if (input.FlatField && output.Petzval.PetzvalDefocusWaves > 1)
            {
                baseElements += 2;
                est.ElementBreakdown["Field flattener"] = 2;
            }
            
            if (input.TelecentricObject)
            {
                baseElements += 1;
                est.ElementBreakdown["Object telecentricity"] = 1;
            }
            
            if (input.TelecentricImage)
            {
                baseElements += 1;
                est.ElementBreakdown["Image telecentricity"] = 1;
            }
            
            if (input.DistortionType != DistortionType.None)
            {
                baseElements += 1;
                est.ElementBreakdown["Distortion control"] = 1;
            }
            
            if (input.IsReversible)
            {
                baseElements = 2 * ((baseElements + 1) / 2);
                est.DominantDrivers.Add("Symmetry (reversible)");
            }
            
            if (output.Derived.WDRatio > 2)
            {
                baseElements += 1;
                est.ElementBreakdown["WD range flexibility"] = 1;
            }
            
            // Chromatic correction
            if (output.Chromatic.MinDoubletsRequired > 0)
            {
                baseElements += output.Chromatic.MinDoubletsRequired;
                est.ElementBreakdown["Chromatic correction"] = output.Chromatic.MinDoubletsRequired;
            }
            
            est.MinElements = Math.Max(3, baseElements - 2);
            est.MaxElements = baseElements + 3;
            est.EstimatedGroups = (baseElements + 2) / 3;
            est.Confidence = "±30% (pre-design estimate)";
            
            // Dominant drivers
            if (input.NAImage > 0.3)
                est.DominantDrivers.Add("High NA");
            if (fieldAngleDeg > 20)
                est.DominantDrivers.Add("Wide field");
            if (input.FlatField)
                est.DominantDrivers.Add("Flat field");
            if (output.Chromatic.NeedsApochromat)
                est.DominantDrivers.Add("Apochromatic correction");
            
            // Similar designs
            if (Math.Abs(input.Magnification) > 2 && input.TelecentricImage)
            {
                est.SimilarDesigns.Add("Microscope objective");
            }
            if (input.DistortionType == DistortionType.FTheta)
            {
                est.SimilarDesigns.Add("F-theta scan lens");
            }
            if (input.IsReversible && Math.Abs(Math.Abs(input.Magnification) - 1) < 0.1)
            {
                est.SimilarDesigns.Add("Symmetric relay");
            }
            if (input.NAImage < 0.1 && fieldAngleDeg > 20)
            {
                est.SimilarDesigns.Add("Wide-angle projection lens");
            }
            if (est.SimilarDesigns.Count == 0)
            {
                est.SimilarDesigns.Add("Custom finite conjugate imaging lens");
            }
            
            // Cost class
            if (baseElements > 10 || output.Chromatic.NeedsAnomalousDispersion)
                est.CostClass = 5;
            else if (baseElements > 7 || output.Chromatic.NeedsApochromat)
                est.CostClass = 4;
            else if (baseElements > 5)
                est.CostClass = 3;
            else if (baseElements > 3)
                est.CostClass = 2;
            else
                est.CostClass = 1;
            
            return est;
        }
        
        /// <summary>
        /// Compute overall feasibility score
        /// </summary>
        private void ComputeFeasibilityScore(FeasibilityOutput output)
        {
            double score = 100;
            
            // Sine condition penalty
            if (!output.SineCondition.Compatible)
            {
                score -= 20 * output.SineCondition.MaxBetaDiscrepancy;
                output.LimitingFactors.Add("Distortion/coma tradeoff");
            }
            
            // Reversibility penalty
            if (output.Reversibility.Constrained && 
                output.Reversibility.UnavoidableRMSWaves > 0.07)
            {
                score -= 100 * output.Reversibility.UnavoidableRMSWaves;
                output.LimitingFactors.Add("Reversibility-induced aberration");
            }
            
            // Working distance penalty
            if (output.WorkingDistance.NeedsRefocusCompensation)
            {
                score -= 10;
                output.LimitingFactors.Add("WD range requires focus compensation");
            }
            if (output.WorkingDistance.WDRatio > 2)
            {
                score -= 5 * (output.WorkingDistance.WDRatio - 2);
            }
            
            // Petzval penalty
            if (output.Input.FlatField && output.Petzval.PetzvalDefocusWaves > 2)
            {
                score -= 5 * output.Petzval.PetzvalDefocusWaves;
                output.LimitingFactors.Add("Demanding flat field correction");
            }
            
            // Telecentricity penalty
            if (output.Input.TelecentricObject && output.Input.TelecentricImage &&
                !output.Telecentricity.DoubleTelecentricFeasible)
            {
                score -= 30;
                output.LimitingFactors.Add("Double telecentricity not achievable");
            }
            
            // Chromatic penalty
            if (output.Chromatic.NeedsAnomalousDispersion)
            {
                score -= 15;
                output.LimitingFactors.Add("Requires anomalous dispersion glasses");
            }
            else if (output.Chromatic.NeedsApochromat)
            {
                score -= 5;
            }
            
            // Tolerance penalty
            if (output.Tolerance.ToleranceWavefrontRMS > output.Input.TargetRMSWavefront)
            {
                score -= 20;
                output.LimitingFactors.Add("Tolerance budget exceeds target");
            }
            else if (output.Tolerance.DesignMarginRequired < 0.02)
            {
                score -= 10;
                output.LimitingFactors.Add("Very tight tolerance budget");
            }
            
            // Aberration floor penalty
            if (!output.AberrationFloors.DiffractionLimited)
            {
                score -= 20;
                output.LimitingFactors.Add($"Aberration floor > diffraction limit");
            }

            // Relative illumination penalty
            if (output.RelativeIllumination.MinRelativeIllumination < 0.5)
            {
                score -= 15;
                output.LimitingFactors.Add($"Severe relative illumination drop ({output.RelativeIllumination.MinRelativeIllumination:P0} at edge)");
            }
            else if (output.RelativeIllumination.MinRelativeIllumination < 0.7)
            {
                score -= 5;
                output.LimitingFactors.Add("Moderate relative illumination non-uniformity");
            }

            // Differential distortion penalty
            if (output.DifferentialDistortion.MaxDifferentialDistortion > 2.0)
            {
                score -= 10;
                output.LimitingFactors.Add($"High differential distortion ({output.DifferentialDistortion.MaxDifferentialDistortion:F1}%)");
            }

            // Illumination-distortion conflict penalty
            if (output.IlluminationDistortionConflict.HasConflict)
            {
                score -= output.IlluminationDistortionConflict.ConflictScore * 0.3;
                if (output.IlluminationDistortionConflict.ConflictSeverity == "Severe")
                {
                    output.LimitingFactors.Add("Severe illumination/distortion/coma conflicts");
                }
                else if (output.IlluminationDistortionConflict.ConflictSeverity == "Moderate")
                {
                    output.LimitingFactors.Add("Moderate illumination/distortion constraints");
                }
            }

            output.FeasibilityScore = Math.Max(0, Math.Min(100, score));
            output.IsFeasible = output.FeasibilityScore >= 50;
        }
        
        /// <summary>
        /// Generate recommendations
        /// </summary>
        private void GenerateRecommendations(FeasibilityOutput output)
        {
            // Design recommendations
            if (output.IsFeasible)
            {
                output.Recommendations.Add(
                    $"Design appears feasible (score: {output.FeasibilityScore:F0}/100)");
            }
            else
            {
                output.Recommendations.Add(
                    $"Design may be challenging (score: {output.FeasibilityScore:F0}/100). " +
                    "Consider relaxing constraints.");
            }
            
            // Starting form suggestions
            if (output.Complexity.SimilarDesigns.Count > 0)
            {
                output.Recommendations.Add(
                    $"Starting forms to consider: {string.Join(", ", output.Complexity.SimilarDesigns)}");
            }
            
            // Specific recommendations from analysis modules
            if (output.AberrationFloors.TotalRMS < output.Input.TargetRMSWavefront * 0.7)
            {
                output.Recommendations.Add(
                    $"Target {output.Tolerance.DesignMarginRequired:F3}λ RMS nominal " +
                    "to meet as-built target with tolerance margin.");
            }
            
            if (output.Chromatic.MinDoubletsRequired > 0)
            {
                output.Recommendations.Add(
                    $"Include {output.Chromatic.MinDoubletsRequired} achromatic doublet(s) " +
                    "for chromatic correction.");
            }
            
            if (output.Input.FlatField && output.Petzval.PetzvalDefocusWaves > 1)
            {
                output.Recommendations.Add(
                    "Consider field flattener group near image plane.");
            }

            // Relative illumination recommendations
            if (output.RelativeIllumination.Severity == "Severe" ||
                output.RelativeIllumination.Severity == "Moderate")
            {
                output.Recommendations.Add(
                    $"Address relative illumination ({output.RelativeIllumination.MinRelativeIllumination:P0} at edge): " +
                    "consider image-side telecentricity or vignetting control.");
            }

            // Differential distortion recommendations
            if (!output.DifferentialDistortion.IsUniform)
            {
                output.Recommendations.Add(
                    $"Differential distortion ({output.DifferentialDistortion.MaxDifferentialDistortion:F1}%) " +
                    "may affect measurement applications. Consider calibration requirements.");
            }

            // Illumination-distortion conflict recommendations
            if (output.IlluminationDistortionConflict.HasConflict)
            {
                foreach (var rec in output.IlluminationDistortionConflict.Recommendations)
                {
                    output.Recommendations.Add(rec);
                }
            }

            // Relaxation suggestions if not feasible
            if (!output.IsFeasible)
            {
                if (output.LimitingFactors.Contains("Demanding flat field correction"))
                {
                    output.RelaxationSuggestions.Add(
                        "Allow curved field or reduced field size");
                }
                
                if (output.LimitingFactors.Contains("Reversibility-induced aberration"))
                {
                    output.RelaxationSuggestions.Add(
                        "Consider asymmetric design (remove reversibility)");
                }
                
                if (output.LimitingFactors.Contains("Distortion/coma tradeoff"))
                {
                    output.RelaxationSuggestions.Add(
                        "Relax distortion requirement or accept edge-of-field coma");
                }
                
                if (output.LimitingFactors.Contains("Double telecentricity not achievable"))
                {
                    output.RelaxationSuggestions.Add(
                        "Use single-side telecentricity or quasi-telecentric design");
                }
                
                if (output.LimitingFactors.Contains("WD range requires focus compensation"))
                {
                    output.RelaxationSuggestions.Add(
                        "Reduce WD range or add motorized focus");
                }

                if (output.RelativeIllumination.MinRelativeIllumination < 0.7)
                {
                    output.RelaxationSuggestions.Add(
                        "Accept non-uniform illumination or add image-side telecentricity");
                }

                if (output.DifferentialDistortion.MaxDifferentialDistortion > 1.0)
                {
                    output.RelaxationSuggestions.Add(
                        "Use rectilinear projection or accept calibration requirements");
                }

                if (output.IlluminationDistortionConflict.HasConflict)
                {
                    output.RelaxationSuggestions.Add(
                        "Prioritize either illumination uniformity OR distortion type - both may not be achievable");
                }

                output.RelaxationSuggestions.Add(
                    "Reduce NA (easiest path to simplification)");
                output.RelaxationSuggestions.Add(
                    "Reduce field size");
            }
        }
    }
}
