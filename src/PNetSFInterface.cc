#include "Tools/Tools/interface/PNetSFInterface.h"
#include <iostream> // FIXME: remove

// Constructor
PNetSFInterface::PNetSFInterface(const std::string period) :
    // First initialize all SF vectors to 1...
    scaleFactorsHH_lowPt_ (9, 1.),
    scaleFactorsHH_medPt_ (9, 1.),
    scaleFactorsHH_highPt_(9, 1.),
    scaleFactorsTT_lowPt_ (9, 1.),
    scaleFactorsTT_medPt_ (9, 1.),
    scaleFactorsTT_highPt_(9, 1.),
    scaleFactorsDY_lowPt_ (9, 1.),
    scaleFactorsDY_medPt_ (9, 1.),
    scaleFactorsDY_highPt_(9, 1.)
{
    // ...then fill them with actual values depending on Period
    FillHHlikeSFs(period);
    FillTTlikeSFs(period);
    FillDYlikeSFs(period);
}

// Destructor
PNetSFInterface::~PNetSFInterface() {}

// Fill SFs for HH-like samples
void PNetSFInterface::FillHHlikeSFs(const std::string period)
{
    // Set scale factor values depending on the Era
    if (period == "2016preVFP")
    {
        // pT < 500
        scaleFactorsHH_lowPt_[0] = (1.054);            // HP central
        scaleFactorsHH_lowPt_[1] = (1.054 + ( 0.080)); // HP up
        scaleFactorsHH_lowPt_[2] = (1.054 + (-0.077)); // HP down
        scaleFactorsHH_lowPt_[3] = (1.052);            // MP central
        scaleFactorsHH_lowPt_[4] = (1.052 + ( 0.087)); // MP up
        scaleFactorsHH_lowPt_[5] = (1.052 + (-0.081)); // MP down
        scaleFactorsHH_lowPt_[6] = (1.032);            // LP central
        scaleFactorsHH_lowPt_[7] = (1.032 + ( 0.096)); // LP up
        scaleFactorsHH_lowPt_[8] = (1.032 + (-0.090)); // LP down

        // 500 <= pT < 600
        scaleFactorsHH_medPt_[0] = (1.139);            // HP central
        scaleFactorsHH_medPt_[1] = (1.139 + ( 0.083)); // HP up
        scaleFactorsHH_medPt_[2] = (1.139 + (-0.081)); // HP down
        scaleFactorsHH_medPt_[3] = (1.068);            // MP central
        scaleFactorsHH_medPt_[4] = (1.068 + ( 0.078)); // MP up
        scaleFactorsHH_medPt_[5] = (1.068 + (-0.073)); // MP down
        scaleFactorsHH_medPt_[6] = (1.062);            // LP central
        scaleFactorsHH_medPt_[7] = (1.062 + ( 0.092)); // LP up
        scaleFactorsHH_medPt_[8] = (1.062 + (-0.082)); // LP down

        // pT >= 600
        scaleFactorsHH_highPt_[0] = (1.049);            // HP central
        scaleFactorsHH_highPt_[1] = (1.049 + ( 0.133)); // HP up
        scaleFactorsHH_highPt_[2] = (1.049 + (-0.130)); // HP down
        scaleFactorsHH_highPt_[3] = (0.996);            // MP central
        scaleFactorsHH_highPt_[4] = (0.996 + ( 0.101)); // MP up
        scaleFactorsHH_highPt_[5] = (0.996 + (-0.097)); // MP down
        scaleFactorsHH_highPt_[6] = (1.002);            // LP central
        scaleFactorsHH_highPt_[7] = (1.002 + ( 0.106)); // LP up
        scaleFactorsHH_highPt_[8] = (1.002 + (-0.101)); // LP down
    }
    else if (period == "2016postVFP")
    {
        // pT < 500
        scaleFactorsHH_lowPt_[0] = (1.031);            // HP central
        scaleFactorsHH_lowPt_[1] = (1.031 + ( 0.050)); // HP up
        scaleFactorsHH_lowPt_[2] = (1.031 + (-0.046)); // HP down
        scaleFactorsHH_lowPt_[3] = (1.029);            // MP central
        scaleFactorsHH_lowPt_[4] = (1.029 + ( 0.051)); // MP up
        scaleFactorsHH_lowPt_[5] = (1.029 + (-0.045)); // MP down
        scaleFactorsHH_lowPt_[6] = (1.031);            // LP central
        scaleFactorsHH_lowPt_[7] = (1.031 + ( 0.058)); // LP up
        scaleFactorsHH_lowPt_[8] = (1.031 + (-0.050)); // LP down

        // 500 <= pT < 600
        scaleFactorsHH_medPt_[0] = (1.055);            // HP central
        scaleFactorsHH_medPt_[1] = (1.055 + ( 0.069)); // HP up
        scaleFactorsHH_medPt_[2] = (1.055 + (-0.067)); // HP down
        scaleFactorsHH_medPt_[3] = (1.070);            // MP central
        scaleFactorsHH_medPt_[4] = (1.070 + ( 0.066)); // MP up
        scaleFactorsHH_medPt_[5] = (1.070 + (-0.062)); // MP down
        scaleFactorsHH_medPt_[6] = (1.089);            // LP central
        scaleFactorsHH_medPt_[7] = (1.089 + ( 0.076)); // LP up
        scaleFactorsHH_medPt_[8] = (1.089 + (-0.068)); // LP down

        // pT >= 600
        scaleFactorsHH_highPt_[0] = (1.088);            // HP central
        scaleFactorsHH_highPt_[1] = (1.088 + ( 0.076)); // HP up
        scaleFactorsHH_highPt_[2] = (1.088 + (-0.072)); // HP down
        scaleFactorsHH_highPt_[3] = (1.077);            // MP central
        scaleFactorsHH_highPt_[4] = (1.077 + ( 0.067)); // MP up
        scaleFactorsHH_highPt_[5] = (1.077 + (-0.059)); // MP down
        scaleFactorsHH_highPt_[6] = (1.057);            // LP central
        scaleFactorsHH_highPt_[7] = (1.057 + ( 0.077)); // LP up
        scaleFactorsHH_highPt_[8] = (1.057 + (-0.056)); // LP down
    }
    else if (period == "2017")
    {
        // pT < 500
        scaleFactorsHH_lowPt_[0] = (1.055);            // HP central
        scaleFactorsHH_lowPt_[1] = (1.055 + ( 0.057)); // HP up
        scaleFactorsHH_lowPt_[2] = (1.055 + (-0.054)); // HP down
        scaleFactorsHH_lowPt_[3] = (1.006);            // MP central
        scaleFactorsHH_lowPt_[4] = (1.006 + ( 0.052)); // MP up
        scaleFactorsHH_lowPt_[5] = (1.006 + (-0.052)); // MP down
        scaleFactorsHH_lowPt_[6] = (0.966);            // LP central
        scaleFactorsHH_lowPt_[7] = (0.966 + ( 0.055)); // LP up
        scaleFactorsHH_lowPt_[8] = (0.966 + (-0.057)); // LP down

        // 500 <= pT < 600
        scaleFactorsHH_medPt_[0] = (1.067);            // HP central
        scaleFactorsHH_medPt_[1] = (1.067 + ( 0.057)); // HP up
        scaleFactorsHH_medPt_[2] = (1.067 + (-0.055)); // HP down
        scaleFactorsHH_medPt_[3] = (1.051);            // MP central
        scaleFactorsHH_medPt_[4] = (1.051 + ( 0.056)); // MP up
        scaleFactorsHH_medPt_[5] = (1.051 + (-0.055)); // MP down
        scaleFactorsHH_medPt_[6] = (1.021);            // LP central
        scaleFactorsHH_medPt_[7] = (1.021 + ( 0.053)); // LP up
        scaleFactorsHH_medPt_[8] = (1.021 + (-0.052)); // LP down

        // pT >= 600
        scaleFactorsHH_highPt_[0] = (1.045);            // HP central
        scaleFactorsHH_highPt_[1] = (1.045 + ( 0.045)); // HP up
        scaleFactorsHH_highPt_[2] = (1.045 + (-0.046)); // HP down
        scaleFactorsHH_highPt_[3] = (0.991);            // MP central
        scaleFactorsHH_highPt_[4] = (0.991 + ( 0.038)); // MP up
        scaleFactorsHH_highPt_[5] = (0.991 + (-0.043)); // MP down
        scaleFactorsHH_highPt_[6] = (0.979);            // LP central
        scaleFactorsHH_highPt_[7] = (0.979 + ( 0.035)); // LP up
        scaleFactorsHH_highPt_[8] = (0.979 + (-0.038)); // LP down
    }
    else if (period == "2018" || period == "2022" || period == "2023")
    {
        // pT < 500
        scaleFactorsHH_lowPt_[0] = (0.994);            // HP central
        scaleFactorsHH_lowPt_[1] = (0.994 + ( 0.064)); // HP up
        scaleFactorsHH_lowPt_[2] = (0.994 + (-0.064)); // HP down
        scaleFactorsHH_lowPt_[3] = (0.966);            // MP central
        scaleFactorsHH_lowPt_[4] = (0.966 + ( 0.056)); // MP up
        scaleFactorsHH_lowPt_[5] = (0.966 + (-0.057)); // MP down
        scaleFactorsHH_lowPt_[6] = (0.921);            // LP central
        scaleFactorsHH_lowPt_[7] = (0.921 + ( 0.071)); // LP up
        scaleFactorsHH_lowPt_[8] = (0.921 + (-0.077)); // LP down

        // 500 <= pT < 600
        scaleFactorsHH_medPt_[0] = (1.072);            // HP central
        scaleFactorsHH_medPt_[1] = (1.072 + ( 0.041)); // HP up
        scaleFactorsHH_medPt_[2] = (1.072 + (-0.036)); // HP down
        scaleFactorsHH_medPt_[3] = (1.033);            // MP central
        scaleFactorsHH_medPt_[4] = (1.033 + ( 0.030)); // MP up
        scaleFactorsHH_medPt_[5] = (1.033 + (-0.025)); // MP down
        scaleFactorsHH_medPt_[6] = (1.006);            // LP central
        scaleFactorsHH_medPt_[7] = (1.006 + ( 0.024)); // LP up
        scaleFactorsHH_medPt_[8] = (1.006 + (-0.026)); // LP down

        // pT >= 600
        scaleFactorsHH_highPt_[0] = (1.046);            // HP central
        scaleFactorsHH_highPt_[1] = (1.046 + ( 0.038)); // HP up
        scaleFactorsHH_highPt_[2] = (1.046 + (-0.038)); // HP down
        scaleFactorsHH_highPt_[3] = (1.010);            // MP central
        scaleFactorsHH_highPt_[4] = (1.010 + ( 0.030)); // MP up
        scaleFactorsHH_highPt_[5] = (1.010 + (-0.035)); // MP down
        scaleFactorsHH_highPt_[6] = (1.001);            // LP central
        scaleFactorsHH_highPt_[7] = (1.001 + ( 0.035)); // LP up
        scaleFactorsHH_highPt_[8] = (1.001 + (-0.037)); // LP down
    }
    else
    {
        std::string errorMessage = "PNseSFInterface - Invalid period: "
                                 + period
                                 + " [options are: 2016preVFP/2016postVFP/2017/2018/2022/2023]";
        throw std::logic_error(errorMessage);
    }
}

// Fill SFs for TT-like samples
void PNetSFInterface::FillTTlikeSFs(const std::string period)
{
    // Set scale factor values depending on the Era
    if (period == "2016preVFP")
    {
        // pT < 300
        // HP and MP SFs not yet available - default to 1
        scaleFactorsTT_lowPt_[6] = (0.5942130);                 // LP central
        scaleFactorsTT_lowPt_[7] = (0.5942130 + ( 0.27098274)); // LP up
        scaleFactorsTT_lowPt_[8] = (0.5942130 + (-0.27098274)); // LP down

        // 300 <= pT < 400
        // HP and MP SFs not yet available - default to 1
        scaleFactorsTT_medPt_[6] = (0.76523553);                 // LP central
        scaleFactorsTT_medPt_[7] = (0.76523553 + ( 0.27098274)); // LP up
        scaleFactorsTT_medPt_[8] = (0.76523553 + (-0.27098274)); // LP down

        // pT >= 400
        // HP and MP SFs not yet available - default to 1
        scaleFactorsTT_highPt_[6] = (0.66079032);                 // LP central
        scaleFactorsTT_highPt_[7] = (0.66079032 + ( 0.37631591)); // LP up
        scaleFactorsTT_highPt_[8] = (0.66079032 + (-0.37631591)); // LP down
    }
    else if (period == "2016postVFP")
    {
        // Do nothing: SFs not yet available - default to 1 for all cases
    }
    else if (period == "2017")
    {
        // Do nothing: SFs not yet available - default to 1 for all cases
    }
    else if (period == "2018" || period == "2022" || period == "2023")
    {
        // Do nothing: SFs not yet available - default to 1 for all cases
    }
    else
    {
        std::string errorMessage = "PNseSFInterface - Invalid period: "
                                 + period
                                 + " [options are: 2016preVFP/2016postVFP/2017/2018/2022/2023]";
        throw std::logic_error(errorMessage);
    }
}

// Fill SFs for DY-like samples
void PNetSFInterface::FillDYlikeSFs(const std::string period)
{
    // Set scale factor values depending on the Era
    if (period == "2016preVFP")
    {
        // pT < 300
        // HP and MP SFs not yet available - default to 1
        scaleFactorsDY_lowPt_[6] = (0.79969335);                // LP central
        scaleFactorsDY_lowPt_[7] = (0.79969335 + ( 0.1880848)); // LP up
        scaleFactorsDY_lowPt_[8] = (0.79969335 + (-0.1880848)); // LP down

        // 300 <= pT < 400
        // HP and MP SFs not yet available - default to 1
        scaleFactorsDY_medPt_[6] = (1.27936957);                 // LP central
        scaleFactorsDY_medPt_[7] = (1.27936957 + ( 0.29244537)); // LP up
        scaleFactorsDY_medPt_[8] = (1.27936957 + (-0.29244537)); // LP down

        // pT >= 400
        scaleFactorsDY_highPt_[6] = (1.57300599);                 // LP central
        scaleFactorsDY_highPt_[7] = (1.57300599 + ( 0.47988435)); // LP up
        scaleFactorsDY_highPt_[8] = (1.57300599 + (-0.47988435)); // LP down
    }
    else if (period == "2016postVFP")
    {
        // Do nothing: SFs not yet available - default to 1 for all cases
    }
    else if (period == "2017")
    {
        // Do nothing: SFs not yet available - default to 1 for all cases
    }
    else if (period == "2018" || period == "2022" || period == "2023")
    {
        // Do nothing: SFs not yet available - default to 1 for all cases
    }
    else
    {
        std::string errorMessage = "PNseSFInterface - Invalid period: "
                                 + period
                                 + " [options are: 2016preVFP/2016postVFP/2017/2018/2022/2023]";
        throw std::logic_error(errorMessage);
    }
}

// Get the SF vector based on jet pT and sampleType
std::vector<float> PNetSFInterface::getSFvec(const ROOT::VecOps::RVec<float> FatJet_pt, const int fatjet_idx,
                                             const std::vector<bool> genAk8_Zbb_matches, const std::vector<bool> genAk8_Hbb_matches,
                                             const std::string sampleType)
{

    // Initialize output SF vector
    std::vector<float> SFvec(9, 1.);

    // Check that the AK8 jet exists
    if (fatjet_idx == -1)
    {
        return SFvec;
    }

    // "Other" backgrounds
    if (sampleType == "other")
    {
        return SFvec;
    }

    // Get pT of selected fatjet
    float pT = FatJet_pt.at(fatjet_idx);

    // HH-like samples SF
    if (sampleType == "HHlike")
    {
        // If the AK8 jet is not matched to a gen Z/H->bb resonance, return 1.
        if ( !genAk8_Zbb_matches.at(fatjet_idx) && !genAk8_Hbb_matches.at(fatjet_idx) )
        {
            return SFvec;
        }

        // Else set the pT-dependent SF
        if (pT < 500)
        {
            SFvec = scaleFactorsHH_lowPt_;
        }
        else if (pT >= 500 && pT < 600)
        {
            SFvec = scaleFactorsHH_medPt_;
        }
        else /* pT >= 600 */
        {
            SFvec = scaleFactorsHH_highPt_;
        }
    }

    // TT-like samples SF
    if (sampleType == "TTlike")
    {
        // Set the pT-dependent SF
        if (pT < 300)
        {
            SFvec = scaleFactorsTT_lowPt_;
        }
        else if (pT >= 300 && pT < 400)
        {
            SFvec = scaleFactorsTT_medPt_;
        }
        else /* pT >= 400 */
        {
            SFvec = scaleFactorsTT_highPt_;
        }
    }

    // DY-like samples SF
    if (sampleType == "DYlike")
    {
        // Set the pT-dependent SF
        if (pT < 300)
        {
            SFvec = scaleFactorsDY_lowPt_;
        }
        else if (pT >= 300 && pT < 400)
        {
            SFvec = scaleFactorsDY_medPt_;
        }
        else /* pT >= 400 */
        {
            SFvec = scaleFactorsDY_highPt_;
        }
    }

    // Return final SF vector
    return SFvec;
}
