#include "Tools/Tools/interface/AK8GenInterface.h"
#include <iostream>
#include <algorithm>
#include <TLorentzVector.h>

// Constructor
AK8GenInterface::AK8GenInterface() {}

// Destructor
AK8GenInterface::~AK8GenInterface() {}

// Get ak8 match info
gen_match_output AK8GenInterface::get_ak8_genmatch_info(
    iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
    iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
    fRVec GenPart_pt, fRVec GenPart_eta, fRVec GenPart_phi, fRVec GenPart_mass,
    fRVec FatJet_pt, fRVec FatJet_eta, fRVec FatJet_phi, fRVec FatJet_mass)
{
    // Declare output
    gen_match_output genMatchOut;

    // Set all ouput to false
    for (size_t i = 0; i < FatJet_pt.size(); i++)
    {
        genMatchOut.Ak8_Zbb_matches.push_back(false);
        genMatchOut.Ak8_Hbb_matches.push_back(false);
    }

    // Declare vectors of gen H/Z->bb idxs
    std::vector<int> genHbb_idxs, genZbb_idxs;

    // Loop on gen particles to fill the genH and genZ idx vectors
    for (size_t idx = 0; idx < GenPart_pdgId.size(); idx++)
    {
        // Particle gen info
        int pdg    = GenPart_pdgId[idx];
        int mthIdx = GenPart_genPartIdxMother[idx];
        int mthPdg = (mthIdx > -1) ? GenPart_pdgId[mthIdx] : 0;

        bool isFirst       = (GenPart_statusFlags[idx] & (1<<12)) ? true : false;
        bool isHardProcess = (GenPart_statusFlags[idx] & (1<<7))  ? true : false;

        // If it's a b-quark
        if (abs(pdg) == 5 && isHardProcess && isFirst)
        {
            // Check whether it comes from H...
            if (abs(mthPdg) == 25)
            {
                // If motherH idx already saved (due to other b quark), continue
                if (std::find(genHbb_idxs.begin(), genHbb_idxs.end(), mthIdx) != genHbb_idxs.end())
                    continue;
                // else, save it
                else
                    genHbb_idxs.push_back(mthIdx);
            }
            // ...or Z
            else if (abs(mthPdg) == 23)
            {
                // If motherZ idx already saved (due to other b quark), continue
                if (std::find(genZbb_idxs.begin(), genZbb_idxs.end(), mthIdx) != genZbb_idxs.end())
                    continue;
                // else, save it
                else
                    genZbb_idxs.push_back(mthIdx);
            }
        }
    }

    // Now loop on FatJets and do the geometric matching (dR < 0.8) with the gen H/Z->bb
    for (size_t i = 0; i < FatJet_pt.size(); i++)
    {
        // Reco ak8 jet
        TLorentzVector fatjet_tlv;
        fatjet_tlv.SetPtEtaPhiM(FatJet_pt[i], FatJet_eta[i], FatJet_phi[i], FatJet_mass[i]);

        // Matching with H->bb
        for (auto idx : genHbb_idxs)
        {
            TLorentzVector gen_Hbb;
            gen_Hbb.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
            if (fatjet_tlv.DeltaR(gen_Hbb) < 0.8)
                genMatchOut.Ak8_Hbb_matches[i] = true;
        }

        // Matching with Z->bb
        for (auto idx : genZbb_idxs)
        {
            TLorentzVector gen_Zbb;
            gen_Zbb.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
            if (fatjet_tlv.DeltaR(gen_Zbb) < 0.8)
                genMatchOut.Ak8_Zbb_matches[i] = true;
        }
    }

    return genMatchOut;
}
