#include "Tools/Tools/interface/HHJetsInterface.h"

// Constructor
HHJetsInterface::HHJetsInterface (std::string model_0, std::string model_1, int year):
  HHbtagger_(std::array<std::string, 2> { {model_0, model_1} })
{
  year_ = year;
}


// Destructor
HHJetsInterface::~HHJetsInterface() {}



std::vector<float> HHJetsInterface::GetHHJets(
    unsigned long long int event, int pairType,
    fRVec Jet_pt, fRVec Jet_eta, fRVec Jet_phi, fRVec Jet_mass,
    iRVec Jet_puId, fRVec Jet_jetId, fRVec Jet_btagDeepFlavB,
    float dau1_pt, float dau1_eta, float dau1_phi, float dau1_mass,
    float dau2_pt, float dau2_eta, float dau2_phi, float dau2_mass,
    float met_pt, float met_phi)
{
  std::vector<float> all_HHbtag_scores;
  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    all_HHbtag_scores.push_back(-999.);
  }
  
  auto dau1_tlv = TLorentzVector();
  auto dau2_tlv = TLorentzVector();
  auto met_tlv = TLorentzVector();
  dau1_tlv.SetPtEtaPhiM(dau1_pt, dau1_eta, dau1_phi, dau1_mass);
  dau2_tlv.SetPtEtaPhiM(dau2_pt, dau2_eta, dau2_phi, dau2_mass);
  met_tlv.SetPxPyPzE(met_pt * cos(met_phi), met_pt * sin(met_phi), 0, met_pt);

  std::vector <jet_idx_btag> jet_indexes;
  std::vector <int> all_jet_indexes;
  for (size_t ijet = 0; ijet < Jet_pt.size(); ijet++) {
    if ((Jet_puId[ijet] < 4 && Jet_pt[ijet] <= 50) || Jet_jetId[ijet] < 2)
      continue;
    auto jet_tlv = TLorentzVector();
    jet_tlv.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], Jet_mass[ijet]);
    if (jet_tlv.DeltaR(dau1_tlv) < 0.5 || jet_tlv.DeltaR(dau2_tlv) < 0.5)
      continue;
    if (Jet_pt[ijet] > 20 && fabs(Jet_eta[ijet]) < 2.4)
      jet_indexes.push_back(jet_idx_btag({ijet, Jet_btagDeepFlavB[ijet]}));
    if (Jet_pt[ijet] > 20 && fabs(Jet_eta[ijet]) < 4.7)
      all_jet_indexes.push_back(ijet);
  }
  if (jet_indexes.size() >= 2) {
    std::stable_sort(jet_indexes.begin(), jet_indexes.end(), jetSort);

    auto htt_tlv = dau1_tlv + dau2_tlv;
    std::vector <float> HHbtag_jet_pt_, HHbtag_jet_eta_, HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_;
    std::vector <float> HHbtag_jet_htt_deta_, HHbtag_jet_htt_dphi_, HHbtag_jet_deepFlavour_;

    for (auto &jet : jet_indexes) {
      auto jet_tlv = TLorentzVector();
      jet_tlv.SetPtEtaPhiM(Jet_pt[jet.idx], Jet_eta[jet.idx], Jet_phi[jet.idx], Jet_mass[jet.idx]);
      HHbtag_jet_pt_.push_back(jet_tlv.Pt());
      HHbtag_jet_eta_.push_back(jet_tlv.Eta());
      HHbtag_rel_jet_M_pt_.push_back(jet_tlv.M() / jet_tlv.Pt());
      HHbtag_rel_jet_E_pt_.push_back(jet_tlv.E() / jet_tlv.Pt());
      HHbtag_jet_htt_deta_.push_back(htt_tlv.Eta() - jet_tlv.Eta());
      HHbtag_jet_htt_dphi_.push_back(ROOT::Math::VectorUtil::DeltaPhi(htt_tlv, jet_tlv));
      HHbtag_jet_deepFlavour_.push_back(jet.btag);
    }
    
    auto HHbtag_htt_met_dphi_ = (float) ROOT::Math::VectorUtil::DeltaPhi(htt_tlv, met_tlv);
    auto HHbtag_htt_scalar_pt_ = (float) (dau1_tlv.Pt() + dau2_tlv.Pt());
    auto HHbtag_rel_met_pt_htt_pt_ = (float) met_tlv.Pt() / HHbtag_htt_scalar_pt_;
    auto HHbtag_htt_pt_ = (float) htt_tlv.Pt();
    auto HHbtag_htt_eta_ = (float) htt_tlv.Eta();
    int HHbtag_channel_ = -1;

    if (pairType == 0)
      HHbtag_channel_ = 1;
    else if (pairType == 1)
      HHbtag_channel_ = 0;
    else if (pairType == 2)
      HHbtag_channel_ = 1;

    auto HHbtag_scores = HHbtagger_.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_,
      HHbtag_rel_jet_M_pt_, HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_,
      HHbtag_jet_deepFlavour_, HHbtag_jet_htt_dphi_, year_, HHbtag_channel_,
      HHbtag_htt_pt_, HHbtag_htt_eta_, HHbtag_htt_met_dphi_,
      HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, event);

    for (size_t ijet = 0; ijet < jet_indexes.size(); ijet++) {
      all_HHbtag_scores[jet_indexes[ijet].idx] = HHbtag_scores[ijet];
    }
  }
  return all_HHbtag_scores;
}


// GetScore
std::vector<float> HHJetsInterface::GetScore(
  std::vector<float> HHbtag_jet_pt_, std::vector<float> HHbtag_jet_eta_, std::vector<float> HHbtag_rel_jet_M_pt_,
  std::vector<float> HHbtag_rel_jet_E_pt_, std::vector<float> HHbtag_jet_htt_deta_, std::vector<float> HHbtag_jet_deepFlavour_,
  std::vector<float> HHbtag_jet_htt_dphi_, int HHbtag_year_, int HHbtag_channel_, float HHbtag_tauH_pt_, float HHbtag_tauH_eta_,
  float HHbtag_htt_met_dphi_, float HHbtag_rel_met_pt_htt_pt_, float HHbtag_htt_scalar_pt_, unsigned long long int HHbtag_evt_)
{
  // Get HHbtag score
  auto HHbtag_scores = HHbtagger_.GetScore(HHbtag_jet_pt_, HHbtag_jet_eta_, HHbtag_rel_jet_M_pt_,
      HHbtag_rel_jet_E_pt_, HHbtag_jet_htt_deta_, HHbtag_jet_deepFlavour_, HHbtag_jet_htt_dphi_,
      HHbtag_year_, HHbtag_channel_, HHbtag_tauH_pt_, HHbtag_tauH_eta_, HHbtag_htt_met_dphi_,
      HHbtag_rel_met_pt_htt_pt_, HHbtag_htt_scalar_pt_, HHbtag_evt_);

  // Store HHbtag scores in a map<jet_idx,HHbtag_score>
  std::vector<float> jets_and_HHbtag;
  for (unsigned int i = 0; i < HHbtag_jet_pt_.size(); i++)
  {
    jets_and_HHbtag.push_back(HHbtag_scores.at(i));
  }

  return jets_and_HHbtag;
}
