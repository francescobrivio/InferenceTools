#ifndef AK8GenInterface_h
#define AK8GenInterface_h

// ------------------------------------------------------------ //
//                                                              //
//   class AK8GenInterface                                      //
//                                                              //
//                                                              //
//   Author: Francesco Brivio (INFN Milano-Bicocca)             //
//   Date  : June 2024                                          //
//                                                              //
// ------------------------------------------------------------ //

#include <vector>
#include <ROOT/RVec.hxx>

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<int> iRVec;

struct gen_match_output {
    std::vector<bool> Ak8_Zbb_matches;
    std::vector<bool> Ak8_Hbb_matches;
};

class AK8GenInterface {

  public:
    AK8GenInterface();
    ~AK8GenInterface();

    gen_match_output get_ak8_genmatch_info(
        iRVec GenPart_statusFlags, iRVec GenPart_pdgId,
        iRVec GenPart_status, iRVec GenPart_genPartIdxMother,
        fRVec GenPart_pt, fRVec GenPart_eta, fRVec GenPart_phi, fRVec GenPart_mass,
        fRVec FatJet_pt, fRVec FatJet_eta, fRVec FatJet_phi, fRVec FatJet_mass);
};

#endif // AK8GenInterface_h
