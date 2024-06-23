import os

from analysis_tools.utils import import_root

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from Tools.Tools.tau_utils import LeptonTauPair, TriggerChecker, lepton_veto
from Base.Modules.baseModules import JetLepMetSyst, JetLepMetModule

ROOT = import_root()

class HHLeptonRDFProducer(JetLepMetSyst):
    def __init__(self, isMC, year, runEra, pairType_filter, *args, **kwargs):
        super(HHLeptonRDFProducer, self).__init__(isMC=isMC, *args, **kwargs)
        self.isMC = isMC
        self.year = year
        self.runEra = runEra
        self.isRun3 = kwargs.pop("isRun3", False)
        self.pairType_filter = pairType_filter
        self.deeptau_version = kwargs.pop("deeptau_version", "2017v2p1")
        self.isV10 = kwargs.pop("isV10", False)
        vvvl_vsjet = kwargs.pop("vvvl_vsjet")
        vl_vse = kwargs.pop("vl_vse")
        vvl_vse = kwargs.pop("vvl_vse")
        t_vsmu = kwargs.pop("t_vsmu")
        vl_vsmu = kwargs.pop("vl_vsmu")
        self.skip_etau_ele_off = kwargs.pop("skip_etau_ele_off", 0)
        self.skip_etau_tau_off = kwargs.pop("skip_etau_tau_off", 0)
        self.skip_etau_dR = kwargs.pop("skip_etau_dR", 0)
        self.skip_etau_trg = kwargs.pop("skip_etau_trg", 0)
        self.skip_etau_veto = kwargs.pop("skip_etau_veto", 0)
        self.skip_mutau_mu_off = kwargs.pop("skip_mutau_mu_off", 0)
        self.skip_mutau_tau_off = kwargs.pop("skip_mutau_tau_off", 0)
        self.skip_mutau_dR = kwargs.pop("skip_mutau_dR", 0)
        self.skip_mutau_trg = kwargs.pop("skip_mutau_trg", 0)
        self.skip_mutau_veto = kwargs.pop("skip_mutau_veto", 0)
        self.skip_tautau_tau_off = kwargs.pop("skip_tautau_tau_off", 0)
        self.skip_tautau_dR = kwargs.pop("skip_tautau_dR", 0)
        self.skip_tautau_trg = kwargs.pop("skip_tautau_trg", 0)
        self.skip_tautau_veto = kwargs.pop("skip_tautau_veto", 0)
        self.skip_btautau_tau_off = kwargs.pop("skip_btautau_tau_off", 0)
        self.skip_btautau_trg = kwargs.pop("skip_btautau_trg", 0)

        if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gSystem.Load("libToolsTools.so")

        # DOCUMENTATION
        # The list of triggers here is checked against the list of branches in NanoAOD, non-existent ones are replaced by false.
        # A "Vbool triggers" array is built with the HLT_* trigger paths and passed to get_*tau_triggers functions
        # Then a trig_req object is built. pass = HLT_* branch (or false in case it does not exist). 
        # struct trig_req {bool pass; float pt1; float eta1; float pt2; float eta2; std::vector<std::vector<int>> bits; };
        # the values here are used to make cuts on *offline* objects, nb1 is electron/muon/hadTau, nb2 is hadTau

        self.mutau_triggers = ["HLT_IsoMu22", "HLT_IsoMu22_eta2p1",
            "HLT_IsoTkMu24", "HLT_IsoTkMu22_eta2p1", "HLT_IsoMu24", "HLT_IsoMu27",
            "HLT_IsoMu19_eta2p1_LooseIsoPFTau20", "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
            "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1",
            "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1"]
        self.etau_triggers = ["HLT_Ele25_eta2p1_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
            "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf",
            "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1",
            "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1"]
        self.tautau_triggers = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg",
            "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg",
            "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg",
            "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg",
            "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1"]
        self.boostedtau_triggers = ["HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1",
            "HLT_AK8PFJet400_TrimMass30"]
        self.tautaujet_triggers = ["HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60"]
        self.vbf_triggers = ["HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg",
            "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1",
            "HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1"]

        if not os.getenv("_HHLepton"):
            os.environ["_HHLepton"] = "_HHLepton"

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/HHLeptonInterface.h".format(base))
            ROOT.gInterpreter.Declare("""
                auto HHLepton = HHLeptonInterface(%s, %s, %s, %s, %s);
            """ % (vvvl_vsjet, vl_vse, vvl_vse, t_vsmu, vl_vsmu))

            if self.year == 2018:
                ROOT.gInterpreter.Declare("""
                    using Vbool = const ROOT::RVec<Bool_t>&;
                    std::vector<trig_req> get_mutau_triggers(
                            Vbool triggers, bool isMC, int run, int runEra) {
                        std::vector<trig_req> trigger_reqs;
                        trigger_reqs.push_back(trig_req({triggers[4], 26, 2.4, 20, 2.3, 24, 0, {{2, 8}, {}}})); // HLT_IsoMu24 (not prescaled for 2018)
                        // trigger_reqs.push_back(trig_req({triggers[5], 29, 2.4, 20, 2.3, 27, 0, {{2, 8}, {}}})); // HLT_IsoMu27 -> not to be used for 2018
                        // https://twiki.cern.ch/twiki/bin/view/CMS/TauTrigger
                        // The TWiki suggests using HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 rather than HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 but our SFs are presumably for the latter (trigger used for gamma gamma ->tautau)
                        // Muon filter bits : 2(->4): *OverlapFilterIsoMu*PFTau*, 6(->64): hlt*OverlapFilterIsoMu*PFTau* (the 2 seem identical ?)
                        // Tau filter bits : 0(->1): LooseChargedIso, 5(->32): *Hps*, 9(->512): hlt*OverlapFilterIsoMu*PFTau* (mutau)
                        // Requiring HPS is not in  Tau POG wiki and probably is not needed as all paths matching hlt*OverlapFilterIsoMu*PFTau* are HPS when that is available
                        // Requiring LooseChargedIso is also not in Tau POG twiki
                        if (!isMC && run < 317509) {
                            trigger_reqs.push_back(trig_req({triggers[8], 22, 2.1, 32, 2.1, 0, 27, {{}, {}}})); // HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1
                        }
                        else {
                            // used to be {1, 32} for tau leg which is wrong I think ? that would require every tau LooseIso HPS trigger
                            trigger_reqs.push_back(trig_req({triggers[9], 22, 2.1, 32, 2.1, 0, 27, {{}, {}}})); // HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1
                        }
                        return trigger_reqs;
                    }
                    std::vector<trig_req> get_etau_triggers(
                            Vbool triggers, bool isMC, int run, int runEra) {
                        std::vector<trig_req> trigger_reqs;
                        trigger_reqs.push_back(trig_req({triggers[2], 33, 2.5, 20, 2.3, 32, 0, {{2}, {}}})); // HLT_Ele32_WPTight_Gsf
                        // trigger_reqs.push_back(trig_req({triggers[3], 36, 2.1, 20, 2.3, {{2}, {}}})); // HLT_Ele35_WPTight_Gsf
                        // gg->tautau analysis does OR with HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 but extremly few events pass this and not our trigger (<0.01% in ttbar)
                        if (!isMC && run < 317509) {
                            trigger_reqs.push_back(trig_req({triggers[4], 25, 2.1, 35, 2.1, 0, 0, {{}, {}}})); // HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1
                        }
                        else {
                            trigger_reqs.push_back(trig_req({triggers[5], 25, 2.1, 35, 2.1, 0, 0, {{}, {}}})); // HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1
                        }
                        return trigger_reqs;
                    }
                    std::vector<trig_req> get_tautau_triggers(
                            Vbool triggers, bool isMC, bool isRun3, int run, int runEra) {
                        std::vector<trig_req> trigger_reqs;
                        // Tau POG doc for ditau : (same as 2017)
                        // TrigObj_id == 15 && (TrigObj_filterBits&64)!=0 && (
                        //        ( (TrigObj_filterBits&4)!=0 && (TrigObj_filterBits&16)!=0 ) ||
                        //        ( TrigObj_pt > 40 && ( (TrigObj_filterBits&2)!=0 && (TrigObj_filterBits&16)!=0 ) 
                        //    || (TrigObj_filterBits&4)!=0 ) ) 
                        // Summary : (4&16&64) || (2&16&64&Trig_pt>40) || (4&64)
                        // Filter bits : 1(->2):MediumChargedIso, 2(->4):TightChargedIso, 4(->16):TightID OOSC photons (*TightOOSCPhotons*), 6(->64):charged iso di-tau (hlt*DoublePFTau*TrackPt1*ChargedIsolation*Dz02*)
                        if (!isMC && run < 317509) {
                            trigger_reqs.push_back(trig_req({triggers[2], 40, 2.1, 40, 2.1, 0, 0, {{}, {}}})); // HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg
                            // the offline req is indeed 40 GeV according to tau pog twiki, + need to put online 40 GeV req.
                            trigger_reqs.push_back(trig_req({triggers[3], 40, 2.1, 40, 2.1, 40, 40, {{}, {}}})); // HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg
                            trigger_reqs.push_back(trig_req({triggers[4], 40, 2.1, 40, 2.1, 0, 0, {{}, {}}})); // HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg
                        }
                        else {
                            trigger_reqs.push_back(trig_req({triggers[5], 40, 2.1, 40, 2.1, 0, 0, {{}, {}}})); // HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg
                        }
                        return trigger_reqs;
                    }
                    std::vector<trig_req> get_boosted_tau_triggers(
                            Vbool triggers, bool isMC, bool isRun3, int run, int runEra) {
                        std::vector<trig_req> trigger_reqs;
                        trigger_reqs.push_back(trig_req({triggers[0], 185, 2.1, 20, 2.1, 180, 0, {{}, {}}})); // HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1
                        return trigger_reqs;
                    }
                    std::vector<trig_req> get_tautaujet_triggers(Vbool triggers, bool isRun3) {
                        std::vector<trig_req> trigger_reqs;
                        return trigger_reqs;
                    }
                    std::vector<trig_req> get_vbf_triggers(
                            Vbool triggers, bool isMC, int run, int runEra) {
                        std::vector<trig_req> trigger_reqs;
                        // if (!isMC && run < 317509)
                        //     trigger_reqs.push_back(trig_req({triggers[1], 25, 2.1, 25, 2.1, {{1}, {1}}})); // HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1
                        // else
                        //     trigger_reqs.push_back(trig_req({triggers[2], 25, 2.1, 25, 2.1, {{1, 32}, {1, 32}}})); // HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1
                        return trigger_reqs;
                    }
                """)

    def run(self, df):
        variables = ["pairType", "dau1_index", "dau2_index",
            "isTauTauJetTrigger", "isVBFtrigger", "isOS",
            "dau1_eta", "dau1_phi", "dau1_iso", "dau1_decayMode",
            "dau1_idDeepTauVSe", "dau1_idDeepTauVSmu",
            "dau1_idDeepTauVSjet",
            "dau2_eta", "dau2_phi", "dau2_decayMode",
            "dau2_idDeepTauVSe", "dau2_idDeepTauVSmu",
            "dau2_idDeepTauVSjet"
        ]                  

        all_branches = df.GetColumnNames()
        for ib, branch in enumerate(self.mutau_triggers):
            if branch not in all_branches:
                self.mutau_triggers[ib] = "false"
        for ib, branch in enumerate(self.etau_triggers):
            if branch not in all_branches:
                self.etau_triggers[ib] = "false"
        for ib, branch in enumerate(self.tautau_triggers):
            if branch not in all_branches:
                self.tautau_triggers[ib] = "false"
        for ib, branch in enumerate(self.tautaujet_triggers):
            if branch not in all_branches:
                self.tautaujet_triggers[ib] = "false"
        for ib, branch in enumerate(self.vbf_triggers):
            if branch not in all_branches:
                self.vbf_triggers[ib] = "false"

        runEras = ["dum", "A", "B", "C", "D", "E", "F", "G", "H"]
        runEra = None
        for _irun, _runEra in enumerate(runEras):
            if self.runEra == _runEra:
                runEra = _irun
                break
        assert runEra != None

        skip = ""
        is_cutflow = (self.skip_etau_ele_off or self.skip_etau_tau_off or self.skip_etau_dR or self.skip_etau_trg or self.skip_etau_veto or \
                self.skip_mutau_mu_off or self.skip_mutau_tau_off or self.skip_mutau_dR or self.skip_mutau_trg or self.skip_mutau_veto or \
                self.skip_tautau_tau_off or self.skip_tautau_dR or self.skip_tautau_trg or self.skip_tautau_veto or \
                self.skip_btautau_tau_off or self.skip_btautau_trg)
        if is_cutflow:
            skip = "_skip" + self.skip_etau_ele_off * "_ETau_eleOff" + self.skip_etau_tau_off * "_ETau_tauOff" + \
                self.skip_etau_dR * "_ETau_dR" + self.skip_etau_trg * "_ETau_Trg" + self.skip_etau_veto * "_ETau_LepVeto" + \
                self.skip_mutau_mu_off * "_MuTau_muOff" + self.skip_mutau_tau_off * "_MuTau_tauOff" + \
                self.skip_mutau_dR * "_MuTau_dR" + self.skip_mutau_trg * "_MuTau_Trg" + self.skip_mutau_veto * "_MuTau_LepVeto" + \
                self.skip_tautau_tau_off * "_TauTau_tauOff" + self.skip_tautau_dR * "_TauTau_dR" + \
                self.skip_tautau_trg * "_TauTau_Trg" + self.skip_tautau_veto * "_TauTau_LepVeto" + \
                self.skip_btautau_tau_off * "_BTauTau_tauOff" + self.skip_btautau_trg * "_BTauTau_Trg"

        df = df.Define(f"mutau_triggers{skip}", "get_mutau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.mutau_triggers), ("true" if self.isMC else "false"), runEra))
        df = df.Define(f"etau_triggers{skip}", "get_etau_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.etau_triggers), ("true" if self.isMC else "false"), runEra))
        df = df.Define(f"tautau_triggers{skip}", "get_tautau_triggers({%s}, %s, %s, run, %s)" % (
            ", ".join(self.tautau_triggers), ("true" if self.isMC else "false"),
            ("true" if self.isRun3 else "false"), runEra))
        df = df.Define(f"boosted_tau_triggers{skip}", "get_boosted_tau_triggers({%s}, %s, %s, run, %s)" % (
            ", ".join(self.boostedtau_triggers), ("true" if self.isMC else "false"),
            ("true" if self.isRun3 else "false"), runEra))
        df = df.Define(f"tautaujet_triggers{skip}", "get_tautaujet_triggers({%s}, %s)" % (
            ", ".join(self.tautaujet_triggers), ("true" if self.isRun3 else "false")))
        df = df.Define(f"vbf_triggers{skip}", "get_vbf_triggers({%s}, %s, run, %s)" % (
            ", ".join(self.vbf_triggers), ("true" if self.isMC else "false"), runEra))

        Electron_mvaIso_WP80 = "Electron_mvaIso_WP80"
        if Electron_mvaIso_WP80 not in all_branches:
            Electron_mvaIso_WP80 = "Electron_mvaFall17V2Iso_WP80"
        Electron_mvaIso_WP90 = "Electron_mvaIso_WP90"
        if Electron_mvaIso_WP90 not in all_branches:
            Electron_mvaIso_WP90 = "Electron_mvaFall17V2Iso_WP90"
        Electron_mvaNoIso_WP90 = "Electron_mvaNoIso_WP90"
        if Electron_mvaNoIso_WP90 not in all_branches:
            Electron_mvaNoIso_WP90 = "Electron_mvaFall17V2noIso_WP90"

        df = df.Define(f"hh_lepton_results{skip}", "HHLepton.get_dau_indexes("
            "Muon_pt{0}, Muon_eta, Muon_phi, Muon_mass{0}, "
            "Muon_pfRelIso04_all, Muon_dxy, Muon_dz, Muon_mediumId, Muon_tightId, Muon_charge, "
            "Electron_pt{1}, Electron_eta, Electron_phi, Electron_mass{1}, "
            "{3}, {4}, {5}, Electron_pfRelIso03_all, "
            "Electron_dxy, Electron_dz, Electron_charge, "
            "Tau_pt{2}, Tau_eta, Tau_phi, Tau_mass{2}, "
            "Tau_idDeepTau{6}VSmu, Tau_idDeepTau{6}VSe, "
            "Tau_idDeepTau{6}VSjet, Tau_rawDeepTau{6}VSjet, "
            "Tau_dz, Tau_decayMode, Tau_charge, "
            "TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, "
            "mutau_triggers{7}, etau_triggers{7}, tautau_triggers{7}, tautaujet_triggers{7}, vbf_triggers{7}, "
            "{8}, {9}, {10}, {11}, {12}, "
            "{13}, {14}, {15}, {16}, {17}, "
            "{18}, {19}, {20}, {21}"
        ")".format(self.muon_syst, self.electron_syst, self.tau_syst,
            Electron_mvaIso_WP80, Electron_mvaNoIso_WP90, Electron_mvaIso_WP90,
            self.deeptau_version, skip, 
            self.skip_etau_ele_off, self.skip_etau_tau_off, self.skip_etau_dR, self.skip_etau_trg, self.skip_etau_veto,
            self.skip_mutau_mu_off, self.skip_mutau_tau_off, self.skip_mutau_dR, self.skip_mutau_trg, self.skip_mutau_veto,
            self.skip_tautau_tau_off, self.skip_tautau_dR, self.skip_tautau_trg, self.skip_tautau_veto))

        df = df.Define(f"hh_boosted_tau_results{skip}", "HHLepton.get_boosted_tau_indexes("
            # "boostedTau_pt{0}, boostedTau_eta, boostedTau_phi, boostedTau_mass{0}, " # [FIXME] no systs for the moment
            "boostedTau_pt, boostedTau_eta, boostedTau_phi, boostedTau_mass, "
            "boostedTau_idAntiEle2018, boostedTau_idAntiMu, "
            "boostedTau_idMVAnewDM2017v2, boostedTau_charge, "
            "TrigObj_id, TrigObj_filterBits, TrigObj_pt, TrigObj_eta, TrigObj_phi, "
            "boosted_tau_triggers{0}, 15, "
            "{1}, {2}"
        ")".format(
            # self.tau_syst # [FIXME] no systs for the moment
            skip, self.skip_btautau_tau_off, self.skip_btautau_trg))

        branches = []
        for var in variables:
            branchName = var
            if "DeepTau" in branchName:
                branchName = var[:var.index("VS")] + self.deeptau_version + var[var.index("VS"):]
            df = df.Define(branchName + skip, "hh_lepton_results%s.pairType >= 0 ? hh_lepton_results%s.%s : hh_boosted_tau_results%s.%s" % (skip, skip, var, skip, var))
            branches.append(branchName + skip)

        if self.pairType_filter:
            df = df.Filter(f"pairType{skip} >= 0", f"HHLeptonRDF{skip}")

        return df, branches

def HHLeptonRDF(**kwargs):
    """
    Returns the index of the two selected taus + several of their variables not affected by
    systematics.

    Lepton systematics (used for pt and mass variables) can be modified using the parameters from 
    :ref:`BaseModules_JetLepMetSyst`.

    :param runEra: run period in caps (data only)
    :type runEra: str

    :param isV10: whether the input sample is from nanoaodV10 (default: ``False``)
    :type isV10: bool

    :param deeptau_version: version of the DeepTau discriminator (default: ``2017v2p1``)
    :type deeptau_version: str

    :param vvvl_vsjet: VVVLoose DeepTauVSjet WP value
    :type vvvl_vsjet: int

    :param vl_vse: VLoose DeepTauVSe WP value
    :type vl_vse: int

    :param vvl_vse: VVLoose DeepTaVSe WP value
    :type vvl_vse: int

    :param t_vsmu: Tight DeepTauVSmu WP value
    :type t_vsmu: int

    :param vl_vsmu: VLoose DeepTauVSmu WP value
    :type vl_vsmu: int

    :param filter: whether to filter out output events if they don't have 2 lepton candidates
    :type filter: bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHLeptonRDF
            path: Tools.Tools.HHLeptonDev
            parameters:
                isMC: self.dataset.process.isMC
                isV10: self.dataset.has_tag("nanoV10")
                year: self.config.year
                runEra: self.dataset.runEra
                runEra: self.dataset.runEra
                vvvl_vsjet: self.config.deeptau.vsjet.VVVLoose
                vl_vse: self.config.deeptau.vse.VLoose
                vvl_vse: self.config.deeptau.vse.VVLoose
                t_vsmu: self.config.deeptau.vsmu.Tight
                vl_vsmu: self.config.deeptau.vsmu.VLoose
                pairType_filter: True

    """
    pairType_filter = kwargs.pop("pairType_filter", False)
    return lambda: HHLeptonRDFProducer(pairType_filter=pairType_filter, **kwargs)
