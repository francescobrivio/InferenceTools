from Corrections.JME.PUjetID_SF import PUjetID_SFRDFProducer
from Corrections.BTV.btag_SF import btag_SFRDFProducer


class HHPUjetID_SFRDFProducer(PUjetID_SFRDFProducer):
    def __init__(self, year, *args, **kwargs):
        super(HHPUjetID_SFRDFProducer, self).__init__(year, *args, **kwargs)
        self.lep_pt = "{dau1_pt%s, dau2_pt%s}" % (self.systs, self.systs)
        self.lep_eta = "{dau1_eta, dau2_eta}"
        self.lep_phi = "{dau1_phi, dau2_phi}"
        self.lep_mass = "{dau1_mass%s, dau2_mass%s}" % (self.systs, self.systs)


def HHPUjetID_SFRDF(**kwargs):
    """
    Module to compute PU Jet Id scale factors for the HH->bbtautau analysis.
    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHPUjetID_SFRDF
            path: Tools.Tools.HH_SF
            parameters:
                year: self.config.year
                isMC: self.dataset.process.isMC
                isUL: self.dataset.has_tag('ul')
                ispreVFP: self.config.get_aux("isPreVFP", False)

    """
    year = kwargs.pop("year")
    return lambda: HHPUjetID_SFRDFProducer(year, **kwargs)


class HHbtag_SFRDFProducer(btag_SFRDFProducer):
    def __init__(self, year, *args, **kwargs):
        super(HHbtag_SFRDFProducer, self).__init__(year, *args, **kwargs)
        self.lep_pt = "{dau1_pt%s, dau2_pt%s}" % (self.systs, self.systs)
        self.lep_eta = "{dau1_eta, dau2_eta}"
        self.lep_phi = "{dau1_phi, dau2_phi}"
        self.lep_mass = "{dau1_mass%s, dau2_mass%s}" % (self.systs, self.systs)


def HHbtag_SFRDF(**kwargs):
    """
    Module to obtain btagging deepJet SFs with their uncertainties
    for the HH->bbtautau analysis.

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HHbtag_SFRDF
            path: Tools.Tools.HH_SF
            parameters:
                isMC: self.dataset.process.isMC
                year: self.config.year
                isUL: self.dataset.has_tag('ul')
                reshape_uncertainties: [central, ...]

    """
    year = kwargs.pop("year")
    return lambda: HHbtag_SFRDFProducer(year, **kwargs)

class HH_pNetSFRDFProduced(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(HH_pNetSFRDFProduced, self).__init__(*args, **kwargs)
        self.isMC = kwargs.pop("isMC")
        self.year = kwargs.pop("year")

        # Check whether HH-like, TT-like, DY-like or other for pNet SF:
        #  - HH-like : ggHH, VBFHH, ZH, WH, ttH, ggH, qqH, ttWH, ttZH
        #  - TT-like : TT fullyLep, semiLep, fullyHad
        #  - DY_like : all DY samples
        #  - other   : all other samples
        dsetName = kwargs.pop("dataset")
        if "HH" in dsetName or "Hto2B" in dsetName or "WH" in dsetName or "ZH" in dsetName:
            sampleType = "HHlike"
        elif "TTto" in dsetName:
            sampleType = "TTlike"
        elif "DY" in dsetName:
            sampleType = "DYlike"
        else:
            sampleType = "other"
        self.sampleType = sampleType

        if self.isMC:
            base = "{}/{}/src/Tools/Tools".format(os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/PNetSFInterface.h".format(base))

            ROOT.gInterpreter.Declare("""
                auto PNetAK8SF = PNetSFInterface("%s");
            """ % (self.year))

    def run(self, df):
        # In case of data: no SFs
        if not self.isMC:
            return df, []

        # In case of MC: compute SFs based on jet pT, sample type and gen-matches variables
        df = df.Define("fatjet_pNet_SF_vec", "PNetAK8SF.getSFvec(FatJet_pt{1}, fatjet_JetIdx, "
                                             " genAk8_Zbb_matches, genAk8_Hbb_matches, "
                                             " \"{0}\")".format(self.sampleType, self.jet_syst))
        df = df.Define("fatjet_pNet_HP_SF",      "fatjet_pNet_SF_vec[0]")
        df = df.Define("fatjet_pNet_HP_SF_up",   "fatjet_pNet_SF_vec[1]")
        df = df.Define("fatjet_pNet_HP_SF_down", "fatjet_pNet_SF_vec[2]")
        df = df.Define("fatjet_pNet_MP_SF",      "fatjet_pNet_SF_vec[3]")
        df = df.Define("fatjet_pNet_MP_SF_up",   "fatjet_pNet_SF_vec[4]")
        df = df.Define("fatjet_pNet_MP_SF_down", "fatjet_pNet_SF_vec[5]")
        df = df.Define("fatjet_pNet_LP_SF",      "fatjet_pNet_SF_vec[6]")
        df = df.Define("fatjet_pNet_LP_SF_up",   "fatjet_pNet_SF_vec[7]")
        df = df.Define("fatjet_pNet_LP_SF_down", "fatjet_pNet_SF_vec[8]")

        return df, ["fatjet_pNet_HP_SF", "fatjet_pNet_HP_SF_up", "fatjet_pNet_HP_SF_down",
                    "fatjet_pNet_MP_SF", "fatjet_pNet_MP_SF_up", "fatjet_pNet_MP_SF_down",
                    "fatjet_pNet_LP_SF", "fatjet_pNet_LP_SF_up", "fatjet_pNet_LP_SF_down"]


def HH_pNetSFRDF(**kwargs):
    """
    Module to obtain pNet SFs for AK8 jets depending on sample type,
    FataJet pT and gen-matching to a H/Z->bb resonance

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: HH_pNetSFRDF
            path: Tools.Tools.HH_SF
            parameters:
                isMC: self.dataset.process.isMC
                dataset: self.dataset.name
                year: self.config.year

    """
    return lambda: HH_pNetSFRDFProduced(**kwargs)
