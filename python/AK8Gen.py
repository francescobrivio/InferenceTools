import os
from analysis_tools.utils import import_root
from Base.Modules.baseModules import JetLepMetSyst

ROOT = import_root()

class AK8GenRDFProducer(JetLepMetSyst):
    def __init__(self, *args, **kwargs):
        super(AK8GenRDFProducer, self).__init__(self, *args, **kwargs)
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libToolsTools.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gSystem.Load("libToolsTools.so")

            base = "{}/{}/src/Tools/Tools".format(
                os.getenv("CMT_CMSSW_BASE"), os.getenv("CMT_CMSSW_VERSION"))
            ROOT.gROOT.ProcessLine(".L {}/interface/AK8GenInterface.h".format(base))
            ROOT.gInterpreter.Declare("""
                auto AK8Gen = AK8GenInterface();
            """)

    def run(self, df):
        if not self.isMC:
            return df, []

        variables = ["Ak8_Zbb_matches", "Ak8_Hbb_matches"]

        df = df.Define("ak8_gen_results", "AK8Gen.get_ak8_genmatch_info("
                            "GenPart_statusFlags, GenPart_pdgId, GenPart_status, "
                            "GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, "
                            "GenPart_phi, GenPart_mass, "
                            "FatJet_pt{0}, FatJet_eta, FatJet_phi, FatJet_mass{0})"
                            .format(self.jet_syst))

        branches = []
        for var in variables:
            df = df.Define(f"gen{var}", f"ak8_gen_results.{var}")
            branches.append(f"gen{var}")

        return df, branches


def AK8GenRDF(**kwargs):
    """
    Module to store genMatch information between AK8 and  H/Z->bb

    :param isMC: flag of the dataset being MC or data
    :type : bool

    YAML sintaxis:

    .. code-block:: yaml

        codename:
            name: AK8GenRDF
            path: Tools.Tools.AK8Gen
            parameters:
                isMC: self.dataset.process.isMC
    """
    return lambda: AK8GenRDFProducer(**kwargs)
