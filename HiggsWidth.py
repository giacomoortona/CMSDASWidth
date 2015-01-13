from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the Higgs width

class Higgswidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.GGsmval = 1.
        self.poiMap = []
        self.pois = {}
        self.verbose = False
        self.is2l2nu = True
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        if process == "ggH_s": return "ggH_s_func"
        elif process == "ggH_b": return "ggH_b_func"
        elif process == "ggH_sbi": return "ggH_sbi_func"
        else: return 1
            
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("GGsmVal="):   
                self.GGsmval = po.replace("GGsmVal=","")
                print 'G/G_SM set to :', self.GGsmval
            if po == "GGsmfixed":
                print "Will fix CMS_zz4l_GGsm to 1 and float mu"
                self.GGsmfixed = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.is2l2nu:
            self.modelBuilder.doVar("CMS_widthH_kbkg[1.,0.,2.]")
            self.modelBuilder.doVar("R[1.,0.,4.]")
            self.modelBuilder.doVar("CMS_zz4l_GGsm[1.,0.,30.]")

        if self.GGsmfixed:
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setVal(1.)
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(True)
            self.modelBuilder.out.var("R").setVal(0.01,6.)
            self.modelBuilder.out.var("R").setConstant(False)
            print "Fixing CMS_zz4l_GGsm"
            poi = "R"
        else:
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setRange(0.0001,50.)
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setVal(float(self.GGsmval))
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(False)
            self.modelBuilder.out.var("R").setVal(1.)
            self.modelBuilder.out.var("R").setConstant(False)
            poi = "CMS_zz4l_GGsm"

        self.modelBuilder.factory_("expr::ggH_s_func(\"@0*@1-sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg)")
        self.modelBuilder.factory_("expr::ggH_b_func(\"@2-sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg)")
        self.modelBuilder.factory_("expr::ggH_sbi_func(\"sqrt(@0*@1*@2)\",R,CMS_zz4l_GGsm,CMS_widthH_kbkg)")

        
        self.modelBuilder.doSet("POI",poi)
        
higgswidth = Higgswidth()
