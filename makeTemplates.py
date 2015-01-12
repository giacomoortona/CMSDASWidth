#! /usr/bin/env python
import sys, pwd, commands
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array

def makeTemplates():

    #includes
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
    ROOT.gSystem.AddIncludePath("-Iinclude/")
    ROOT.gSystem.Load("libRooFit")
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    #Get inputs
    fileS = TFile.Open("ZZ4lAnalysis_Sig.root")
    fileB = TFile.Open("ZZ4lAnalysis_Bkg.root")
    fileSBI = TFile.Open("ZZ4lAnalysis_All.root")

    fileTemplates = TFile.Open("/afs/cern.ch/work/u/usarica/public/CombineTemplates/02_04_2014/8TeV/2mu2e/220/HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_Nominal.root")

    tempS = fileTemplates.Get("T_2D_1").Clone("tempS")
    tempB = fileTemplates.Get("T_2D_2").Clone("tempB")
    tempI = fileTemplates.Get("T_2D_4").Clone("tempI")

    rate_signal_ggzz_Shape = tempS.Integral()
    rate_bkg_ggzz_Shape = tempB.Integral()
    rate_interf_ggzz_Shape = tempI.Integral()

    #Create histograms
    h2S = tempI.Clone("h2S")
    h2B = tempI.Clone("h2B")
    h2I = tempI.Clone("h2SBI")

    treeS = fileS.Get("ZZ4muTree").Get("candTree")
    treeB = fileB.Get("ZZ4muTree").Get("candTree")
    treeSBI = fileSBI.Get("ZZ4muTree").Get("candTree")

    treeS.Draw("p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM):ZZMass>>h2S","ZZMass<1600&&ZZMass>220","COLZ")
    treeB.Draw("p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM):ZZMass>>h2B","ZZMass<1600&&ZZMass>220","COLZ")
    treeSBI.Draw("p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM):ZZMass>>h2I","ZZMass<1600&&ZZMass>220","COLZ")

    h2S.SetName("T_2D_1")
    h2B.SetName("T_2D_2")
    h2I.SetName("T_2D_4")
    h2S.SetTitle("T_2D_1")
    h2B.SetTitle("T_2D_2")
    h2I.SetTitle("T_2D_4")

    #Normalize to desired rates
    h2S.Scale(rate_signal_ggzz_Shape/h2S.Integral())
    h2B.Scale(rate_signal_ggzz_Shape/h2B.Integral())
    #h2I.Scale(rate_signal_ggzz_Shape/h2I.Integral())

    for ix in range(h2I.GetNbinsX()):
        for iy in range(h2I.GetNbinsY()):
            h2I.SetBinContent(ix,iy,h2I.GetBinContent(ix,iy)-h2S.GetBinContent(ix,iy)-h2B.GetBinContent(ix,iy))


    #Save new templates
    fileout = TFile.Open("templates2e2mu.root","RECREATE")
    fileout.cd()
    h2S.Write()
    h2I.Write()
    h2B.Write()

# run the create_RM_cfg() as main()
if __name__ == "__main__":
    makeTemplates()
