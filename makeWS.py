#! /usr/bin/env python
import sys, pwd, commands
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
#from systematicsClass import *
#from inputReader import *

## ------------------------------------
##  card and workspace class
## ------------------------------------

## ------------------------------------
##  ISSUES:
##
##  qqZZ templates should be normalized xbin by xbin, I'm using roberto's
##  Do I still need signal rate normalizations, now that I use Ulash templates? Removed, but will need to be put back if using shape region < template region
##
## ------------------------------------


class width_datacardClass:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True
        
    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
            
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theLowSide, theOutputDir):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = 125.6   ## FIXED
        self.lumi =19.79
        self.inputlumi = 19.79
        self.sqrts = 8
        self.channel = 3
        self.outputDir = theOutputDir

        self.templRange =220
        
        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        ## ---------------- VARIABLES FOR LATER --------------- ##
        
        self.low_M = theLowSide
        self.high_M = 1600
        
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
            
        if (self.channel == self.ID_2e2mu): self.appendNameAlt = '2mu2e'
        else: self.appendNameAlt = self.appendName
        ## -------------------------- SIGNAL SHAPE VARIABLES ---------------------- ##
    
        bins = (self.high_M-self.templRange)/20 
        bins2 = (self.high_M-self.low_M)/20

        CMS_zz4l_widthMass_name = "CMS_zz4l_widthMass"
            
        CMS_zz4l_widthMass = ROOT.RooRealVar(CMS_zz4l_widthMass_name,CMS_zz4l_widthMass_name,self.low_M,self.high_M)
        CMS_zz4l_widthMass.setBins(bins2)

        x_name = "CMS_zz4l_GGsm"

        x = ROOT.RooRealVar(x_name,x_name,0.0001,100)
        x.setVal(1)
        x.setBins(100)

        mu_name = "CMS_zz4l_mu"

        mu = ROOT.RooRealVar(mu_name,mu_name,0.93,0.001,10)
        mu.setBins(100)

        mu_name = "CMS_widthH_kbkg"

        kbkg = ROOT.RooRealVar(mu_name,mu_name,0.1,10)
        kbkg.setVal(1.0)
        kbkg.setBins(100)

        D2name = "CMS_zz4l_widthKD"
        CMS_zz4l_widthKD = ROOT.RooRealVar(D2name,D2name,0.,1.)
        CMS_zz4l_widthKD.setBins(20)

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        print '2D signal shapes for Width'
        
        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##

        
        #templates
        templateSigName = "/afs/cern.ch/work/u/usarica/public/CombineTemplates/02_04_2014/{0:.0f}TeV/{1}/{2}/HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_Nominal.root".format(self.sqrts,self.appendNameAlt,self.templRange)
        templateBkgName = templateSigName
        newtemplateSigName = "templates2e2mu.root"
        templateSigNameUp = "/afs/cern.ch/work/u/usarica/public/CombineTemplates/02_04_2014/{0:.0f}TeV/{1}/{2}/HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysUp_ggPDF.root".format(self.sqrts,self.appendNameAlt,self.templRange)
        templateSigNameDown = "/afs/cern.ch/work/u/usarica/public/CombineTemplates/02_04_2014/{0:.0f}TeV/{1}/{2}/HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysDown_ggPDF.root".format(self.sqrts,self.appendNameAlt,self.templRange)

        sigTempFileU = ROOT.TFile(templateSigName)
        sigTempFileG = ROOT.TFile(templateSigName)
        bkgTempFileU = ROOT.TFile(templateBkgName)
        sigTempFileUp = ROOT.TFile(templateSigNameUp)
        sigTempFileDown = ROOT.TFile(templateSigNameDown)

        Sig_T_1 = sigTempFileG.Get("T_2D_2").Clone("mZZ_bkg") 
        Sig_T_2 = sigTempFileU.Get("T_2D_1").Clone("mZZ_sig")
        Sig_T_4 = sigTempFileU.Get("T_2D_4").Clone("mZZ_inter")
        Bkg_T = bkgTempFileU.Get("T_2D_qqZZ_UnConditional").Clone("mZZ_bkg")

        Sig_T_1_Up = sigTempFileUp.Get("T_2D_2").Clone("T_2D_2_Up")
        Sig_T_2_Up = sigTempFileUp.Get("T_2D_1").Clone("T_2D_1_Up")
        Sig_T_4_Up = sigTempFileUp.Get("T_2D_4").Clone("T_2D_4_Up")
        #Bkg_T_Up = sigTempFileUp.Get("T_2D_qqZZ").Clone("T_2D_qqZZ_Up")
        Sig_T_1_Down = sigTempFileDown.Get("T_2D_2").Clone("T_2D_2_Down")
        Sig_T_2_Down = sigTempFileDown.Get("T_2D_1").Clone("T_2D_1_Down")
        Sig_T_4_Down = sigTempFileDown.Get("T_2D_4").Clone("T_2D_4_Down")
        #Bkg_T_Down = sigTempFileDown.Get("T_2D_qqZZ").Clone("T_2D_qqZZ_Down")

        #rates
        totalRateDown = Sig_T_1_Down.Integral("width")+Sig_T_2_Down.Integral("width")+Sig_T_4_Down.Integral("width")
        totalRateUp = Sig_T_1_Up.Integral("width")+Sig_T_2_Up.Integral("width")+Sig_T_4_Up.Integral("width")
        totalRate_ggzz = Sig_T_1.Integral("width")+Sig_T_2.Integral("width")+Sig_T_4.Integral("width")
        totalRateDown = totalRateDown #*2.3
        totalRateUp  = totalRateUp #*2.3
        rate_signal_ggzz_Shape = Sig_T_2.Integral("width")*self.lumi #*2.3
        rate_bkg_ggzz_Shape = Sig_T_1.Integral("width")*self.lumi #*2.3
        rate_interf_ggzz_Shape = Sig_T_4.Integral("width")*self.lumi #*2.3

        ## rates per lumi for scaling
        bkgRate_qqzz = 76.82#theInputs['qqZZ_rate']/theInputs['qqZZ_lumi'] #*1.8

        totalRate_ggzz_Shape = totalRate_ggzz*self.lumi
        bkgRate_qqzz_Shape = bkgRate_qqzz*self.lumi
        #bkgRate_zjets_Shape = bkgRate_zjets*self.lumi

        
        if Sig_T_4.Integral()<0 : #negative interference, turn it positive, the sign will be taken into account later when building the pdf
            print "negative interference"
            for ix in range (1,Sig_T_4.GetXaxis().GetNbins()+1):
                for iy in range (1,Sig_T_4.GetYaxis().GetNbins()+1):
                    Sig_T_4.SetBinContent(ix,iy,-1.0*Sig_T_4.GetBinContent(ix,iy))
                    Sig_T_4_Down.SetBinContent(ix,iy,-1.0*Sig_T_4_Down.GetBinContent(ix,iy))
                    Sig_T_4_Up.SetBinContent(ix,iy,-1.0*Sig_T_4_Up.GetBinContent(ix,iy))
            rate_interf_ggzz_Shape = Sig_T_4.Integral("width")*self.lumi

        #normalization on background and protection against negative fluctuations
        for ix in range(1,Bkg_T.GetXaxis().GetNbins()+1):
            yNorm = Bkg_T.Integral(ix,ix,1,Bkg_T.GetYaxis().GetNbins())
            #print yNorm
            if yNorm == 0: yNorm = 1.0
            for iy in range(1,Bkg_T.GetYaxis().GetNbins()+1):
                Bkg_T.SetBinContent(ix,iy,Bkg_T.GetBinContent(ix,iy)/yNorm)
                if Bkg_T.GetBinContent(ix,iy) == 0: Bkg_T.SetBinContent(ix,iy,0.000001)
                binI = Sig_T_4.GetBinContent(ix,iy)
                if binI > 0 : #check signs, should be < 0 for the template but I changed the sign above (secondo me >0)
                    binS = Sig_T_2.GetBinContent(ix,iy)
                    binB = Sig_T_1.GetBinContent(ix,iy)
                    if binI*binI >= 4*binS*binB:
                        Sig_T_4.SetBinContent(ix,iy,sqrt(abs(4*binS*binB))-0.00001)#check signs (secondo me 4 -0.0)
                        
                binI = Sig_T_4_Up.GetBinContent(ix,iy)
                if binI > 0 : #check signs, should be < 0 for the template but I changed the sign above (secondo me >0)
                    binS = Sig_T_2_Up.GetBinContent(ix,iy)
                    binB = Sig_T_1_Up.GetBinContent(ix,iy)
                    if binI*binI >= 4*binS*binB:
                        Sig_T_4_Up.SetBinContent(ix,iy,sqrt(abs(4*binS*binB))-0.00001)#check signs (secondo me 4 -0.0)

                binI = Sig_T_4_Down.GetBinContent(ix,iy)
                if binI > 0 : #check signs, should be < 0 for the template but I changed the sign above (secondo me >0)
                    binS = Sig_T_2_Down.GetBinContent(ix,iy)
                    binB = Sig_T_1_Down.GetBinContent(ix,iy)
                    if binI*binI >= 4*binS*binB:
                        Sig_T_4_Down.SetBinContent(ix,iy,sqrt(abs(4*binS*binB))-0.00001)#check signs (secondo me 4 -0.0)

        #protection against empty bins
        for ix in range(1,Bkg_T.GetXaxis().GetNbins()+1):
            for iy in range(1,Bkg_T.GetYaxis().GetNbins()+1):
                if Bkg_T.GetBinContent(ix,iy) == 0 : Bkg_T.SetBinContent(ix,iy,0.000001)
                if Sig_T_1.GetBinContent(ix,iy) == 0 : Sig_T_1.SetBinContent(ix,iy,0.000001)
                if Sig_T_2.GetBinContent(ix,iy) == 0 : Sig_T_2.SetBinContent(ix,iy,0.000001)
                if Sig_T_4.GetBinContent(ix,iy) == 0 : Sig_T_4.SetBinContent(ix,iy,0.000001)
                if Sig_T_1_Up.GetBinContent(ix,iy) == 0 : Sig_T_1_Up.SetBinContent(ix,iy,0.000001)
                if Sig_T_2_Up.GetBinContent(ix,iy) == 0 : Sig_T_2_Up.SetBinContent(ix,iy,0.000001)
                if Sig_T_4_Up.GetBinContent(ix,iy) == 0 : Sig_T_4_Up.SetBinContent(ix,iy,0.000001)
                if Sig_T_1_Down.GetBinContent(ix,iy) == 0 : Sig_T_1_Down.SetBinContent(ix,iy,0.000001)
                if Sig_T_2_Down.GetBinContent(ix,iy) == 0 : Sig_T_2_Down.SetBinContent(ix,iy,0.000001)
                if Sig_T_4_Down.GetBinContent(ix,iy) == 0 : Sig_T_4_Down.SetBinContent(ix,iy,0.000001)                                


        dBinsX = Sig_T_1.GetXaxis().GetNbins()
        dBinsY = Sig_T_1.GetYaxis().GetNbins()
        
        CMS_zz4l_widthMass.setBins(dBinsX)

        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)

        ## -------------------------- gg2ZZ SHAPES ---------------------------------- ##
        
        sigRateName = "signal_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigRates = ROOT.RooRealVar(sigRateName,sigRateName,0.0,10000.0)
        bkgRateName = "bkg_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgRates = ROOT.RooRealVar(bkgRateName,bkgRateName,0.0,10000.0)
        interfRateName = "interf_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        interfRates = ROOT.RooRealVar(interfRateName,interfRateName,0.0,10000.0)
        
        sigRateNameNorm = "signalNorm_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigRatesNorm = ROOT.RooFormulaVar(sigRateNameNorm,"@0*@1/(@0*@1-sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))+@2)",ROOT.RooArgList(x,mu,kbkg))
        interfRateNameNorm = "interfNorm_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        interfRatesNorm = ROOT.RooFormulaVar(interfRateNameNorm,"-sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))/(@0*@1-sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))+@2)",ROOT.RooArgList(x,mu,kbkg))
        bkgRateNameNorm = "bkgNorm_ggZZrate_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgRatesNorm = ROOT.RooFormulaVar(bkgRateNameNorm,"@2/(@0*@1-sqrt(@0*@1)*sign(@2)*sqrt(abs(@2))+@2)",ROOT.RooArgList(x,mu,kbkg))

        TemplateName = "ggZZsignal_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZsignal_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZsignal_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_2)
        ggZZsignal_TemplatePdf = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZsignal_TempDataHist)

        TemplateName = "ggZZbkg_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZbkg_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZbkg_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_1)
        ggZZbkg_TemplatePdf = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZbkg_TempDataHist)
            
        TemplateName = "ggZZinterf_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZinterf_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZinterf_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_4)
        ggZZinterf_TemplatePdf = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZinterf_TempDataHist)
            
        ggZZpdfName = "ggZZ_RooWidth_Nominal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZpdf_Nominal = ROOT.RooRealSumPdf(ggZZpdfName,ggZZpdfName,ROOT.RooArgList(ggZZsignal_TemplatePdf,ggZZinterf_TemplatePdf,ggZZbkg_TemplatePdf),ROOT.RooArgList(sigRatesNorm,interfRatesNorm,bkgRatesNorm))

        ## -------------------------- SHAPE Systematic ---------------------------------- ##

        #Up Systematics pdf
        TemplateName = "ggZZsignal_TempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZsignal_TemplatePdf_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZsignal_TempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_2_Up)
        ggZZsignal_TemplatePdf_Up = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZsignal_TempDataHist_Up)
            
        TemplateName = "ggZZbkg_TempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZbkg_TemplatePdf_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        
        ggZZbkg_TempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_1_Up)
        ggZZbkg_TemplatePdf_Up = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZbkg_TempDataHist_Up)
            
        TemplateName = "ggZZinterf_TempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZinterf_TemplatePdf_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
    
        ggZZinterf_TempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_4_Up)
        ggZZinterf_TemplatePdf_Up = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZinterf_TempDataHist_Up)

        ggZZpdfName = "ggZZ_RooWidth_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZpdf_Up = ROOT.RooRealSumPdf(ggZZpdfName,ggZZpdfName,ROOT.RooArgList(ggZZsignal_TemplatePdf_Up,ggZZinterf_TemplatePdf_Up,ggZZbkg_TemplatePdf_Up),ROOT.RooArgList(sigRatesNorm,interfRatesNorm,bkgRatesNorm))

        #Down Systematics pdf        
        TemplateName = "ggZZsignal_TempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZsignal_TemplatePdf_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZsignal_TempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_2_Down)
        ggZZsignal_TemplatePdf_Down = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZsignal_TempDataHist_Down)
            
        TemplateName = "ggZZbkg_TempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZbkg_TemplatePdf_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZbkg_TempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_1_Down)
        ggZZbkg_TemplatePdf_Down = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZbkg_TempDataHist_Down)
        
        TemplateName = "ggZZinterf_TempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        PdfName = "ggZZinterf_TemplatePdf_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        ggZZinterf_TempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Sig_T_4_Down)
        ggZZinterf_TemplatePdf_Down = ROOT.RooHistFunc(PdfName,PdfName,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),ggZZinterf_TempDataHist_Down)

        ggZZpdfName = "ggZZ_RooWidth_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggZZpdf_Down = ROOT.RooRealSumPdf(ggZZpdfName,ggZZpdfName,ROOT.RooArgList(ggZZsignal_TemplatePdf_Down,ggZZinterf_TemplatePdf_Down,ggZZbkg_TemplatePdf_Down),ROOT.RooArgList(sigRatesNorm,interfRatesNorm,bkgRatesNorm))

        CMS_zz4l_APscale_syst = ROOT.RooRealVar("CMS_zz4l_pdf_QCDscale_gg_syst","CMS_zz4l_pdf_QCDscale_gg_syst",0.0,-1,1)
        morphVarListggZZ = ROOT.RooArgList()
        morphVarListggZZ.add(CMS_zz4l_APscale_syst)
        MorphList_ggZZ = ROOT.RooArgList()
        MorphList_ggZZ.add(ggZZpdf_Nominal)
        MorphList_ggZZ.add(ggZZpdf_Up)
        MorphList_ggZZ.add(ggZZpdf_Down)
        
        ggZZpdf = ROOT.VerticalInterpPdf("ggzz","ggzz",MorphList_ggZZ,morphVarListggZZ)
        #ggZZpdf = ggZZpdf_Nominal

        asympowname = "kappalow_ggZZ_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        kappalow = ROOT.RooRealVar(asympowname,asympowname,totalRateDown/totalRate_ggzz)#kappalow = ROOT.RooRealVar(asympowname,asympowname,rateSignal_Down+rateBkg_Down-rateInterf_Down)
        asympowname = "kappahigh_ggZZ_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        kappahigh = ROOT.RooRealVar(asympowname,asympowname,totalRateUp/totalRate_ggzz)#kappahigh = ROOT.RooRealVar(asympowname,asympowname,rateSignal_Up+rateBkg_Up-rateInterf_Up)        
        asympowname = "Asympow_ggZZ_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        thetaSyst_ggZZ = AsymPow(asympowname,asympowname,kappalow,kappahigh,CMS_zz4l_APscale_syst)


        ### -------------------------- OTHER BACKGROUND SHAPES ---------------------------------- ##
     
        TemplateName = "qqzz_TempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqzz_TempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_widthMass,CMS_zz4l_widthKD),Bkg_T)
        PdfName = "qqzz_TemplatePdf_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        qqzz_TemplatePdf = ROOT.RooHistPdf("bkg_qqzz","bkg_qqzz",ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),qqzz_TempDataHist)
        #bkg_qqzz = ROOT.RooProdPdf("bkg_qqzz","bkg_qqzz",ROOT.RooArgSet(bkg_qqzz_mass),ROOT.RooFit.Conditional(ROOT.RooArgSet(qqzz_TemplatePdf),ROOT.RooArgSet(CMS_zz4l_widthKD)))

        #bkg_qqzz = ROOT.RooHistFunc("bkg_qqzz","bkg_qqzz",ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD),bkg_qqzz_mass_TempDataHist_Up)
         ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
         ## ----------------------- SIGNAL AND BACKGROUND RATES ----------------------- ##

        sigRates.setVal(rate_signal_ggzz_Shape)
        sigRates.setConstant(True)
        bkgRates.setVal(rate_bkg_ggzz_Shape)
        bkgRates.setConstant(True)
        interfRates.setVal(rate_interf_ggzz_Shape)
        interfRates.setConstant(True)

        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",low_M," - ",high_M
            CMS_zz4l_widthMass.setRange("shapiro",100.,160.)
            bkgRate_qqzz_shapiro = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_widthMass), ROOT.RooFit.Range("shapiro") ).getVal()
            totalRate_ggzz_shapiro = sclFactorTotal_ggzz * ggZZpdf.createIntegral( ROOT.RooArgSet(CMS_zz4l_widthMass), ROOT.RooFit.Range("shapiro") ).getVal()
            bkgRate_zjets_shapiro = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_widthMass), ROOT.RooFit.Range("shapiro") ).getVal()
            lowmassyield = bkgRate_qqzz_shapiro + totalRate_ggzz_shapiro + bkgRate_zjets_shapiro
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.inputlumi)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        

        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_widthMass,CMS_zz4l_widthKD))
        data_obs_red = data_obs.reduce("CMS_zz4l_widthMass > {0}".format(self.low_M))
            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.low_M
        if ( math.fabs(math.floor(tmpMH)-self.low_M) > 0.001): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""

        if (endsInP5): name_Shape = "{0}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        else: name_Shape = "{0}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        
        if (endsInP5): name_ShapeWS = "{0}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        else: name_ShapeWS = "{0}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.low_M,self.appendName,self.sqrts)
        
        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)
                
                
        getattr(w,'import')(data_obs_red,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?

        #ggZZpdf_nominal.SetNameTitle("ggzz_nominal","ggzz_nominal")
        #getattr(w,'import')(ggZZpdf_nominal, ROOT.RooFit.RecycleConflictNodes())

        ggZZpdf.SetNameTitle("ggzz","ggzz")
        getattr(w,'import')(ggZZpdf, ROOT.RooFit.RecycleConflictNodes())

        #ggZZpdf_Down.SetNameTitle("ggzz_CMS_zz4l_scale_systDown","ggzz_CMS_zz4l_scale_systDown")
        #getattr(w,'import')(ggZZpdf_Down, ROOT.RooFit.RecycleConflictNodes())

        #ggZZpdf_Up.SetNameTitle("ggzz_CMS_zz4l_scale_systUp","ggzz_CMS_zz4l_scale_systUp")
        #getattr(w,'import')(ggZZpdf_Up, ROOT.RooFit.RecycleConflictNodes())


        #getattr(w,'import')(CMS_zz4l_syst, ROOT.RooFit.RecycleConflictNodes())

        #w.factory("EXPR::ggzz('ggzz_nominal*(one+ggzz_CMS_zz4l_scale_systUp*CMS_zz4l_syst)',one[1],ggzz_nominal,ggzz_CMS_zz4l_scale_systUp,CMS_zz4l_syst)")
        
        ggZZpdfNormName = "ggZZ_RooWidth_{0:.0f}_{1:.0f}_norm".format(self.channel,self.sqrts)
        #ggZZpdf_norm = ROOT.RooFormulaVar(ggZZpdfNormName,"@0*@3*@4-@1*sqrt(@3*@4)*sign(@5)*sqrt(abs(@5))+@2*@5",ROOT.RooArgList(sigRates,interfRates,bkgRates,x,mu,kbkg))
        ggZZpdf_norm = ROOT.RooFormulaVar(ggZZpdfNormName,"(@0*@3*@4-@1*sqrt(@3*@4)*sign(@5)*sqrt(abs(@5))+@2*@5)*@6",ROOT.RooArgList(sigRates,interfRates,bkgRates,x,mu,kbkg,thetaSyst_ggZZ))
        ggZZpdf_norm.SetNameTitle("ggzz_norm","ggzz_norm")
        getattr(w,'import')(ggZZpdf_norm, ROOT.RooFit.RecycleConflictNodes())

        qqzz_TemplatePdf.SetNameTitle("bkg_qqzz","bkg_qqzz")
        getattr(w,'import')(qqzz_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        #bkg_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
        #getattr(w,'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())

        #bkg_qqzz_norm.SetNameTitle("bkg_qqzz_norm","bkg_qqzz_norm")
        #getattr(w,'import')(bkg_qqzz_norm, ROOT.RooFit.RecycleConflictNodes())

        ##ggZZsignal_TemplatePdf.SetNameTitle("ggsignalzz","ggsignalzz")
        ##ggZZbkg_TemplatePdf.SetNameTitle("ggbkgzz","ggbkgzz")
        ##ggZZinterf_TemplatePdf.SetNameTitle("gginterfzz","gginterfzz")
        ##getattr(w,'import')(ggZZsignal_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        ##getattr(w,'import')(ggZZbkg_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
        ##getattr(w,'import')(ggZZinterf_TemplatePdf, ROOT.RooFit.RecycleConflictNodes())
  
        w.writeToFile(name_ShapeWS)

    def WriteDatacard(self,file,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)
        
        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS w:$PROCESS_$SYSTEMATIC\n".format(nameWS))
        file.write("------------\n")
        

        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("## signal,bkg,interf,tot rates [{0:.4f}, {1:.4f}, -{2:.4f}, {3:.4f}] \n".format(theRates["ggZZ_signal"],theRates["ggZZ_bkg"],theRates["ggZZ_interf"],theRates["ggZZ_tot"]))
        file.write("bin ")        

        channelList=['ggZZ','qqZZ','zjets'] 

        channelName=['ggzz','vbf_offshell','bkg_qqzz','bkg_zjets'] 
         
        for chan in channelList:
            if theInputs[chan]:
                file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0

        for chan in channelList:
            #print 'checking if ',chan,' is in the list of to-do'
            #print "{0} ".format(channelName[i])
            if theInputs[chan]:
                file.write("{0} ".format(channelName[i]))
                #print 'writing in card index=',i,'  chan=',chan
                #print "{0} ".format(channelName[i])
            i+=1

        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        for chan in channelList:
            if theInputs[chan]:
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


# run the create_RM_cfg() as main()
if __name__ == "__main__":
    dc =width_datacardClass()
    dc.loadIncludes()
    cmd = 'mkdir -p provaX/'
    status, output = commands.getstatusoutput(cmd)    
    dc.makeCardsWorkspaces(220,"provaX")
