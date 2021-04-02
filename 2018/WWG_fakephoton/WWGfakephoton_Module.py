#!/usr/bin/env python
# Analyzer for WWG Analysis based on nanoAOD tools

import os, sys
import math
import ROOT
from math import sin, cos, sqrt
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import countHistogramsProducer

class WWG_Producer(Module):
    def __init__(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree

        self.out.branch("event",  "i")
        self.out.branch("run",  "i")
        self.out.branch("lumi",  "i")
	self.out.branch("channel",  "I");
        self.out.branch("pass_selection1",  "B");
        self.out.branch("pass_selection2",  "B");

        self.out.branch("lep1_pid",  "I")
        self.out.branch("lep2_pid",  "I")
        self.out.branch("lep1pt",  "F")
        self.out.branch("lep2pt",  "F")
        self.out.branch("lep1eta",  "F")
        self.out.branch("lep2eta",  "F")
        self.out.branch("lep1phi",  "F")
        self.out.branch("lep2phi",  "F")
        self.out.branch("lepton1_isprompt", "I")
        self.out.branch("lepton2_isprompt", "I")
        self.out.branch("photonet",  "F")
        self.out.branch("photoneta",  "F")
        self.out.branch("photonphi",  "F")
        self.out.branch("photon_isprompt", "I")
        self.out.branch("photon_gen_matching", "I")
        self.out.branch("mll",  "F")
        self.out.branch("mllg",  "F")
        self.out.branch("ptll",  "F")
        self.out.branch("met",  "F")
        self.out.branch("metup",  "F")
        self.out.branch("puppimet","F")
        self.out.branch("puppimetphi","F")
        self.out.branch("rawmet","F")
        self.out.branch("rawmetphi","F")
        self.out.branch("metphi","F")
        self.out.branch("gen_weight","F")
        self.out.branch("npu",  "I");
        self.out.branch("ntruepu",  "F");
        self.out.branch("n_pos", "I")
        self.out.branch("n_minus", "I")
        self.out.branch("n_num", "I")
        self.out.branch("MET_pass","I")
        self.out.branch("npvs","I")
        self.out.branch("njets","I")
        self.out.branch("njets50","I")
        self.out.branch("njets40","I")
        self.out.branch("njets30","I")
        self.out.branch("njets20","I")
        self.out.branch("njets15","I")
        self.out.branch("HLT_Ele1","I")
        self.out.branch("HLT_Ele2","I")
        self.out.branch("HLT_Mu1","I")
        self.out.branch("HLT_Mu2","I")
        self.out.branch("HLT_emu1","I")
        self.out.branch("HLT_emu2","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.out.fillBranch("event",event.event)
        self.out.fillBranch("lumi",event.luminosityBlock)
        self.out.fillBranch("run",event.run)
#        print event.event,event.luminosityBlock,event.run
        if hasattr(event,'Generator_weight'):
            if event.Generator_weight > 0 :
                n_pos=1
                n_minus=0
            else:
                n_minus=1
                n_pos=0
            self.out.fillBranch("gen_weight",event.Generator_weight)
            self.out.fillBranch("n_pos",n_pos)
            self.out.fillBranch("n_minus",n_minus)
        else:    
            self.out.fillBranch("gen_weight",0)
            self.out.fillBranch("n_pos",0)
            self.out.fillBranch("n_minus",0)

        HLT_Ele1 = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
        HLT_Ele2 = event.HLT_Ele35_WPTight_Gsf

        HLT_Mu1 = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
        HLT_Mu2 = event.HLT_IsoMu24

        HLT_emu1 = event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
        HLT_emu2 = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ

        pass_selection1 = False
        pass_selection2 = False

        if not (HLT_Ele1 or HLT_Ele2 or HLT_Mu1 or HLT_Mu2 or HLT_emu1 or HLT_emu2):
           self.out.fillBranch("pass_selection",0)
           return False

        self.out.fillBranch("HLT_Ele1",HLT_Ele1)
        self.out.fillBranch("HLT_Ele2",HLT_Ele2)
        self.out.fillBranch("HLT_Mu1",HLT_Mu1)
        self.out.fillBranch("HLT_Mu2",HLT_Mu2)
        self.out.fillBranch("HLT_emu1",HLT_emu1)
        self.out.fillBranch("HLT_emu2",HLT_emu2)

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")
        jets = Collection(event, "Jet")
	if hasattr(event, 'nGenPart'):
           genparts = Collection(event, "GenPart")

        jet_select = [] 
        dileptonp4 = ROOT.TLorentzVector()
        selected_medium_or_control_photons = []
        selected_fake_template_photons = []
        electrons_select = []
        muons_select = [] 
        jets_select = []
        leptons_select=[]

        #selection on muons
        muon_pass =0
        for i in range(0,len(muons)):
            if muons[i].pt < 20:
                continue
            if abs(muons[i].eta) > 2.5:
                continue
            if muons[i].pfRelIso04_all > 0.25:
                continue   
            if muons[i].mediumId == True:
                muons_select.append(i)
                muon_pass += 1
                leptons_select.append(i)

        # selection on electrons
        electron_pass=0
        for i in range(0,len(electrons)):
            if electrons[i].pt < 20:
                continue
            if abs(electrons[i].eta + electrons[i].deltaEtaSC) > 2.5:
                continue
            if (abs(electrons[i].eta + electrons[i].deltaEtaSC) < 1.479 and abs(electrons[i].dz) < 0.1 and abs(electrons[i].dxy) < 0.05) or (abs(electrons[i].eta + electrons[i].deltaEtaSC) > 1.479 and abs(electrons[i].dz) < 0.2 and abs(electrons[i].dxy) < 0.1):
                if electrons[i].cutBased >= 3:
                    electrons_select.append(i)
                    electron_pass += 1
                    leptons_select.append(i)

#        print 'the number of leptons: ',len(electrons_select)+len(muons_select)
        if len(electrons_select)+len(muons_select) != 2:      #reject event if there are not exactly two leptons
	   return False

        # select medium photons
	photon_pass=0
        for i in range(0,len(photons)):
            if photons[i].pt < 20:
                continue
            if abs(photons[i].eta) > 2.5:
                continue
            if not (photons[i].isScEtaEE or photons[i].isScEtaEB):
                continue
            if photons[i].pixelSeed:
                continue

            #| pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho |
            mask1 = '0b10101010101010' # full medium ID
            mask2 = '0b00101010101010' # remove Isopho
            mask3 = '0b10001010101010' # remove IsoNeu
            mask4 = '0b10100010101010' # remove Isoch
            mask5 = '0b10101000101010' # remove sigma ieie

            bitmap = photons[i].vidNestedWPBitmap & mask1

            #after adding the photons that pass the full ID, add the photons that pass the inverted ID
            if not (bitmap == mask1):
                continue

            #this is redundant after the if statement above
            if not((bitmap == mask1) or (bitmap == mask2) or (bitmap == mask3) or (bitmap == mask4) or (bitmap == mask5)):
                continue

            pass_lepton_dr_cut = True
            for j in range(0,len(muons_select)):
                if deltaR(muons[muons_select[j]].eta,muons[muons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            for j in range(0,len(electrons_select)):
                if deltaR(electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False

            if not pass_lepton_dr_cut:
                continue

            photon_pass += 1
            selected_medium_or_control_photons.append(i)  #append the medium photons

        # select control photons
        for i in range(0,len(photons)):
            if photons[i].pt < 20:
                continue
            if abs(photons[i].eta) > 2.5:
                continue
            if not (photons[i].isScEtaEE or photons[i].isScEtaEB):
                continue
            if photons[i].pixelSeed:
                continue

            #| pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho |
            mask1 = '0b10101010101010' # full medium ID
            mask2 = '0b00101010101010' # remove Isopho
            mask3 = '0b10001010101010' # remove IsoNeu
            mask4 = '0b10100010101010' # remove Isoch
            mask5 = '0b10101000101010' # remove sigma ieie

            bitmap = photons[i].vidNestedWPBitmap & mask1

            #not pass the full ID
            if (bitmap == mask1):
                continue

            #pass ID that one of cut removed, which means this cut is inverted
            if not ((bitmap == mask1) or (bitmap == mask2) or (bitmap == mask3) or (bitmap == mask4) or (bitmap == mask5)):
                continue

            pass_lepton_dr_cut = True
            for j in range(0,len(muons_select)):
                if deltaR(muons[muons_select[j]].eta,muons[muons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            for j in range(0,len(electrons_select)):
                if deltaR(electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False

            if not pass_lepton_dr_cut:
                continue

            photon_pass += 1
            selected_medium_or_control_photons.append(i)  # append the control photons

        # select fake photons
	fakephoton_pass=0
        for i in range(0,len(photons)):
            if photons[i].pt < 20:
                continue
            if abs(photons[i].eta) > 2.5:
                continue
            if not (photons[i].isScEtaEE or photons[i].isScEtaEB):
                continue
            if photons[i].pixelSeed:
                continue

            #| pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho |
            mask1 = '0b10100000101010' # remove the Isoch and sigma ieie

            bitmap = photons[i].vidNestedWPBitmap & mask1

            #save photons pass the ID without sigma ieie and Isoch
            if not (bitmap == mask1):
                continue

            pass_lepton_dr_cut = True
            for j in range(0,len(muons_select)):
                if deltaR(muons[muons_select[j]].eta,muons[muons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            for j in range(0,len(electrons_select)):
                if deltaR(electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False

            if not pass_lepton_dr_cut:
                continue

            fakephoton_pass += 1
            selected_fake_template_photons.append(i) #for fake template from data

        if not (len(selected_medium_or_control_photons) >= 1 or len(selected_fake_template_photons) == 1):#only one jet fake photon 
            return False

        isprompt_mask = (1 << 0) #isPrompt
        isdirectprompttaudecayproduct_mask = (1 << 5) #isDirectPromptTauDecayProduct
        isdirecttaudecayproduct_mask = (1 << 4) #isDirectTauDecayProduct
        isprompttaudecayproduct = (1 << 3) #isPromptTauDecayProduct
        isfromhardprocess_mask = (1 << 8) #isPrompt

        channel = 0 
        # emu:     1
        # ee:      2
        # mumu:    3

        # emu
        mll = -10
        ptll = -10
        lepton1_isprompt=-10
        lepton2_isprompt=-10
        if len(muons_select)==1 and len(electrons_select)==1:  # emu channel 
            if deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
               return False 
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,muons[muons_select[0]].eta,muons[muons_select[0]].phi) < 0.5:
               return False
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
               return False
            if muons[muons_select[0]].charge * (electrons[electrons_select[0]].charge) >= 0:
               return False
#           print 'test emu channel',len(genparts)
            if hasattr(event, 'nGenPart'):
                print 'calculate the lepton flag in channel emu'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1
                       break 
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1 
                       break 
            channel = 1

            if len(selected_medium_or_control_photons) >= 1:
                if abs((electrons[electrons_select[0]].p4() + photons[selected_medium_or_control_photons[0]].p4()).M() - 91.2) >= 10:
                    pass_selection1 = True

            if len(selected_fake_template_photons) == 1:
                if abs((electrons[tight_electrons[0]].p4() + photons[selected_fake_template_photons[0]].p4()).M() - 91.2) >= 10:
                    pass_selection2 = True

            self.out.fillBranch("channel",channel)
	    self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
	    self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lep1_pid",13)
            self.out.fillBranch("lep2_pid",11)
            self.out.fillBranch("lep1pt",muons[muons_select[0]].pt)
            self.out.fillBranch("lep1eta",muons[muons_select[0]].eta)
            self.out.fillBranch("lep1phi",muons[muons_select[0]].phi)
            self.out.fillBranch("lep2pt",electrons[electrons_select[0]].pt)
            self.out.fillBranch("lep2eta",electrons[electrons_select[0]].eta)
            self.out.fillBranch("lep2phi",electrons[electrons_select[0]].phi)
            self.out.fillBranch("mll",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()).M())
            self.out.fillBranch("ptll",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()).Pt())
            self.out.fillBranch("photonet",photons[photons_select[0]].pt)
            self.out.fillBranch("photoneta",photons[photons_select[0]].eta)
            self.out.fillBranch("photonphi",photons[photons_select[0]].phi)
            self.out.fillBranch("mllg",(muons[muons_select[0]].p4() + electrons[electrons_select[0]].p4()+photons[photons_select[0]].p4()).M())
        # ee
        elif len(muons_select)==0 and len(electrons_select)==2:
            if deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi)<0.5:
               return False 
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi) < 0.5:
               return False
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi) < 0.5:
               return False
            if electrons[electrons_select[0]].charge * electrons[electrons_select[1]].charge >=0:
               return False 
#           print 'test',len(genparts)
	    if hasattr(event, 'nGenPart'):
                print 'calculate the lepton flag in channel ee'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electrons_select[0]].eta,electrons[electrons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1 
                       break 
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electrons_select[1]].eta,electrons[electrons_select[1]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1 
                       break 
            channel = 2
            if len(selected_medium_or_control_photons) >= 1:
                if (abs((electrons[electrons_select[0]].p4() + photons[selected_medium_or_control_photons[0]].p4()).M() - 91.2) >= 10) and (abs((electrons[electrons_select[1]].p4() + photons[selected_medium_or_control_photons[0]].p4()).M() - 91.2) >= 10):
                    pass_selection1 = True

            if len(selected_fake_template_photons) == 1:
                if (abs((electrons[tight_electrons[0]].p4() + photons[selected_fake_template_photons[0]].p4()).M() - 91.2) >= 10) and (abs((electrons[electrons_select[1]].p4() + photons[selected_medium_or_control_photons[0]].p4()).M() - 91.2) >= 10):
                    pass_selection2 = True

            self.out.fillBranch("channel",channel)
	    self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
	    self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lep1_pid",11)
            self.out.fillBranch("lep2_pid",11)
            self.out.fillBranch("lep1pt",electrons[electrons_select[0]].pt)
            self.out.fillBranch("lep1eta",electrons[electrons_select[0]].eta)
            self.out.fillBranch("lep1phi",electrons[electrons_select[0]].phi)
            self.out.fillBranch("lep2pt",electrons[electrons_select[1]].pt)
            self.out.fillBranch("lep2eta",electrons[electrons_select[1]].eta)
            self.out.fillBranch("lep2phi",electrons[electrons_select[1]].phi)
            self.out.fillBranch("mll",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()).M())
            self.out.fillBranch("ptll",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()).Pt())
            self.out.fillBranch("photonet",photons[photons_select[0]].pt)
            self.out.fillBranch("photoneta",photons[photons_select[0]].eta)
            self.out.fillBranch("photonphi",photons[photons_select[0]].phi)
            self.out.fillBranch("mllg",(electrons[electrons_select[0]].p4() + electrons[electrons_select[1]].p4()+photons[photons_select[0]].p4()).M())


        # mumu 
        elif len(electrons_select)==0 and len(muons_select)==2:
            if deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi)<0.5:
               return False 
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,muons[muons_select[0]].eta,muons[muons_select[0]].phi) < 0.5:
               return False
            if deltaR(photons[photons_select[0]].eta,photons[photons_select[0]].phi,muons[muons_select[1]].eta,muons[muons_select[1]].phi) < 0.5:
               return False
            if muons[muons_select[0]].charge * muons[muons_select[1]].charge >= 0:
               return False 
	    if hasattr(event, 'nGenPart'):
                print 'calculate the lepton flag in channel mumu'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muons_select[0]].eta,muons[muons_select[0]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton1_isprompt=1
                       break 
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muons_select[1]].eta,muons[muons_select[1]].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton2_isprompt=1
                       break 
            channel = 3
            pass_selection1 = len(selected_medium_or_control_photons) >= 1 # select medium and control photons
            pass_selection2 = len(selected_fake_template_photons) == 1     # select fake photons
            self.out.fillBranch("channel",channel)
	    self.out.fillBranch("lepton1_isprompt",lepton1_isprompt)
	    self.out.fillBranch("lepton2_isprompt",lepton2_isprompt)
            self.out.fillBranch("lep1_pid",13)
            self.out.fillBranch("lep2_pid",13)
            self.out.fillBranch("lep1pt",muons[muons_select[0]].pt)
            self.out.fillBranch("lep1eta",muons[muons_select[0]].eta)
            self.out.fillBranch("lep1phi",muons[muons_select[0]].phi)
            self.out.fillBranch("lep2pt",muons[muons_select[1]].pt)
            self.out.fillBranch("lep2eta",muons[muons_select[1]].eta)
            self.out.fillBranch("lep2phi",muons[muons_select[1]].phi)
            self.out.fillBranch("mll",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()).M())
            self.out.fillBranch("ptll",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()).Pt())
            self.out.fillBranch("photonet",photons[photons_select[0]].pt)
            self.out.fillBranch("photoneta",photons[photons_select[0]].eta)
            self.out.fillBranch("photonphi",photons[photons_select[0]].phi)
            self.out.fillBranch("mllg",(muons[muons_select[0]].p4() + muons[muons_select[1]].p4()+photons[photons_select[0]].p4()).M())

        else:
            return False

        if not pass_selection1 and not pass_selection2:
            return False

        njets = 0
        njets50 = 0
        njets40 = 0
        njets30 = 0
        njets20 = 0
        njets15 = 0
        njets_fake_template = 0
        njets_fake = 0
        pass_lepton_dr_cut = True
        for i in range(0,len(jets)):
            if jets[i].btagDeepB > 0.4184 and i<=6 :  # DeepCSVM
               return False
            if abs(jets[i].eta) > 4.7:
               continue
            if jets[i].pt<30:
               continue
	    if deltaR(jets[i].eta,jets[i].phi,photons[photons_select[0]].eta,photons[photons_select[0]].phi) < 0.5:
	       continue;
            for j in range(0,len(electrons_select)):
                if deltaR(jets[i].eta,jets[i].phi,electrons[electrons_select[j]].eta,electrons[electrons_select[j]].phi) < 0.5:
                   pass_lepton_dr_cut = False
            for j in range(0,len(muons_select)):
                if deltaR(jets[i].eta,jets[i].phi,muons[muons_select[j]].eta,muons[muons_select[j]].phi) < 0.5:
                   pass_lepton_dr_cut = False

            if  not pass_lepton_dr_cut == True:
	        continue
            if jets[i].jetId >> 1 & 1:
               jets_select.append(i)
               njets += 1
               if jets[i].pt > 50:
                   njets50+=1
               if jets[i].pt > 40:
                   njets40+=1
               if jets[i].pt > 30:
                   njets30+=1
               if jets[i].pt > 20:
                   njets20+=1
               if jets[i].pt > 15:
                   njets15+=1
            if len(selected_medium_or_control_photons) >= 1:

                if deltaR(photons[selected_medium_or_control_photons[0]].eta,photons[selected_medium_or_control_photons[0]].phi,jets[i].eta,jets[i].phi) > 0.5:
                    njets_fake+=1

            if len(selected_fake_template_photons) == 1:
                if deltaR(photons[selected_fake_template_photons[0]].eta,photons[selected_fake_template_photons[0]].phi,jets[i].eta,jets[i].phi) > 0.5:
                    njets_fake_template+=1
#        print len(jets),("njets",njets)
#        if njets >=2 :
#            self.out.fillBranch("pass_selection",0)
#            return True
        self.out.fillBranch("pass_selection1",pass_selection1) # select medium and control photons
        self.out.fillBranch("pass_selection2",pass_selection2) # select fake photons
####### from here 
        if pass_selection1:
           photon_gen_matching=-10
           photon_isprompt =-10
           if hasattr(photons[selected_medium_or_control_photons[0]],'genPartIdx') :
               print 'calculate the photon flag'
               if photons[selected_medium_or_control_photons[0]].genPartIdx >= 0 and genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].pdgId  == 22: 
                   if ((genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) and (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask):
                       photon_gen_matching = 6
                   elif ((genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):       
                       if (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].genPartIdxMother >= 0 and (abs(genparts[genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].genPartIdxMother].pdgId) == 11 or abs(genparts[genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].genPartIdxMother].pdgId) == 13 or abs(genparts[genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].genPartIdxMother].pdgId) == 15)):
                           photon_gen_matching = 4
                       else:    
                           photon_gen_matching = 5
                   else:
                       photon_gen_matching = 3
               elif photons[selected_medium_or_control_photons[0]].genPartIdx >= 0 and abs(genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].pdgId) == 11:     
                   if ((genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):  
                       photon_gen_matching = 1
                   else:
                       photon_gen_matching = 2
                       
               else:
                   assert(photons[selected_medium_or_control_photons[0]].genPartFlav == 0)
                   photon_gen_matching = 0
           if hasattr(event, 'nGenPart') :
               for j, genpart in enumerate(genparts):
	           if photons[selected_medium_or_control_photons[0]].genPartIdx >=0 and genpart.pt > 5 and abs(genpart.pdgId) == 22 and ((genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) or (genparts[photons[selected_medium_or_control_photons[0]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) and deltaR(photons[selected_medium_or_control_photons[0]].eta,photons[selected_medium_or_control_photons[0]].phi,genpart.eta,genpart.phi) < 0.3:
                      photon_isprompt =1
                      break
           mask1 = '0b10101010101010' # full medium ID
           mask2 = '0b00101010101010' # remove Isopho
           mask3 = '0b10001010101010' # remove IsoNeu
           mask4 = '0b10100010101010' # remove Isoch
           mask5 = '0b10101000101010' # remove sigma ieie
        
	   bitmap = photons[selected_medium_or_control_photons[0]].vidNestedWPBitmap & mask1   
           if (bitmap == mask1):
               self.out.fillBranch("photon_selection",1) #all cuts applied
           elif (bitmap == mask2):
               self.out.fillBranch("photon_selection",2) # remove Isopho
           elif (bitmap == mask3):
               self.out.fillBranch("photon_selection",3) # remove IsoNeu
           elif (bitmap == mask4):
               self.out.fillBranch("photon_selection",4) # remove Isoch
           elif (bitmap == mask5):
               self.out.fillBranch("photon_selection",5) # remove sigma ieie
           #pass_selection1 && photon_selection==5 -> build ture template from MC and data template from data
           #pass_selection1 && ((photon_selection!=1 && photon_selection==2) || (!1 && ==3) || (!1 && ==4) || (!=1 && ==5)) -> build fake photon enriched sample from data
           else:
               assert(0)
           self.out.fillBranch("photonet",photons[selected_medium_or_control_photons[0]].pt)
           self.out.fillBranch("photoneta",photons[selected_medium_or_control_photons[0]].eta)
           self.out.fillBranch("photonphi",photons[selected_medium_or_control_photons[0]].phi)
           self.out.fillBranch("photon_gen_matching",photon_gen_matching)
           self.out.fillBranch("photon_isprompt",photon_isprompt)

        elif pass_selection2:
           photon_gen_matching=-10
           photon_isprompt =-10
           if hasattr(photons[selected_fake_template_photons[0]],'genPartIdx') :
               print 'calculate the photon flag'
               if photons[selected_fake_template_photons[0]].genPartIdx >= 0 and genparts[photons[selected_fake_template_photons[0]].genPartIdx].pdgId  == 22: 
                   if ((genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) and (genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask):
                       photon_gen_matching = 6
                   elif ((genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):       
                       if (genparts[photons[selected_fake_template_photons[0]].genPartIdx].genPartIdxMother >= 0 and (abs(genparts[genparts[photons[selected_fake_template_photons[0]].genPartIdx].genPartIdxMother].pdgId) == 11 or abs(genparts[genparts[photons[selected_fake_template_photons[0]].genPartIdx].genPartIdxMother].pdgId) == 13 or abs(genparts[genparts[photons[selected_fake_template_photons[0]].genPartIdx].genPartIdxMother].pdgId) == 15)):
                           photon_gen_matching = 4
                       else:    
                           photon_gen_matching = 5
                   else:
                       photon_gen_matching = 3
               elif photons[selected_fake_template_photons[0]].genPartIdx >= 0 and abs(genparts[photons[selected_fake_template_photons[0]].genPartIdx].pdgId) == 11:     
                   if ((genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):  
                       photon_gen_matching = 1
                   else:
                       photon_gen_matching = 2
                       
               else:
                   assert(photons[selected_fake_template_photons[0]].genPartFlav == 0)
                   photon_gen_matching = 0
           if hasattr(event, 'nGenPart') :
               for j, genpart in enumerate(genparts):
	           if photons[selected_fake_template_photons[0]].genPartIdx >=0 and genpart.pt > 5 and abs(genpart.pdgId) == 22 and ((genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) or (genparts[photons[selected_fake_template_photons[0]].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) and deltaR(photons[selected_fake_template_photons[0]].eta,photons[selected_fake_template_photons[0]].phi,genpart.eta,genpart.phi) < 0.3:
                      photon_isprompt =1
                      break
           self.out.fillBranch("photonet",photons[selected_fake_template_photons[0]].pt)
           self.out.fillBranch("photoneta",photons[selected_fake_template_photons[0]].eta)
           self.out.fillBranch("photonphi",photons[selected_fake_template_photons[0]].phi)
           self.out.fillBranch("photon_gen_matching",photon_gen_matching)
           self.out.fillBranch("photon_isprompt",photon_isprompt)
           #pass_selection2  -> build fake template from data

        if hasattr(event,'Pileup_nPU'):    
            self.out.fillBranch("npu",event.Pileup_nPU)
        else:
            self.out.fillBranch("npu",0)
    
        if hasattr(event,'Pileup_nTrueInt'):    
            self.out.fillBranch("ntruepu",event.Pileup_nTrueInt)
        else:
            self.out.fillBranch("ntruepu",0)

        print 'channel', channel,'mu_pass:',muon_pass,' ele_pass:',electron_pass,' photon_pass:',photon_pass,' is lepton1 real ',lepton1_isprompt,' is lepton2 real ',lepton2_isprompt,' is photon real ',photon_isprompt,' or ',photon_gen_matching

        self.out.fillBranch("njets50",njets50)
        self.out.fillBranch("njets40",njets40)
        self.out.fillBranch("njets30",njets30)
        self.out.fillBranch("njets20",njets20)
        self.out.fillBranch("njets15",njets15)
        self.out.fillBranch("npvs",event.PV_npvs)
        self.out.fillBranch("met",event.MET_pt)
        self.out.fillBranch("metup",sqrt(pow(event.MET_pt*cos(event.MET_phi) + event.MET_MetUnclustEnUpDeltaX,2) + pow(event.MET_pt*sin(event.MET_phi) + event.MET_MetUnclustEnUpDeltaY,2)))
        self.out.fillBranch("puppimet",event.PuppiMET_pt)
        self.out.fillBranch("puppimetphi",event.PuppiMET_phi)
        self.out.fillBranch("rawmet",event.RawMET_pt)
        self.out.fillBranch("rawmetphi",event.RawMET_phi)
        self.out.fillBranch("metphi",event.MET_phi)
        self.out.fillBranch("pass_selection",1)
        return True

WWG_Module = lambda: WWG_Producer()

