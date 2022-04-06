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

        self.out.branch("event",  "F")
        self.out.branch("run",  "F")
        self.out.branch("lumi",  "F")

        self.out.branch("lepton_pid",  "I")
        self.out.branch("lepton_pt",  "F")
        self.out.branch("lepton_eta",  "F")
        self.out.branch("lepton_phi",  "F")
        self.out.branch("lepton_mvaTTH",  "F")
        self.out.branch("lepton_miniISO",  "F")
        self.out.branch("is_lepton_tight", "I")
	self.out.branch("lepton_isprompt", "I")
        self.out.branch("n_bjets","I")
        self.out.branch("mt",  "F")
        self.out.branch("puppimt",  "F")
        self.out.branch("met",  "F")
        self.out.branch("puppimet","F")
        self.out.branch("gen_weight","F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        self.out.fillBranch("event",event.event)
        self.out.fillBranch("lumi",event.luminosityBlock)
        self.out.fillBranch("run",event.run)
#        print event.event,event.luminosityBlock,event.run
        if hasattr(event,'Generator_weight'):
            self.out.fillBranch("gen_weight",event.Generator_weight)
        else:    
            self.out.fillBranch("gen_weight",0)

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        photons = Collection(event, "Photon")

	isprompt_mask = (1 << 0) #isPrompt
	isdirectprompttaudecayproduct_mask = (1 << 5) #isDirectPromptTauDecayProduct
	isdirecttaudecayproduct_mask = (1 << 4) #isDirectTauDecayProduct
	isprompttaudecayproduct = (1 << 3) #isPromptTauDecayProduct
	isfromhardprocess_mask = (1 << 8) #isPrompt
	   
        jets = Collection(event, "Jet")
	if hasattr(event, 'nGenPart'):
           genparts = Collection(event, "GenPart")

        jet_select = [] 
        dileptonp4 = ROOT.TLorentzVector()
        photons_select = []
        electrons_select = []
        muons_select = [] 
        jets_select = []
        leptons_select=[]
        loose_but_not_tight_muons = []
        loose_but_not_tight_electrons = []

        #selection on muons
        muon_pass =0
        for i in range(0,len(muons)):
            if muons[i].pt < 10:
                continue
            if abs(muons[i].eta) > 2.5:
                continue
            if muons[i].mediumId == True and muons[i].mvaTTH > -0.2 and muons[i].miniPFRelIso_all < 0.4:
                muons_select.append(i)
                muon_pass += 1
                leptons_select.append(i)
            elif muons[i].mediumId == True and muons[i].mvaTTH < -0.2 and muons[i].miniPFRelIso_all < 0.4:
                 loose_but_not_tight_muons.append(i)


        # selection on electrons
        electron_pass=0
        loose_electron_pass=0
        for i in range(0,len(electrons)):
            if electrons[i].pt < 10:
                continue
            if abs(electrons[i].eta + electrons[i].deltaEtaSC) > 2.5:
                continue
            if (abs(electrons[i].eta + electrons[i].deltaEtaSC) < 1.479 and abs(electrons[i].dz) < 0.1 and abs(electrons[i].dxy) < 0.05) or (abs(electrons[i].eta + electrons[i].deltaEtaSC) > 1.479 and abs(electrons[i].dz) < 0.2 and abs(electrons[i].dxy) < 0.1):
                if electrons[i].mvaFall17V2Iso_WP80==True:
                    electrons_select.append(i)
                    electron_pass += 1
                    leptons_select.append(i)
                elif electrons[i].cutBased >= 1:
                    loose_but_not_tight_electrons.append(i)

	lepton_isprompt=-10
        if len(muons_select) + len(loose_but_not_tight_muons) == 1 and len(electrons_select) + len(loose_but_not_tight_electrons) == 0:
	    if len(muons_select) == 1:
                muon_index = muons_select[0]
                self.out.fillBranch("is_lepton_tight",1)
            if len(loose_but_not_tight_muons) == 1:
                muon_index = loose_but_not_tight_muons[0]
                self.out.fillBranch("is_lepton_tight",0)
            if hasattr(event, 'nGenPart'):
	       for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 13 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(muons[muon_index].eta,muons[muon_index].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton_isprompt=1
                       break

            njets=0
            n_bjets=0
            pass_lepton_dr_cut = True
            for i in range(0,len(jets)):
                if jets[i].btagDeepB > 0.4168 and i<=6 :  # DeepCSVM
                   n_bjets += 1
                if abs(jets[i].eta) > 4.7:
                   continue
                if jets[i].pt<20:
                   continue
		if deltaR(jets[i].eta,jets[i].phi,muons[muon_index].eta,muons[muon_index].phi) < 0.3:
                       pass_lepton_dr_cut = False

                if  not pass_lepton_dr_cut == True:
	            continue
                if jets[i].jetId >> 1 & 1:
                   jets_select.append(i)
                   njets += 1
#           print 'fake muon, the number of jets ',njets
            if njets <1 :
               return False
            self.out.fillBranch("n_bjets",n_bjets)
            self.out.fillBranch("mt",sqrt(2*muons[muon_index].pt*event.MET_pt*(1 - cos(event.MET_phi - muons[muon_index].phi))))
            self.out.fillBranch("puppimt",sqrt(2*muons[muon_index].pt*event.PuppiMET_pt*(1 - cos(event.PuppiMET_phi - muons[muon_index].phi))))
            self.out.fillBranch("lepton_pt",muons[muon_index].pt)
            self.out.fillBranch("lepton_eta",muons[muon_index].eta)
            self.out.fillBranch("lepton_phi",muons[muon_index].phi)
            self.out.fillBranch("lepton_pid",muons[muon_index].pdgId)
            self.out.fillBranch("lepton_isprompt",lepton_isprompt)
            self.out.fillBranch("lepton_mvaTTH",muons[muon_index].mvaTTH)
            self.out.fillBranch("lepton_miniISO",muons[muon_index].miniPFRelIso_all)

	elif (len(electrons_select) + len(loose_but_not_tight_electrons) == 1) and (len(muons_select) + len(loose_but_not_tight_muons) == 0):

            if len(electrons_select) == 1:
                electron_index = electrons_select[0]
                self.out.fillBranch("is_lepton_tight",1)

            if len(loose_but_not_tight_electrons) == 1:
                electron_index = loose_but_not_tight_electrons[0]
                self.out.fillBranch("is_lepton_tight",0)

            if hasattr(event, 'nGenPart'):
	       for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == 11 and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(electrons[electron_index].eta,electrons[electron_index].phi,genparts[i].eta,genparts[i].phi) < 0.3:
                       lepton_isprompt=1
                       break

            njets=0
            n_bjets=0
            pass_lepton_dr_cut = True
            for i in range(0,len(jets)):
                if jets[i].btagDeepB > 0.4184 and i<=6 :  # DeepCSVM
                   n_bjets +=1       
                if abs(jets[i].eta) > 4.7:
                   continue
                if jets[i].pt<20:
                   continue
		if deltaR(jets[i].eta,jets[i].phi,electrons[electron_index].eta,electrons[electron_index].phi) < 0.3:
                       pass_lepton_dr_cut = False

                if  not pass_lepton_dr_cut == True:
	            continue
                if jets[i].jetId >> 1 & 1:
                   jets_select.append(i)
                   njets += 1
#           print 'fake ele, the number of jets ',njets
            if njets <1 :
               return False
            self.out.fillBranch("n_bjets",n_bjets)
            self.out.fillBranch("mt",sqrt(2*electrons[electron_index].pt*event.MET_pt*(1 - cos(event.MET_phi - electrons[electron_index].phi))))
            self.out.fillBranch("puppimt",sqrt(2*electrons[electron_index].pt*event.PuppiMET_pt*(1 - cos(event.PuppiMET_phi - electrons[electron_index].phi))))

            self.out.fillBranch("lepton_pt",electrons[electron_index].pt)
            self.out.fillBranch("lepton_eta",electrons[electron_index].eta)
            self.out.fillBranch("lepton_phi",electrons[electron_index].phi)
            self.out.fillBranch("lepton_pid",electrons[electron_index].pdgId)
            self.out.fillBranch("lepton_isprompt",lepton_isprompt)
            self.out.fillBranch("lepton_mvaTTH",electrons[electron_index].mvaTTH)
            self.out.fillBranch("lepton_miniISO",electrons[electron_index].miniPFRelIso_all)
	else:
	    return False

        self.out.fillBranch("met",event.MET_pt)
        self.out.fillBranch("puppimet",event.PuppiMET_pt)
#        print 'lepton is prompt', lepton_isprompt,' met',event.MET_pt,' the number of jets ',njets,'-> this event is saved'
        return True
WWGfakelepton_Module = lambda: WWG_Producer()
