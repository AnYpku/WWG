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

_rootLeafType2rootBranchType = {
    'UChar_t': 'b',
    'Char_t': 'B',
    'UInt_t': 'i',
    'Int_t': 'I',
    'Float_t': 'F',
    'Double_t': 'D',
    'ULong64_t': 'l',
    'Long64_t': 'L',
    'Bool_t': 'O'
}

class WWG_Producer(Module):
    def __init__(self,isdata=False):
	self.isdata = isdata
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        # Find list of activated branches in input tree
        _brlist_in = inputTree.GetListOfBranches()
        branches_in = set(
            [_brlist_in.At(i) for i in range(_brlist_in.GetEntries())])
        branches_in = [
            x for x in branches_in if inputTree.GetBranchStatus(x.GetName())
        ]
        # Find list of activated branches in output tree
        _brlist_out = wrappedOutputTree._tree.GetListOfBranches()
        branches_out = set(
            [_brlist_out.At(i) for i in range(_brlist_out.GetEntries())])
        branches_out = [
            x for x in branches_out
            if wrappedOutputTree._tree.GetBranchStatus(x.GetName())]
        # Use both
        branches = branches_in + branches_out
        self.isdata = not bool(inputTree.GetBranch("nGenPart"))

        # Only keep branches with right collection name
        self.brlist_lep, self.branchType_lep = self.filterBranchNames(branches, 'Lepton')
        self.brlist_jet, self.branchType_jet = self.filterBranchNames(branches, 'Jet')
        self.brlist_photon, self.branchType_photon = self.filterBranchNames(branches, 'Photon')

        # Create output branches
        self.out = wrappedOutputTree

        # basic
        self.out.branch("event",  "F")
        self.out.branch("run",  "F")
        self.out.branch("lumi",  "F")
	self.out.branch("photon_selection",  "I");
	self.out.branch("channel",  "I");
        self.out.branch("n_loose_mu", "I")
        self.out.branch("n_loose_ele", "I")
        self.out.branch("n_leptons", "I")

        # lepton
        for br in self.brlist_lep:
            self.out.branch("%s_%s" % ('lepton', br),
                            _rootLeafType2rootBranchType[self.branchType_lep[br]],
                            lenVar="nlepton")
        self.out.branch("lepton_is_tight", "O", lenVar="nlepton")
        if not self.isdata:
            self.out.branch("lepton_is_real", "O", lenVar="nlepton")

        # photon
        for br in self.brlist_photon:
            self.out.branch("%s_%s" % ('photon', br),
                            _rootLeafType2rootBranchType[self.branchType_photon[br]],
                            lenVar="nphoton")

        # ak4 jet
        for br in self.brlist_jet:
            self.out.branch("%s_%s" % ('jet', br),
                            _rootLeafType2rootBranchType[self.branchType_jet[br]],
                            lenVar="njet")

        self.out.branch("n_photon", "I")
        self.out.branch("photon_isprompt", "I")
        self.out.branch("photon_gen_matching", "I")
        self.out.branch("gen_weight","F")
        self.out.branch("npu",  "I");
        self.out.branch("ntruepu",  "F");
        self.out.branch("npvs","I")
        self.out.branch("njets50","I")
        self.out.branch("njets40","I")
        self.out.branch("njets30","I")
        self.out.branch("njets20","I")
        self.out.branch("njets15","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	pass

    def is_clean_with_collection(self, pat1, pat2_col, cone_size):
        for pat in pat2_col:
            if not deltaR(pat1.eta, pat1.phi, pat.eta, pat.phi) > cone_size:
                return False
        return True

    def filterBranchNames(self, branches, collection):
        out = []
        branchType = {}
        for br in branches:
            name = br.GetName()
            if not name.startswith(collection + '_'):
                continue
            out.append(name.replace(collection + '_', ''))
            branchType[out[-1]] = br.FindLeaf(br.GetName()).GetTypeName()
        return out, branchType

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.out.fillBranch("event",event.event)
        self.out.fillBranch("lumi",event.luminosityBlock)
        self.out.fillBranch("run",event.run)
#       print event.event,event.luminosityBlock,event.run
        if hasattr(event,'Generator_weight'):
            self.out.fillBranch("gen_weight",event.Generator_weight)
        else:    
            self.out.fillBranch("gen_weight",0)

	leptons = Collection(event, "Lepton")
        photons = Collection(event, "Photon")
        jets = Collection(event, "Jet")
	if hasattr(event, 'nGenPart'):
           genparts = Collection(event, "GenPart")

        veto_lepton_idx = []
        loose_leptons = []
        loose_leptons_pt = []
        loose_leptons_idx = []
        lepton_is_tight = []
        muons = []
        good_photon = []
        photons_select = []

	loose_lepton_pass=0
        # Find muons in leptons
        for i, lep in enumerate(leptons):
            if abs(lep.pdgId) == 13:
                if lep.looseId and lep.corrected_pt > 10 and abs(lep.eta) < 2.4:
                    veto_lepton_idx.append(i)
                    if lep.mediumId and lep.miniPFRelIso_all < 0.4:
                        muons.append(lep)
                        loose_leptons.append(lep)
                        # assume passed module: muonScaleResProducer, which will add rochester correction
                        loose_leptons_pt.append(lep.corrected_pt)
                        loose_leptons_idx.append(i)
                        lepton_is_tight.append(lep.mvaTTH > -0.2)
        # Find electrons in leptons
        for i, lep in enumerate(leptons):
            if abs(lep.pdgId) == 11:
                if lep.pt > 10 and abs(lep.eta + lep.deltaEtaSC) < 2.5:
                    if (abs(lep.eta + lep.deltaEtaSC) < 1.479 and abs(lep.dz) < 0.1 and abs(lep.dxy) < 0.05) or (
                            abs(lep.eta + lep.deltaEtaSC) > 1.479 and abs(lep.dz) < 0.2 and abs(lep.dxy) < 0.1):
                        if lep.cutBased >= 1:
                            veto_lepton_idx.append(i)
                            if self.is_clean_with_collection(lep, muons, 0.4):
                                loose_leptons.append(lep)
                                loose_leptons_pt.append(lep.pt)
                                loose_leptons_idx.append(i)
                                lepton_is_tight.append(lep.mvaFall17V2Iso_WP80)  # and lep.tightCharge==2

        # check if veto leptons are more than loose (tight+fakeable) leptons, if yes: there are additional leptons, veto this event
        if len(veto_lepton_idx) > len(loose_leptons_idx): return False
        # allowed lepton numbers: 1 (W), 2(WW/Z) and 3(WZ)
        if len(loose_leptons) < 1 or len(loose_leptons) > 3: return False

        # sort the leptons: according the pt
        zipped = zip(loose_leptons, loose_leptons_pt, loose_leptons_idx, lepton_is_tight)
        resort_zipped = sorted(zipped, key=lambda x: x[1], reverse=True)
        loose_leptons, loose_leptons_pt, loose_leptons_idx, lepton_is_tight = zip(
            *resort_zipped)  # now they are tuples, cannot change
        # leading leptons cut
        # if loose_leptons_pt[0].pt < 20 or loose_leptons_pt[1].pt < 15: return False

	loose_electron_pass = 0
	loose_muon_pass = 0
	for i, lep in enumerate(loose_leptons):
            if abs(lep.pdgId) == 13:
               loose_muon_pass += 1
            elif abs(lep.pdgId) == 11:
               loose_electron_pass +=1

        # Select medium photons
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
            pass_lepton_dr_cut = True
            for j in range(0,len(loose_leptons)):
                if deltaR(loose_leptons[j].eta,loose_leptons[j].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            if not pass_lepton_dr_cut:
                continue

            #| pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho |
            mask1 = 0b10101010101010 # full medium ID
            mask2 = 0b00101010101010 # fail Isopho
            mask3 = 0b10001010101010 # fail IsoNeu
            mask4 = 0b10100010101010 # fail Isoch
            mask5 = 0b10101000101010 # fail sigma ieie

            bitmap = photons[i].vidNestedWPBitmap & mask1

            #the photon pass the full ID
            if not (bitmap == mask1):
                continue

            #this is redundant after the if statement above
            if not ((bitmap == mask1) or (bitmap == mask2) or (bitmap == mask3) or (bitmap == mask4) or (bitmap == mask5)):
                continue
            photons_select.append(i)
	    good_photon.append(photons[i])
            photon_pass += 1

        # Select control photons
        for i in range(0,len(photons)):
            if photons[i].pt < 20:
                continue
            if abs(photons[i].eta) > 2.5:
                continue
            if not (photons[i].isScEtaEE or photons[i].isScEtaEB):
                continue
            if photons[i].pixelSeed:
                continue
            pass_lepton_dr_cut = True

            for j in range(0,len(loose_leptons)):
                if deltaR(loose_leptons[j].eta,loose_leptons[j].phi,photons[i].eta,photons[i].phi) < 0.5:
                    pass_lepton_dr_cut = False
            if not pass_lepton_dr_cut:
                continue

            #| pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho |
            mask1 = 0b10101010101010 # full medium ID
            mask2 = 0b00101010101010 # fail Isopho
            mask3 = 0b10001010101010 # fail IsoNeu
            mask4 = 0b10100010101010 # fail Isoch
            mask5 = 0b10101000101010 # fail sigma ieie

            bitmap = photons[i].vidNestedWPBitmap & mask1

            #not pass the full ID
            if (bitmap == mask1):
                continue

            #fail one of varaible in the ID
            if not((bitmap == mask1) or (bitmap == mask2) or (bitmap == mask3) or (bitmap == mask4) or (bitmap == mask5)):
                continue
            photons_select.append(i)
	    good_photon.append(photons[i])
            photon_pass += 1

        self.out.fillBranch("n_photon",photon_pass) #the number of medium photons and control photons


        isprompt_mask = (1 << 0) #isPrompt used for lepton
        isdirectprompttaudecayproduct_mask = (1 << 5) #isDirectPromptTauDecayProduct used for photon
        isprompttaudecayproduct = (1 << 3) #isPromptTauDecayProduct used for lepton
        isfromhardprocess_mask = (1 << 8) #isPrompt  used for photon


        lepton_is_real = []
        channel=0

        for j, lep in enumerate(loose_leptons):
            is_real_flag = False
            if hasattr(event, 'nGenPart'):
#               print 'calculate the lepton flag in channel emu'
                for i in range(0,len(genparts)):
		   if genparts[i].pt > 5 and abs(genparts[i].pdgId) == abs(lep.pdgId) and ((genparts[i].statusFlags & isprompt_mask == isprompt_mask) or (genparts[i].statusFlags & isprompttaudecayproduct == isprompttaudecayproduct)) and deltaR(lep.eta,lep.phi,genparts[i].eta,genparts[i].phi) < 0.3:
	               is_real_flag = True
                       break 
	        lepton_is_real.append(is_real_flag)
            else:
		lepton_is_real = [False] * len(loose_leptons)

        if len(loose_leptons) == 2:
	    if loose_muon_pass ==1 and loose_electron_pass==1:
               channel = 1
	    if loose_muon_pass ==2 and loose_electron_pass==0:
               channel = 2
	    if loose_muon_pass ==2 and loose_electron_pass==0:
               channel = 3
        else:
	    channel=4
        self.out.fillBranch("channel",channel)
        self.out.fillBranch("n_leptons",len(loose_leptons))
        self.out.fillBranch("n_loose_ele", loose_electron_pass)
        self.out.fillBranch("n_loose_mu", loose_muon_pass)
        self.out.fillBranch("lepton_is_tight", lepton_is_tight)
	self.out.fillBranch("lepton_is_real", lepton_is_real)

        for i in self.brlist_lep:
            if not i.endswith('corrected_pt'):
                self.out.fillBranch('lepton_' + i, [iobj[i] for iobj in loose_leptons])
            else:
                self.out.fillBranch('lepton_' + i, loose_leptons_pt)

        photon_gen_matching=-10
        photon_isprompt =-10
        if photon_pass>0:
	   photon_index=photons_select[0] 
           if hasattr(photons[photon_index],'genPartIdx') :
#               print 'calculate the photon flag'
               if photons[photon_index].genPartIdx >= 0 and genparts[photons[photon_index].genPartIdx].pdgId  == 22: 
                   if ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) and (genparts[photons[photon_index].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask):
                       photon_gen_matching = 6
                   elif ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):       
                       if (genparts[photons[photon_index].genPartIdx].genPartIdxMother >= 0 and (abs(genparts[genparts[photons[photon_index].genPartIdx].genPartIdxMother].pdgId) == 11 or abs(genparts[genparts[photons[photon_index].genPartIdx].genPartIdxMother].pdgId) == 13 or abs(genparts[genparts[photons[photon_index].genPartIdx].genPartIdxMother].pdgId) == 15)):
                           photon_gen_matching = 4
                       else:    
                           photon_gen_matching = 5
                   else:
                       photon_gen_matching = 3
               elif photons[photon_index].genPartIdx >= 0 and abs(genparts[photons[photon_index].genPartIdx].pdgId) == 11:     
                   if ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)):  
                       photon_gen_matching = 1
                   else:
                       photon_gen_matching = 2
                       
               else:
                   assert(photons[photon_index].genPartFlav == 0)
                   photon_gen_matching = 0
           if hasattr(event, 'nGenPart') and len(photons_select)>0 :
               for j, genpart in enumerate(genparts):
	           if photons[photon_index].genPartIdx >=0 and genpart.pt > 5 and abs(genpart.pdgId) == 22 and ((genparts[photons[photon_index].genPartIdx].statusFlags & isprompt_mask == isprompt_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) or (genparts[photons[photon_index].genPartIdx].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) and deltaR(photons[photon_index].eta,photons[photon_index].phi,genpart.eta,genpart.phi) < 0.3:
                      photon_isprompt =1
                      break
           mask1 = 0b10101010101010 # full medium ID
           mask2 = 0b00101010101010 # fail Isopho
           mask3 = 0b10001010101010 # fail IsoNeu
           mask4 = 0b10100010101010 # fail Isoch
           mask5 = 0b10101000101010 # fail sigma ieie
        
	   bitmap = photons[photon_index].vidNestedWPBitmap & mask1   
           if (bitmap == mask1):
               self.out.fillBranch("photon_selection",1) #all cuts applied
           elif (bitmap == mask2):
               self.out.fillBranch("photon_selection",2) # fail Isopho
           elif (bitmap == mask3):
               self.out.fillBranch("photon_selection",3) # fail IsoNeu
           elif (bitmap == mask4):
               self.out.fillBranch("photon_selection",4) # fail Isoch
           elif (bitmap == mask5):
               self.out.fillBranch("photon_selection",5) # fail sigma ieie
	   else:
               assert(0)
           #(photon_selection==2 || photon_selection==3 || photon_selection==4 || photon_selection ==5 )->build fake photon enriched sample 
	   for i in self.brlist_photon:
	       self.out.fillBranch('photon_' + i, [iobj[i] for iobj in good_photon])
        else:
           self.out.fillBranch("photon_selection",0) #if there is no photons selected
	   for i in self.brlist_photon:
	       self.out.fillBranch('photon_' + i,[ -10 for iobj in good_photon])
           
        self.out.fillBranch("photon_gen_matching",photon_gen_matching)
        self.out.fillBranch("photon_isprompt",photon_isprompt)

        pass_dr_cut = True
        njets50 = 0
        njets40 = 0
        njets30 = 0
        njets20 = 0
        njets15 = 0

        good_jets = []
        good_jets_idx = []
        for i in range(0,len(jets)):
            if abs(jets[i].eta) > 4.7:
               continue

            if photon_pass>0:
	       pass_dr_cut = deltaR(jets[i].eta,jets[i].phi,photons[photon_index].eta,photons[photon_index].phi) > 0.4
            for j in range(0,len(loose_leptons)):
                if deltaR(loose_leptons[j].eta,loose_leptons[j].phi,jets[i].eta,jets[i].phi) < 0.4:
                   pass_dr_cut = False

            if pass_dr_cut == False:
               continue

            if jets[i].jetId >> 1 & 1:
               good_jets.append(jets[i])
               good_jets_idx.append(i)
               if jets[i].pt_nom > 50:
                   njets50 +=1
               if jets[i].pt_nom > 40:
                   njets40 +=1
               if jets[i].pt_nom > 30:
                   njets30 +=1
               if jets[i].pt_nom > 20:
                   njets20 +=1
               if jets[i].pt_nom > 15:
                   njets15 +=1

        if hasattr(event,'Pileup_nPU'):    
            self.out.fillBranch("npu",event.Pileup_nPU)
        else:
            self.out.fillBranch("npu",0)
    
        if hasattr(event,'Pileup_nTrueInt'):    
            self.out.fillBranch("ntruepu",event.Pileup_nTrueInt)
        else:
            self.out.fillBranch("ntruepu",0)

        self.out.fillBranch("njets50",njets50)
        self.out.fillBranch("njets40",njets40)
        self.out.fillBranch("njets30",njets30)
        self.out.fillBranch("njets20",njets20)
        self.out.fillBranch("njets15",njets15)
        self.out.fillBranch("npvs",event.PV_npvs)
	if len(good_jets)>0:
           for i in self.brlist_jet:
                self.out.fillBranch('jet_' + i, [iobj[i] for iobj in good_jets])
	else:
           for i in self.brlist_jet:
                self.out.fillBranch('jet_' + i, [ -10 for iobj in good_jets])

        print 'channel',channel,', mu_pass:',loose_muon_pass,', ele_pass:',loose_electron_pass,', photon_pass:',photon_pass,', is photon real ',photon_isprompt,' or ',photon_gen_matching
#        print '------\n'

        return True

WWG_Module = lambda: WWG_Producer()

