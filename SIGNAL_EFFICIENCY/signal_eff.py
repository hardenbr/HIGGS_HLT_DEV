import os
import itertools
import array
import ROOT as rt
import rootlogon
rootlogon.style()

pileup20x25 = rt.TFile("/afs/cern.ch/work/h/hardenbr/2014/HIGGS_DIPHOTON_HLT/PROD_20x25/res/tree.root")
pileup40x25 = rt.TFile("/afs/cern.ch/work/h/hardenbr/2014/HIGGS_DIPHOTON_HLT/PROD_40x25/res/tree.root")
pileup40x50 = rt.TFile("/afs/cern.ch/work/h/hardenbr/2014/HIGGS_DIPHOTON_HLT/PROD_40x50/res/tree.root")
pileup80x50 = rt.TFile("/afs/cern.ch/work/h/hardenbr/2014/HIGGS_DIPHOTON_HLT/PROD_80x50/res/tree.root")
#data = rt.TFile("hlt_tree_run208390.root")

    
#declare the tree vectors
sieie_seed = rt.vector('float')()
ecal_iso_seed = rt.vector('float')()
h_for_hoe_seed = rt.vector('float')()
hoe_seed  = rt.vector('float')()
hcal_iso_seed  = rt.vector('float')()
r9_seed  = rt.vector('float')()

#no seed requirement
sieie  = rt.vector('float')()
ecal_iso  = rt.vector('float')()
h_for_hoe = rt.vector('float')()
hcal_iso = rt.vector('float')()
track_iso = rt.vector('float')()
r9   = rt.vector('float')()

#kinematics
pt = rt.vector('float')()
eta = rt.vector('float')()
phi = rt.vector('float')()
energy = rt.vector('float')()

#kinematics
pt_seed = rt.vector('float')()
eta_seed = rt.vector('float')()
phi_seed = rt.vector('float')()
energy_seed = rt.vector('float')()

#get the selected pair of photons passing id
def getSelectionPair(cut):
    good_id = []

    #we must have a seed
    if not hasSeed(cut): return (-99,-99,-99,-99)

    #make a list of photons which satisfy the photon id
    for ii in range(len(pt)):
        if photonPassID(ii,cut): good_id.append(ii)

    #take a all possible combinations
    pairs_id = list(itertools.product(good_id,good_id))

    #make a list of pairs of photons that pass id
    good_pairs = []

    for pair in pairs_id:
        x = pair[0]
        y = pair[1]

        #cant be two of the same photon
        if x == y: continue

        #only look at half the events and
        #require the first index to have higher pt
        if pt[x] < pt[y]: continue

        #check the kinematic requirement
        pt_cut = pt[x] > 26 and pt[y] > 18
        if not pt_cut: continue

        #check the mass
        v1 = rt.TLorentzVector()
        v2 = rt.TLorentzVector()

        #set the kinematic
        v1.SetPtEtaPhiM(pt[x],eta[x],phi[x],0)        
        v2.SetPtEtaPhiM(pt[y],eta[y],phi[y],0)

        mass = (v1 + v2).Mag()
        
        #check the mass cut
        if mass > 70.0: good_pairs.append(pair)

    highest_sum_pt = -99
    best_pair = (-99,-99)
    pt_high = -99
    pt_low = -99
    
    for pair in good_pairs:
        x = pair[0]
        y = pair[1]

        sum = pt[x] + pt[y]
        #find the highest sum pt pair
        if sum > highest_sum_pt:
            best_pair = pair
            highest_sum_pt = sum
            pt_high = pt[pair[0]]
            pt_low = pt[pair[1]]

    return (pt_high,pt_low) + best_pair

def hasSeed(cut):
    found = False

    #loop over the number of seed photons
    for index in range(len(pt_seed)):
        if abs(eta_seed[index]) > 3.0: continue
        if pt_seed[index] < 26: continue
        
        hoe_cut = (h_for_hoe_seed[index] / energy_seed[index]) < .1
        r9_cut = r9_seed[index] > cut
        
        #sieie varies by eta
        sieie_cut = False
        if abs(eta_seed[index]) < 1.5:
            sieie_cut = sieie_seed[index] < .014
        else:
            sieie_cut = sieie_seed[index] < .035
        
        ecal_cut = ecal_iso_seed[index] < (5 + .012 * pt_seed[index])
        hcal_cut = hcal_iso_seed[index] < (5 + .005 * pt_seed[index])
        #track_cut = track_iso[index] < (5 + .002 * pt[index])
        
        calo_cut = hoe_cut and sieie_cut
        iso_cut = ecal_cut and hcal_cut
        
        if r9_cut or (calo_cut and iso_cut):
            found = True
            break

    return found

#photon id 
def photonPassID(index,cut):
    if abs(eta[index]) > 3.0: return False

    hoe_cut = (h_for_hoe[index] / energy[index]) < .1
    r9_cut = r9[index] > cut
    
    #sieie varies by eta
    sieie_cut = False
    if abs(eta[index]) < 1.5:
        sieie_cut = sieie[index] < .014
    else:
        sieie_cut = sieie[index] < .035
        
    ecal_cut = ecal_iso[index] < (5 + .012 * pt[index])
    hcal_cut = hcal_iso[index] < (5 + .005 * pt[index])
    track_cut = track_iso[index] < (5 + .002 * pt[index])

    calo_cut = hoe_cut and sieie_cut
    iso_cut = ecal_cut and hcal_cut and track_cut

    return r9_cut or (calo_cut and iso_cut)

def setBranchAddresses(tree):
    tree.SetBranchAddress("hltL1SeededHLTClusterShape", sieie_seed)
    tree.SetBranchAddress("hltL1SeededPhotonEcalIso", ecal_iso_seed)
    tree.SetBranchAddress("hltL1SeededPhotonHcalForHE", h_for_hoe_seed)
    tree.SetBranchAddress("hltL1SeededPhotonHcalIso", hcal_iso_seed)
    tree.SetBranchAddress("hltL1SeededR9ID", r9_seed)

    tree.SetBranchAddress("hltActivityPhotonClusterShape", sieie)
    tree.SetBranchAddress("hltActivityPhotonEcalIso", ecal_iso)
    tree.SetBranchAddress("hltActivityPhotonHcalForHE", h_for_hoe)
    tree.SetBranchAddress("hltActivityPhotonHcalIso", hcal_iso) 
    tree.SetBranchAddress("hltActivityPhotonHollowTrackIsoWithId", track_iso)
    tree.SetBranchAddress("hltActivityR9ID", r9)

    tree.SetBranchAddress("pt", pt)             
    tree.SetBranchAddress("eta", eta)            
    tree.SetBranchAddress("phi", phi)            
    tree.SetBranchAddress("energy", energy)

    tree.SetBranchAddress("pt_seed", pt_seed)             
    tree.SetBranchAddress("eta_seed", eta_seed)            
    tree.SetBranchAddress("phi_seed", phi_seed)            
    tree.SetBranchAddress("energy_seed", energy_seed)           

cut_vals = [.7, .75, .8, .85, .9, .95]

#treedata =  data.Get("hltTree")
tree20x25 = pileup20x25.Get("hltTree")
tree40x25 = pileup40x25.Get("hltTree")
tree40x50 = pileup40x50.Get("hltTree")
tree80x50 = pileup80x50.Get("hltTree")



outfile =  rt.TFile("outfile.root","RECREATE")

hist_20x25 = rt.TH1F("h20x25","PU: 20 bx: 25 ns", len(cut_vals)-1, array.array("d",cut_vals))
hist_40x25 = rt.TH1F("h40x25","PU: 40 bx: 25 ns",len(cut_vals)-1, array.array("d",cut_vals))
hist_40x50 = rt.TH1F("h40x50","PU: 40 bx: 50 ns",len(cut_vals)-1, array.array("d",cut_vals))
hist_80x50 = rt.TH1F("h80x50","PU: 80 bx: 50 ns",len(cut_vals)-1, array.array("d",cut_vals))
#hist_data = rt.TH1F("hdata","8 TeV Run 208390",len(cut_vals)-1, array.array("d",cut_vals))

trees = [tree20x25, tree40x25, tree40x50, tree80x50]#, treedata]
hists = [hist_20x25, hist_40x25, hist_40x50, hist_80x50]#, hist_data]

#loop ove rthe pairs of trees and histograms
for pair in zip(trees,hists):
    tt = pair[0]
    hh = pair[1]

    print tt

    #set the branch addresses for the given tree    
    setBranchAddresses(tt)
    #fill the histograms
    for ii in range(1000):
        if ii % 100 == 0 : print "Scanning Event...", ii
        tt.GetEntry(ii)
        for cut in cut_vals[:-1]:
            (high, low, ihigh, ilow) = getSelectionPair(cut)
            if ihigh != -99 and ilow != -99: #both photons pass
                hh.Fill(cut, 1/1000.)

for h in hists: h.Write()
outfile.Close()

os.system("python draw.py")
