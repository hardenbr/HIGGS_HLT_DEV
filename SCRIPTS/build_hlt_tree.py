from  optparse  import OptionParser
import ROOT as rt
from DataFormats.FWLite import Events, Handle
import os, sys

parser = OptionParser()

default_file = "/afs/cern.ch/work/h/hardenbr/2013/HIGGS_DIPHOTON_HLT/fscan_mass100_testfile/res/temp_27_19.root"

parser.add_option("-f", "--file", dest="filename",
                  help="hlt.root file name to analyze FILE",default=default_file,
                  action="store",type="string")

parser.add_option("-o", "--outfile", dest="outfilename",
                  help="tree.root file name to output",default="hlt_tree.root",
                  action="store",type="string")

parser.add_option("-t", "--tag", dest="tag",
                  help="which HLT Process",
                  action="store",type="string", default="HLT3")

(options, args) = parser.parse_args()

parser.print_help()

HLT_TAG = "HLT3"

###WRITE OUT ALL OF THE HANDLE INFORMATION###
cluster = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltActivityPhotonClusterShape"   ,"",                HLT_TAG)

ecal_iso = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltActivityPhotonEcalIso"   ,"",                HLT_TAG)

hoe = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltActivityPhotonHcalForHE"   ,"",                HLT_TAG)

hcal_iso = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltActivityPhotonHcalIso"   ,"",                HLT_TAG)

track_iso = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltActivityPhotonHollowTrackIsoWithId"   ,"",                HLT_TAG)

r9 = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltActivityR9ID"           ,"",                HLT_TAG)

###SEEDED INFORMATION###
cluster_seed = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltL1SeededHLTClusterShape","", HLT_TAG)

ecal_iso_seed = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltL1SeededPhotonEcalIso"   ,"",                HLT_TAG)

hoe_seed = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltL1SeededPhotonHcalForHE"   ,"",                HLT_TAG)

hcal_iso_seed = ("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltL1SeededPhotonHcalIso"   ,"",                HLT_TAG)

r9_seed =("edm::AssociationMap<edm::OneToValue<vector<reco::RecoEcalCandidate>,float,unsigned int> >"    ,"hltL1SeededR9ID"           ,"",                HLT_TAG)


##SUMMARY
seed = (cluster_seed, ecal_iso_seed, hoe_seed, hcal_iso_seed, r9_seed)
non_seed = (cluster, ecal_iso, hoe, hcal_iso, track_iso, r9)
handles= []

#build the vectors

#fill the handle array (label, handle)
for ii in seed+non_seed:
    handle = Handle(ii[0])
    label = ii[1]
    middle = ii[2]
    process = ii[3]

    handles.append([(label,middle,process),handle])
        
rt.gROOT.SetBatch()
counter = 0

#build the vectors for the tree
vectors = []
vectors_cand = []
for ii in handles:
    vector = rt.vector("float")()
    vectors.append((ii[0][0],vector)) #(label, vector)

candidates = ["pt","eta","phi","energy"]

#kinematic, vector pairs
for ii in candidates:
    vector = rt.vector("float")()
    vectors_cand.append((ii,vector))

events = Events(options.filename)
tree = rt.TTree("hltTree", "hltTree")

#set the branches in the tree
for ii in vectors:
    tree.Branch(ii[0],ii[1])

for ii in vectors_cand:
    tree.Branch(ii[0],ii[1])


# Create a struct for the run information
rt.gROOT.ProcessLine(\
      "struct MyStruct{\
      Int_t run;\
      Int_t ls;\
      Int_t event;\
      };")

from ROOT import MyStruct

# Create branches in the tree
s = MyStruct()

tree.Branch("run",rt.AddressOf(s,"run"),"run/I")
tree.Branch("ls",rt.AddressOf(s,"ls"),"ls/I")
tree.Branch("event",rt.AddressOf(s,"event"),"event/I")

#retreives the lables for a handle for an event
def getAllLabels(event):
    for ii in handles:
        input_tag = ii[0] #triplet (label, middle, process(
        handle = ii[1] 
        event.getByLabel(input_tag, handle)

#gets the maps for all of the handles and returns a list
def getMaps():
    maps = []

    for ii in handles:
        product = ii[1].product() #handle.product
        name = ii[0][0] #label        
        maps.append((name,product))

    return maps

#get the events from the file 
for event in events:    
    if counter % 1000 == 0: print "Analyzing " + str(counter) + "..."

    getAllLabels(event)
    maps = getMaps() # list of (name, product)

    aux = event.eventAuxiliary()
    s.run = aux.run()
    s.ls = aux.luminosityBlock()
    s.event = aux.event()
    
    for map in maps:
        name = map[0]        
        values = map[1].values()
        cand = map[1].refProd().key.product()        
        handle_ii = maps.index(map)
        
        for ii in range(values.size()):
            vectors[handle_ii][1].push_back(values[ii])

        if cand != None:
            for ii in range(values.size()): # for each entry
                for vec in vectors_cand: #loop over pt, eta, phi ect..
                    name = vec[0]
                    val = eval("cand[%i].%s()" % (ii,name)) # fill the coressponding val                
                    if maps.index(map) == 9: vec[1].push_back(val)  #9 is an index for activity handle
                    # if we used a different index we would only get kinematics for those bits

    #fill tree and clear the vectors out
    tree.Fill()
    
    for ii in vectors:
        ii[1].clear()

    for ii in vectors_cand:
        ii[1].clear()
                        
    counter+=1
    continue

outfile = rt.TFile(options.outfilename,"RECREATE")
tree.CloneTree().Write()
outfile.Close()



