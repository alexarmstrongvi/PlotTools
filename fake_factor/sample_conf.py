from PlotTools.sample import Sample, MCsample, Data, Background, Signal
import ROOT

import global_variables as g
# Path to directories with flat ntuples
bkg_ntuple_dir    = g.mc_ntuples
signal_ntuple_dir = g.mc_ntuples
data_ntuple_dir   = g.data_ntuples
fake_ntuple_dir   = g.data_ntuples

NUM_STR = "num"
DEN_STR = "den"

################################################################################
# Samples
# - Create all the samples that will be put into plots (e.g. backgrounds,
#   signal, data)
################################################################################
SAMPLES = []

# Luminosity options
# Description | 2015-16 (ipb) |  2015-16 (ifb) | 2015 (ipb) | 2016 (ipb) | dummy
# Value       |     36180     |      36.01     |    320     |    32971   |   1
lumi = 36180

MCsample.weight_str = 'eventweight'
MCsample.scale_factor = lumi
Sample.input_file_treename = 'superNt'

for num_den in [NUM_STR, DEN_STR]:
    # Data
    data = Data('data_%s'%num_den)
    SAMPLES.append(data)

    #######################################
    # Initialize all backgrounds
    # Higgs (Htt + HWW + HV + ttH)
    higgs = Background("higgs_%s"%num_den, "Higgs")
    higgs.color = ROOT.kRed - 7
    SAMPLES.append(higgs)

    # ttbar
    ttbar = Background("ttbar_%s"%num_den, "t#bar{t}")
    ttbar.color = ROOT.kOrange+2
    SAMPLES.append(ttbar)

    # singletop
    stop = Background("st_%s"%num_den, "Single top")
    stop.color = ROOT.kOrange+1
    SAMPLES.append(stop)

    # W+top
    wtop = Background("wt_%s"%num_den, "Wt")
    wtop.color = ROOT.kOrange+8
    SAMPLES.append(wtop)

    # ttbar + X
    ttbar_x = Background("ttbar_x_%s"%num_den, "t#bar{t} + X")
    ttbar_x.color = ROOT.kOrange+5
    SAMPLES.append(ttbar_x)

    # Top (ttbar + singletop + Wt)
    top = Background("top_%s"%num_den, "Top")
    top.color = ROOT.kYellow+1
    SAMPLES.append(top)
    
    # VV combined
    VV = Background("VV_%s"%num_den, "VV")
    VV.color = ROOT.kSpring-7
    VV.scale_factor *= 1.05
    SAMPLES.append(VV)

    # VVV combined
    VVV = Background("VVV_%s"%num_den, "VVV")
    VVV.color = ROOT.kSpring-7
    SAMPLES.append(VVV)

    # Zll
    zll = Background("zll_%s"%num_den, "Zll")
    zll.color = ROOT.kAzure-7
    SAMPLES.append(zll)

    # Zee
    zee = Background("zee_%s"%num_den, "Zee")
    zee.color = ROOT.kAzure-7
    SAMPLES.append(zee)

    # Zmumu
    zmumu = Background("zmumu_%s"%num_den, "Zmumu")
    zmumu.color = ROOT.kAzure-9
    SAMPLES.append(zmumu)

    # Ztt
    ztt = Background("ztt_%s"%num_den, "Z#tau#tau")
    ztt.color = ROOT.kAzure+8
    SAMPLES.append(ztt)

    # Wjets
    wjets = Background("wjets_%s"%num_den, "W+jets")
    wjets.color = ROOT.kOrange + 2
    SAMPLES.append(wjets)

    # W+gamma
    wgamma = Background("wgamma_%s"%num_den, "W+gamma")
    wgamma.color = ROOT.kOrange-1
    SAMPLES.append(wgamma)

    # Z+gamma
    zgamma = Background("zgamma_%s"%num_den, "Z+gamma")
    zgamma.color = ROOT.kOrange-3
    SAMPLES.append(zgamma)

    #######################################
    ## Build the TChain/TTree for each sample
    # To remove sample from plot, comment out the line setting its TChain
    # Samples with empty TChains get removed below
    data_ntuple_dir2 = data_ntuple_dir[:-1] + '_' + num_den + "/"
    bkg_ntuple_dir2 = bkg_ntuple_dir[:-1] + '_' + num_den + "/"

    data.set_chain_from_dsid_list(g.groups['data'], data_ntuple_dir2)
    #top.set_chain_from_dsid_list(g.groups['top'], bkg_ntuple_dir2)
    ###ttbar.set_chain_from_dsid_list(g.groups['ttbar'], bkg_ntuple_dir2)
    ###stop.set_chain_from_dsid_list(g.groups['singletop'], bkg_ntuple_dir2)
    ###wtop.set_chain_from_dsid_list(g.groups['Wt'], bkg_ntuple_dir2)
    #ttbar_x.set_chain_from_dsid_list(g.groups['ttbarX'], bkg_ntuple_dir2)
    #VV.set_chain_from_dsid_list(g.groups['VV'], bkg_ntuple_dir2)
    VVV.set_chain_from_dsid_list(g.groups['VVV'], bkg_ntuple_dir2)
    #zll.set_chain_from_dsid_list(g.groups['zll'], bkg_ntuple_dir2)
    zee.set_chain_from_dsid_list(g.groups['zee'], bkg_ntuple_dir2)
    ###zmumu.set_chain_from_dsid_list(g.groups['zmumu'], bkg_ntuple_dir2)
    #ztt.set_chain_from_dsid_list(g.groups['ztt'], bkg_ntuple_dir2)
    #wjets.set_chain_from_dsid_list(g.groups['wjets'], bkg_ntuple_dir2)
    #wgamma.set_chain_from_dsid_list(g.groups['wgamma'], bkg_ntuple_dir2)
    #zgamma.set_chain_from_dsid_list(g.groups['zgamma'], bkg_ntuple_dir2)
    #higgs.set_chain_from_dsid_list(g.groups['higgs'], bkg_ntuple_dir2)

SAMPLES = [s for s in SAMPLES if s.is_setup()]
assert SAMPLES, "ERROR :: No samples are setup"

