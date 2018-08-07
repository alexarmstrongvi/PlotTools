"""
================================================================================
An example configuration file for plotter.py

Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with lots of ideas borrowed from Danny Antrim <dantrim@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""

# General python
import sys, os
from copy import deepcopy

# Root data analysis framework
import ROOT
# Prevent root from printing anything to the screen
ROOT.gROOT.SetBatch(ROOT.kTRUE)

################################################################################
# Useful globals
################################################################################
import global_variables as g
# Toggles
run_fakes = False
add_truth_den = False
add_truth_num = False
run_den = False
run_num = True

assert not run_fakes or (add_truth_num)
assert not (run_num and run_den)
assert not (add_truth_num and add_truth_den)

# Path to directories with flat ntuples
if run_num:
    bkg_ntuple_dir    = g.mc_num_ntuples
    signal_ntuple_dir = g.mc_num_ntuples
    data_ntuple_dir   = g.data_num_ntuples
elif run_den:
    bkg_ntuple_dir    = g.mc_den_ntuples
    signal_ntuple_dir = g.mc_den_ntuples
    data_ntuple_dir   = g.data_den_ntuples
else:
    bkg_ntuple_dir    = g.mc_ntuples
    signal_ntuple_dir = g.mc_ntuples
    data_ntuple_dir   = g.data_ntuples

if run_fakes:
    fake_ntuple_dir   = g.fake_ntuples



# Luminosity options
# Description | 2015-16 (ipb) |  2015-16 (ifb) | 2015 (ipb) | 2016 (ipb) | dummy
# Value       |     36180     |      36.01     |    320     |    32971   |   1
lumi = 36180

from example_sample_conf import *
#######################################
## Build the TChain/TTree for each sample
# To remove sample from plot, comment out the line setting its TChain
# Samples with empty TChains get removed below
data_dsids = g.groups['data15']+g.groups['data16']

data.set_chain_from_dsid_list(data_dsids, data_ntuple_dir)
top.set_chain_from_dsid_list(g.groups['top'], bkg_ntuple_dir)
#ttbar.set_chain_from_dsid_list(g.groups['ttbar'], bkg_ntuple_dir)
#stop.set_chain_from_dsid_list(g.groups['singletop'], bkg_ntuple_dir)
#wtop.set_chain_from_dsid_list(g.groups['Wt'], bkg_ntuple_dir)
ttbarlep.set_chain_from_dsid_list(g.groups['ttbar_lep'], bkg_ntuple_dir)
VV.set_chain_from_dsid_list(g.groups['VV'], bkg_ntuple_dir)
VVV.set_chain_from_dsid_list(g.groups['VVV'], bkg_ntuple_dir)
#zll.set_chain_from_dsid_list(g.groups['zll'], bkg_ntuple_dir)
zee.set_chain_from_dsid_list(g.groups['zee'], bkg_ntuple_dir)
zmumu.set_chain_from_dsid_list(g.groups['zmumu'], bkg_ntuple_dir)
ztt.set_chain_from_dsid_list(g.groups['ztt'], bkg_ntuple_dir)
wjets.set_chain_from_dsid_list(g.groups['wjets'], bkg_ntuple_dir)
wgamma.set_chain_from_dsid_list(g.groups['wgamma'], bkg_ntuple_dir)
zgamma.set_chain_from_dsid_list(g.groups['zgamma'], bkg_ntuple_dir)
higgs.set_chain_from_dsid_list(g.groups['higgs'], bkg_ntuple_dir)
#signal.set_chain_from_dsid_list(g.groups['higgs_lfv'], signal_ntuple_dir)
if run_fakes:
    fakes.set_chain_from_dsid_list(data_dsids, fake_ntuple_dir)

SAMPLES = [s for s in SAMPLES if s.is_setup()]
assert SAMPLES, "ERROR :: No samples are setup"

################################################################################
# Regions
# - Regions defined at runtime can be moved from region_conf to here
################################################################################
from example_region_conf import *

# Fake Factor estimation region: Z+Jets regions
zjets_FF_CR = '1'
zjets_FF_CR += ' && %s'%singlelep_trig_pT
zjets_FF_CR += ' && fabs(lep_d0sigBSCorr[0]) < 15 && fabs(lep_d0sigBSCorr[1]) < 15 && fabs(lep_d0sigBSCorr[2]) < 15'
zjets_FF_CR += ' && fabs(lep_z0SinTheta[0]) < 15 && fabs(lep_z0SinTheta[1]) < 15 && fabs(lep_z0SinTheta[2]) < 15'
zjets_FF_CR += ' && (80 < Z_MLL && Z_MLL < 100)'
zjets_FF_CR += ' && nBJets == 0'
zjets_FF_CR += ' && l_mT[2] < 50'
zjets_FF_CR += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
#zjets_FF_CR += ' && MET < 50'
if run_fakes:
    zjets_FF_CR += ' && (!isMC || nLepID == 2 || (0 < l_truthClass[0] && l_truthClass[0] <= 2))' #Prompt Leading Lepton
    zjets_FF_CR += ' && (!isMC || nLepID == 2 || (0 < l_truthClass[1] && l_truthClass[1] <= 2))' #Prompt Subleading Lepton
    zjets_FF_CR += ' && (!isMC || nLepID == 2 || (0 < l_truthClass[2] && l_truthClass[2] <= 2))' #Prompt Probe Lepton
if add_truth_den or add_truth_num:
    zjets_FF_CR += ' && (!isMC || (0 < l_truthClass[0] && l_truthClass[0] <= 2))' #Prompt Leading Lepton
    zjets_FF_CR += ' && (!isMC || (0 < l_truthClass[1] && l_truthClass[1] <= 2))' #Prompt Subleading Lepton
    if add_truth_num:
        zjets_FF_CR += ' && (!isMC || (0 < l_truthClass[2] && l_truthClass[2] <= 2))' #Prompt Probe Lepton
    else:
        zjets_FF_CR += ' && (!isMC || (l_truthClass[2] <= 0 || 2 < l_truthClass[2]))' #Fake Probe Lepton

REGIONS.append(Region('zjets_FF_CR', 'Z+jets FF CR'))
REGIONS[-1].tcut = zjets_FF_CR

num_den_list = ['num', 'den']
chan_dict = {'eee' : ['ee','e','Z_dilep_flav==2 && l_flav[2]==0'],
             'mme' : ['mumu','e','Z_dilep_flav==3 && l_flav[2]==0'],
             'eem' : ['ee','m','Z_dilep_flav==2 && l_flav[2]==1'],
             'mmm' : ['mumu','m','Z_dilep_flav==3 && l_flav[2]==1'],
             'm' : ['ll','m','l_flav[2]==1'],
             'e' : ['ll','e','l_flav[2]==0'],
        }
for num_den in num_den_list:
    for chan, ops in chan_dict.iteritems():
        id_or_aid = 'ID' if num_den == 'num' else 'anti-ID'
        chan_name = '%s + %s %s'%(ops[0], id_or_aid, ops[1])

        name = 'zjets_FF_CR%s_%s'%(num_den, chan)
        displayname = 'Z+jets FF CR (%s)'%(chan_name)
        REGIONS.append(Region(name, displayname))
        REGIONS[-1].tcut = ' && '.join([ops[2], zjets_FF_CR])

# Fake Factor validation region: Wjet
wjets_FF_VR = singlelep_trig_pT
wjets_FF_VR += ' && fabs(lep_d0sigBSCorr[0]) < 15 && fabs(lep_d0sigBSCorr[1]) < 15'
wjets_FF_VR += ' && fabs(lep_z0SinTheta[0]) < 15 && fabs(lep_z0SinTheta[1]) < 15'
wjets_FF_VR += ' && MLL > 10'
wjets_FF_VR += ' && nBJets==0' # reject top
wjets_FF_VR += ' && nLJets <= 1'
wjets_FF_VR += ' && l_mT[0] > 50'
wjets_FF_VR += ' && RelMET < 60'
wjets_FF_VR += ' && dpt_ll > 15'
wjets_FF_VR += ' && MLL < 110'
wjets_FF_VR += ' && n_jets == 0'
#wjets_FF_VR += ' && n_jets == 1'

if add_truth_den or add_truth_num:
    wjets_FF_VR += ' && (!isMC || (0 < l_truthClass[0] && l_truthClass[0] <= 2))' #Prompt Leading Lepton
    if add_truth_num:
        wjets_FF_VR += ' && (!isMC || (0 < l_truthClass[1] && l_truthClass[1] <= 2))' #Prompt Probe Lepton
    else:
        wjets_FF_VR += ' && (!isMC || (l_truthClass[1] <= 0 || 2 < l_truthClass[1]))' #Fake Probe Lepton

REGIONS.append(Region('wjets_FF_VR', 'wjets FF VR'))
REGIONS[-1].tcut = wjets_FF_VR

REGIONS.append(Region('wjets_FF_VRden_emu', 'wjets FF VR (anti-ID mu)'))
REGIONS[-1].tcut = wjets_FF_VR +' && ' + emu

REGIONS.append(Region('wjets_FF_VRden_mue', 'wjets FF VR (anti-ID el)'))
REGIONS[-1].tcut = wjets_FF_VR + ' && ' + mue

REGIONS.append(Region('wjets_FF_VRnum_emu', 'wjets FF VR (ID mu)'))
REGIONS[-1].tcut = wjets_FF_VR +' && ' + emu

REGIONS.append(Region('wjets_FF_VRnum_mue', 'wjets FF VR (ID el)'))
REGIONS[-1].tcut = wjets_FF_VR + ' && ' + mue

region_plots = {}
#region_plots['zjets_FF_CRden_eem'] = {
#  'Z_MLL' : deepcopy(plot_defaults['Z_MLL'])
#}
#region_plots['zjets_FF_CRden_eem']['Z_MLL'].update(bin_range=[80, 100], bin_width=1)

#######################################
# What regions to plot
region_ops = []
if run_den:
    #region_ops += ['wzCR']
    #region_ops += ['zjets_FF_CRden_m'] # Test Region
    #region_ops += ['zjets_FF_CRden_eee'] # Test Region
    #region_ops += ['zjets_FF_CRden_m', 'zjets_FF_CRden_e']
    #region_ops += ['zjets_FF_CRden_eem'] # Test region
    #region_ops += ['zjets_FF_CRden_eem', 'zjets_FF_CRden_mmm']
    #region_ops += ['zjets_FF_CRden_eee', 'zjets_FF_CRden_mme']
    region_ops += ['wjets_FF_VRden_emu']
elif run_num:
    region_ops += ['wzCR']
    #region_ops += ['zjets_FF_CRnum_m']
    #region_ops += ['zjets_FF_CRnum_m', 'zjets_FF_CRnum_e']
    #region_ops += ['zjets_FF_CRnum_eem', 'zjets_FF_CRnum_mmm']
    #region_ops += ['zjets_FF_CRnum_eee', 'zjets_FF_CRnum_mme']
    #region_ops += ['wjets_FF_VRnum_emu'] #, 'wjets_FF_VRnum_mue'
else:
    #region_ops += ['zjets_FF_CR']
    region_ops += ['zCR']


################################################################################
# Toggle options for execution
# - Options expected to change often between different plottings
################################################################################

#######################################
# Yield Table
from tools.YieldTable import YieldTable

YIELD_TBL = YieldTable()
# Add formulas to yield table
# Key values will be the labels on the table
# Formulas should use sample names and these symbols: +, -, *, /, (, ), [0-9]
# 'MC' is also a recognized label for all backgrounds
YIELD_TBL.formulas['Data/MC'] = "data/MC"
YIELD_TBL.formulas['VV/MC'] = "vv/MC"
YIELD_TBL.formulas['normF'] = "(data-(MC-vv))/vv"
YIELD_TBL.formulas['VV/sqrt(MC)'] = "vv/(MC**(0.5))"
#YIELD_TBL.formulas['WZ/sqrt(Data)'] = "wz/(data**(0.5))"
#YIELD_TBL.formulas['Bkg/Data'] = "(MC-wjets)/data"
#YIELD_TBL.formulas['W+Jets/MC'] = "wjets/MC"
#YIELD_TBL.formulas['Fakes/MC'] = "fakes/MC"
#YIELD_TBL.formulas['Fakes/Data'] = "fakes/data"


#######################################
# What variables to plot
# Strings for plotting
from example_plot_conf import *
PlotBase.save_dir = g.plots_dir
Plot1D.auto_set_ylimits = True
Plot1D.doLogY = False
Plot2D.doLogZ = False
Plot2D.auto_set_zlimits = False

#vars_to_plot = ['isMC','l_pt[0]','l_pt[0]:l_pt[2]']
vars_to_plot = ['l_pt[0]:l_pt[2]']

# Remove duplicate names
vars_to_plot = list(set(vars_to_plot))
assert vars_to_plot, ("No plots requested")

################################################################################
# Create plots
################################################################################
PLOTS = []
for var in vars_to_plot:
    for region in region_ops:
        key = var
        # Grab the default plot unless a region specific one is defined
        if region in region_plots and key in region_plots[region]:
            p = deepcopy(region_plots[region][key])
        elif key in plot_defaults:
            p = deepcopy(plot_defaults[key])
        else:
            assert False, ("ERROR :: requested plot not defined:", var)

        # Defing the region and variables
        if p.is2D:
            varx = var.split(':')[1]
            vary = var.split(':')[0]
            p.update(region, varx, vary)
        else:
            p.update(region, var)

        # Set plot type if not already set
        if p.ptype == Types.default:
            n_bkgds = len([s for s in SAMPLES if s.isMC and not s.isSignal])
            n_signal = len([s for s in SAMPLES if s.isMC and s.isSignal])
            n_data = len([s for s in SAMPLES if not s.isMC])
            assert n_data <= 1, "ERROR :: More than one data sample setup"

            if n_bkgds and n_data and not p.doNorm:
                p.ptype = Types.ratio
            elif n_bkgds or n_data:
                 p.ptype = Types.stack
            else:
                assert False, "ERROR :: No samples are setup"

        # Setup correct canvas for plotting
        if p.ptype == Types.ratio:
            p.setRatioPads(p.name)
        elif p.ptype == Types.stack:
            p.setStackPads(p.name)
        elif p.ptype == Types.two_dim:
            p.set2DPads(p.name)
        else:
            print "WARNING :: %s plots are not yet setup"%p.ptype.name
            continue

        PLOTS.append(p)

print "Configuration File End"
