"""
================================================================================
An example configuration file for plotter.py

Author:
    Alex Armstrong <alarmstr@cern.ch>

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
import atlasrootstyle.AtlasStyle
ROOT.SetAtlasStyle()

################################################################################
# Useful globals
################################################################################

eta_gt_1p45 = False
eta_lt_1p45 = True
njet_0 = True
njet_1 = False
njet_ge2 = False

NONCLOSURE_SYS = 0.30

# Strings for plotting


################################################################################
# Samples
################################################################################
from PlotTools.fake_factor.sample_conf import *


################################################################################
# Regions
################################################################################
from PlotTools.fake_factor.region_conf import *

zjets_FF_CR = '1'
zjets_FF_CR += ' && %s'%singlelep_trig_pT
zjets_FF_CR += ' && fabs(lep_d0sigBSCorr[0]) < 15 && fabs(lep_d0sigBSCorr[1]) < 15 && fabs(lep_d0sigBSCorr[2]) < 15'
zjets_FF_CR += ' && fabs(lep_z0SinTheta[0]) < 15 && fabs(lep_z0SinTheta[1]) < 15 && fabs(lep_z0SinTheta[2]) < 15'
zjets_FF_CR += ' && (80 < Z_MLL && Z_MLL < 100)'
zjets_FF_CR += ' && l_mT[2] < 40'
zjets_FF_CR += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
zjets_FF_CR += ' && MET < 60'
#zjets_FF_CR += ' && nBJets == 0'
#zjets_FF_CR += ' && nLJets <= 2'
if eta_lt_1p45:
    zjets_FF_CR += ' && fabs(l_eta[2]) < 1.45'
elif eta_gt_1p45:
    zjets_FF_CR += ' && fabs(l_eta[2]) >= 1.45'

if njet_0:
    zjets_FF_CR += ' && n_jets == 0'
elif njet_1:
    zjets_FF_CR += ' && n_jets == 1'
elif njet_ge2:
    zjets_FF_CR += ' && n_jets >= 2'

zjets_FF_truth_base = '(!isMC || (0 < l_truthClass[0] && l_truthClass[0] <= 2))' #Prompt Leading Lepton
zjets_FF_truth_base += ' && (!isMC || (0 < l_truthClass[1] && l_truthClass[1] <= 2))' #Prompt Subleading Lepton
zjets_FF_truth_num = zjets_FF_truth_base\
                   +' && (!isMC || (0 < l_truthClass[2] && l_truthClass[2] <= 2))' #Prompt Probe Lepton
zjets_FF_truth_den = zjets_FF_truth_base\
                   + ' && (!isMC || (l_truthClass[2] <= 0 || 2 < l_truthClass[2]))' #Fake Probe Lepton
num_den_dict = {'den' : 'nLepID == 2 && nLepAntiID >= 1',
                'num' : 'nLepID == 3'}
chan_dict = {'eee' : ['ee','e','Z_dilep_flav==2 && l_flav[2]==0'],
             'mme' : ['mumu','e','Z_dilep_flav==3 && l_flav[2]==0'],
             'eem' : ['ee','m','Z_dilep_flav==2 && l_flav[2]==1'],
             'mmm' : ['mumu','m','Z_dilep_flav==3 && l_flav[2]==1'],
             'm' : ['ll','m','l_flav[2]==1'],
             'e' : ['ll','e','l_flav[2]==0'],
        }
for num_den, num_den_sel in num_den_dict.iteritems():
    for chan, ops in chan_dict.iteritems():
        id_or_aid = 'ID' if num_den == 'num' else 'anti-ID'
        chan_name = '%s + %s %s'%(ops[0], id_or_aid, ops[1])

        name = 'zjets_FF_CR%s_%s'%(num_den, chan)
        if eta_lt_1p45:
            chan_name += ', |#eta| < 1.4'
        elif eta_gt_1p45:
            chan_name += ', |#eta| #ge 1.4'
        if njet_0:
            chan_name += ', N_{jets} = 0'
        elif njet_1:
            chan_name += ', N_{jets} = 1'
        elif njet_ge2:
            chan_name += ', N_{jets} #ge 2'

        displayname = 'Z+jets CR (%s)'%(chan_name)

        REGIONS.append(Region(name, displayname))
        REGIONS[-1].tcut = ' && '.join([ops[2], zjets_FF_CR])
        REGIONS[-1].truth_fake_sel = zjets_FF_truth_den
        REGIONS[-1].truth_bkg_sel = zjets_FF_truth_num



################################################################################
# Variables
################################################################################
from PlotTools.fake_factor.plot_conf import * #TODO: Dont use import *

################################################################################
# Toggle options for execution
# - Options expected to change often between different plottings
################################################################################
region_plots = {}


#######################################
# Yield Table
from PlotTools.YieldTable import YieldTable
YIELD_TBL = YieldTable()
# Add formulas to yield table
# Key values will be the labels on the table
# Formulas should use sample names and these symbols: +, -, *, /, (, ), [0-9]
#YIELD_TBL.formulas['Data/MC'] = "data/MC"
#YIELD_TBL.formulas['1 - (WZ+ZZ)/Data'] = "1 - (wz + zz)/data"
#YIELD_TBL.formulas['Zll/Data'] = "zll/data"

#######################################
# What regions to plot
region_ops = []
region_ops += ['zjets_FF_CRden_e', 'zjets_FF_CRnum_e']
region_ops += ['zjets_FF_CRden_m', 'zjets_FF_CRnum_m']
#region_ops += ['zjets_FF_CRden_eem', 'zjets_FF_CRnum_eem']
#region_ops += ['zjets_FF_CRden_mmm', 'zjets_FF_CRnum_mmm']
#region_ops += ['zjets_FF_CRden_eee', 'zjets_FF_CRnum eee']
#region_ops += ['zjets_FF_CRden_mme', 'zjets_FF_CRnum_mme']
#region_ops += ['zjets_FF_CRden_eem']
#######################################
# What variables to plot
vars_to_plot = []
if eta_lt_1p45 or eta_gt_1p45 or njet_0 or njet_1 or njet_ge2:
    vars_to_plot += ['l_pt[2]']
else:
    vars_to_plot += ['l_pt[2]']
    #vars_to_plot = ['fabs(l_eta[2]):l_pt[2]']
    #vars_to_plot = ['n_jets:fabs(l_eta[2]):l_pt[2]']

# Remove duplicate names
vars_to_plot = list(set(vars_to_plot))

################################################################################
# Create plots
################################################################################
PLOTS = []
for var in vars_to_plot:

    # Use defualt settings, if needed, when plotting a specific element of a
    # vector (i.e. lepton_flav[1] uses lepton_flav settings if lepton_flav[1]
    # is not defined separately)
    if '[' in var and var not in plot_defaults:
        key = var.split('[')[0]
    else:
        key = var

    # Create plot for each region
    for region in region_ops:

        # Grab the default plot unless a region specific one is defined
        if region in region_plots and key in region_plots[region]:
            p = deepcopy(region_plots[region][key])
        elif key in plot_defaults:
            p = deepcopy(plot_defaults[key])
        else:
            assert False, ("ERROR :: requested plot not defined:", var)

        if p.is3D:
            varz, vary, varx = var.split(':')
            p.update(region, varx, vary, varz)
        elif p.is2D:
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

            if n_bkgds and n_data:
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
        elif p.ptype == Types.three_dim:
            pass #Not intended for plotting
        else:
            print "WARNING :: %s plots are not yet setup"%p.ptype.name
            continue

        PLOTS.append(p)

print "Configuration File End"
