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
add_truth_den = True
add_truth_num = False
run_den = False
run_num = True
run_zjets = True
run_base = False

assert run_den != run_num
assert run_zjets != run_base
assert not run_fakes or add_truth_num
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

data.set_chain_from_dsid_list(g.groups['data'], data_ntuple_dir)
top.set_chain_from_dsid_list(g.groups['top'], bkg_ntuple_dir)
#ttbar.set_chain_from_dsid_list(g.groups['ttbar'], bkg_ntuple_dir)
#stop.set_chain_from_dsid_list(g.groups['singletop'], bkg_ntuple_dir)
#wtop.set_chain_from_dsid_list(g.groups['Wt'], bkg_ntuple_dir)

ttbarX.set_chain_from_dsid_list(g.groups['ttbarX'], bkg_ntuple_dir)

VV.set_chain_from_dsid_list(g.groups['VV'], bkg_ntuple_dir)

VVV.set_chain_from_dsid_list(g.groups['VVV'], bkg_ntuple_dir)

zll.set_chain_from_dsid_list(g.groups['zll'], bkg_ntuple_dir)
#zee.set_chain_from_dsid_list(g.groups['zee'], bkg_ntuple_dir)
#zmumu.set_chain_from_dsid_list(g.groups['zmumu'], bkg_ntuple_dir)

ztt.set_chain_from_dsid_list(g.groups['ztt'], bkg_ntuple_dir)
#ztt_sherpa.set_chain_from_dsid_list(g.groups['ztt_sherpa'], bkg_ntuple_dir)
#ztt_lowMLL.set_chain_from_dsid_list(g.groups['ztt_lowMLL'], bkg_ntuple_dir)
#ztt_2jets.set_chain_from_dsid_list(g.groups['ztt_2jets'], bkg_ntuple_dir)
##ztt_l13l7.set_chain_from_dsid_list(g.groups['ztt_l13l7'], bkg_ntuple_dir)

wjets.set_chain_from_dsid_list(g.groups['wjets'], bkg_ntuple_dir)

wgamma.set_chain_from_dsid_list(g.groups['wgamma'], bkg_ntuple_dir)

zgamma.set_chain_from_dsid_list(g.groups['zgamma'], bkg_ntuple_dir)

higgs.set_chain_from_dsid_list(g.groups['higgs'], bkg_ntuple_dir)
#HWW.set_chain_from_dsid_list(g.groups['HWW'], bkg_ntuple_dir)
#HV.set_chain_from_dsid_list(g.groups['HV'], bkg_ntuple_dir)
#ttH.set_chain_from_dsid_list(g.groups['ttH'], bkg_ntuple_dir)
#Htautau.set_chain_from_dsid_list(g.groups['Htautau'], bkg_ntuple_dir)
#ggF_Htt.set_chain_from_dsid_list(g.groups['ggF_Htt'], bkg_ntuple_dir)
#VBF_Htt.set_chain_from_dsid_list(g.groups['VBF_Htt'], bkg_ntuple_dir)

signal_taum.set_chain_from_dsid_list(g.groups['higgs_lfv_taum'], signal_ntuple_dir)
signal_taue.set_chain_from_dsid_list(g.groups['higgs_lfv_taue'], signal_ntuple_dir)
if run_fakes:
    fakes.set_chain_from_dsid_list(g.groups['data'], fake_ntuple_dir)


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
zjets_FF_CR += ' && l_mT[2] < 40'
zjets_FF_CR += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
zjets_FF_CR += ' && MET < 60'
#zjets_FF_CR += ' && fabs(dR_Z_Fake - 3.1415) < 0.5'
#zjets_FF_CR += ' && nBJets == 0'
#zjets_FF_CR += ' && nLJets <= 2'
#zjets_FF_CR += ' && n_jets == 1'

promptlep0 = '(0 < l_truthClass[0] && l_truthClass[0] <= 2)' 
promptlep1 = '(0 < l_truthClass[1] && l_truthClass[1] <= 2)'
promptlep2 = '(0 < l_truthClass[2] && l_truthClass[2] <= 2)'
fakelep2 = '(l_truthClass[2] <= 0 || 2 < l_truthClass[2])'
prompt_truth = " && ".join([promptlep0, promptlep1, promptlep2])
fake_truth = " && ".join([promptlep0, promptlep1, fakelep2])

if run_fakes:
    zjets_FF_CR += " && (nLepID == 2 || !isMC || %s)" % prompt_truth
elif add_truth_num:
    zjets_FF_CR += " && (!isMC || %s)" % prompt_truth
if add_truth_den:
    zjets_FF_CR += " && (!isMC || %s)" % fake_truth


for num_den, id_aid in [("num","ID"),("den",'anti-ID')]:
    REGIONS.append(Region('zjets_FF_CR_'+num_den, 'Z(->ll)+jets','$Z(\\rightarrow\ell\ell)+\\text{jets}$'))
    REGIONS[-1].tcut = zjets_FF_CR
    base_reg = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('e', 'll + %s e'%id_aid, '$\ell\ell$ + %s $e$'%id_aid, 'l_flav[2]==0'))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_e = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('m', 'll + %s m'%id_aid, '$\ell\ell$ + %s $\mu$'%id_aid, 'l_flav[2]==1'))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_m = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('eee', 'ee + %s e'%id_aid, '$ee$ + %s $e$'%id_aid, 'dilep_flav==2 && l_flav[2]==0'))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_eee = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('mme', 'mm + %s e'%id_aid, '$\mu\mu$ + %s $e$'%id_aid, 'dilep_flav==3 && l_flav[2]==0'))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_mme = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('eem', 'ee + %s m'%id_aid, '$ee$ + %s $\mu$'%id_aid, 'dilep_flav==2 && l_flav[2]==1'))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_eem = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('mmm', 'mm + %s m'%id_aid, '$\mu\mu$ + %s $\mu$'%id_aid, 'dilep_flav==3 && l_flav[2]==1'))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_mmm = REGIONS[-1]
    if not (run_fakes or add_truth_den or add_truth_num):
        base_reg_e.compare_regions.append(base_reg_e.build_channel('prompt', 'Prompt', cuts = prompt_truth))
        base_reg_e.compare_regions.append(base_reg_e.build_channel('fake', 'Fake', cuts = fake_truth))
        base_reg_m.compare_regions.append(base_reg_m.build_channel('prompt', 'Prompt', cuts = prompt_truth))
        base_reg_m.compare_regions.append(base_reg_m.build_channel('fake', 'Fake', cuts = fake_truth))
        base_reg_eee.compare_regions.append(base_reg_eee.build_channel('prompt', 'Prompt', cuts = prompt_truth))
        base_reg_eee.compare_regions.append(base_reg_eee.build_channel('fake', 'Fake', cuts = fake_truth))
        base_reg_mme.compare_regions.append(base_reg_mme.build_channel('prompt', 'Prompt', cuts = prompt_truth))
        base_reg_mme.compare_regions.append(base_reg_mme.build_channel('fake', 'Fake', cuts = fake_truth))
        base_reg_eem.compare_regions.append(base_reg_eem.build_channel('prompt', 'Prompt', cuts = prompt_truth))
        base_reg_eem.compare_regions.append(base_reg_eem.build_channel('fake', 'Fake', cuts = fake_truth))
        base_reg_mmm.compare_regions.append(base_reg_mmm.build_channel('prompt', 'Prompt', cuts = prompt_truth))
        base_reg_mmm.compare_regions.append(base_reg_mmm.build_channel('fake', 'Fake', cuts = fake_truth))

#num_den_list = ['num', 'den']
#chan_dict = {'eee' : ['ee','e','Z_dilep_flav==2 && l_flav[2]==0'],
#             'mme' : ['mumu','e','Z_dilep_flav==3 && l_flav[2]==0'],
#             'eem' : ['ee','m','Z_dilep_flav==2 && l_flav[2]==1'],
#             'mmm' : ['mumu','m','Z_dilep_flav==3 && l_flav[2]==1'],
#             'm' : ['ll','m','l_flav[2]==1'],
#             'e' : ['ll','e','l_flav[2]==0'],
#        }
#for num_den in num_den_list:
#    for chan, ops in chan_dict.iteritems():
#        id_or_aid = 'ID' if num_den == 'num' else 'anti-ID'
#        chan_name = '%s + %s %s'%(ops[0], id_or_aid, ops[1])
#
#        name = 'zjets_FF_CR%s_%s'%(num_den, chan)
#        displayname = 'Z+jets FF CR (%s)'%(chan_name)
#        REGIONS.append(Region(name, displayname))
#        REGIONS[-1].tcut = ' && '.join([ops[2], zjets_FF_CR])

# Fake Factor validation region: Wjet
wjets_FF_VR = preselection_sel
wjets_FF_VR += '&& !(l_pt[0] >= 45 && DphiLep1MET < 1.0 && 30 < MLL && MLL < 150)' # orthgonal to SR
wjets_FF_VR += ' && nBJets==0 && l_pt[1] >= 10 && ( !'+mue+' || el1pT_trackclus_ratio < 1.2) && 50 < l_mT[0] && l_mT[1] < 40' # similar to SR
#wjets_FF_VR += ' && 30 < MLL' #TODO: To be added
wjets_FF_VR += ' && (drll <= 2.5 || dpt_ll > 25)' # Orthogonal to Ztautau
wjets_FF_VR += ' && 25 < RelMET && MET < 60' #Remove additional top+VV+Ztt


wjets_emu_FF_VR = '1'
#wjets_emu_FF_VR += ' && dpt_ll > 15'
#wjets_emu_FF_VR += ' && 2.2 < DphiLep0MET'

wjets_mue_FF_VR = '1'

#wjets_FF_VR += ' && nLJets <= 1'
#wjets_FF_VR += ' && n_jets == 0'
#wjets_FF_VR += ' && n_jets == 1'

promptlep0 = '(0 < l_truthClass[0] && l_truthClass[0] <= 2)' 
promptlep1 = '(0 < l_truthClass[1] && l_truthClass[1] <= 2)'
fakelep1 = '(l_truthClass[1] <= 0 || 2 < l_truthClass[1])'
prompt_truth = " && ".join([promptlep0, promptlep1])
fake_truth = " && ".join([promptlep0, fakelep1])

if run_fakes:
    wjets_FF_VR += " && (nLepID == 1 || !isMC || %s)" % prompt_truth
if add_truth_num:
    wjets_FF_VR += " && (!isMC || %s)" % prompt_truth
if add_truth_den:
    wjets_FF_VR += " && (!isMC || %s)" % fake_truth

for num_den, id_aid in zip(["num","den"],['ID','anti-ID']):
    REGIONS.append(Region('wjets_FF_VR_'+num_den, 'W+jets FF VR'))
    REGIONS[-1].tcut = wjets_FF_VR
    base_reg = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('emu', 'ID e + %s m'%id_aid, 'ID $e$ + %s $\mu$'%id_aid, wjets_emu_FF_VR+' && '+emu))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_emu = REGIONS[-1]
    REGIONS.append(base_reg.build_channel('mue', 'ID m + %s e'%id_aid, 'ID $\mu$ + %s $e$'%id_aid, wjets_mue_FF_VR+' && '+mue))
    base_reg.compare_regions.append(REGIONS[-1])
    base_reg_mue = REGIONS[-1]
    #if not (run_fakes or add_truth_den or add_truth_num):
    #    base_reg_emu.compare_regions.append(base_reg_emu.build_channel('prompt', 'Prompt', cuts = prompt_truth))
    #    base_reg_emu.compare_regions.append(base_reg_emu.build_channel('fake', 'Fake', cuts = fake_truth))
    #    base_reg_mue.compare_regions.append(base_reg_mue.build_channel('prompt', 'Prompt', cuts = prompt_truth))
    #    base_reg_mue.compare_regions.append(base_reg_mue.build_channel('fake', 'Fake', cuts = fake_truth))

region_plots = {}
#region_plots['zjets_FF_CRden_eem'] = {
#  'Z_MLL' : deepcopy(plot_defaults['Z_MLL'])
#}
#region_plots['zjets_FF_CRden_eem']['Z_MLL'].update(bin_range=[80, 100], bin_width=1)

#######################################
# What regions to plot
region_ops = []
truth_op = '_fake' if add_truth_den else '_real' if add_truth_num else "" #TMP HACK
if run_den:
    if run_zjets:
        #region_ops += ['zjets_FF_CRden_eee'] # Test Region
        #region_ops += ['zjets_FF_CR_den']
        region_ops += ['zjets_FF_CR_den_m', 'zjets_FF_CR_den_e']
        #region_ops += ['zjets_FF_CR_den_eem', 'zjets_FF_CR_den_mmm']
        #region_ops += ['zjets_FF_CR_den_eee', 'zjets_FF_CR_den_mme']
    if run_base:
        #region_ops += ['wjets_FF_VR_den']
        region_ops += ['wjets_FF_VR_den_emu', 'wjets_FF_VR_den_mue']
        #region_ops += ['ggh_SR_den'+truth_op]
        region_ops += ['ggh_SR_den_emu'+truth_op, 'ggh_SR_den_mue'+truth_op]
elif run_num:
    if run_zjets:
        #region_ops += ['wzCR']
        #region_ops += ['wzCR_e','wzCR_m']
        #region_ops += ['wzCR_eee', 'wzCR_emm','wzCR_mee', 'wzCR_mmm']
        #region_ops += ['zjets_FF_CR_num']
        region_ops += ['zjets_FF_CR_num_m', 'zjets_FF_CR_num_e']
        #region_ops += ['zjets_FF_CR_num_eem', 'zjets_FF_CR_num_mmm']
        #region_ops += ['zjets_FF_CR_num_eee', 'zjets_FF_CR_num_mme']
    if run_base:
        #region_ops += ['wjets_FF_VR_num']
        region_ops += ['wjets_FF_VR_num_emu', 'wjets_FF_VR_num_mue']
        #region_ops += ['preselection_emu', 'preselection_mue']
        #region_ops += ['baseline_emu', 'baseline_mue']
        #region_ops += ['baseline_emu']
        region_ops += ['ggh_SR_num_emu'+truth_op, 'ggh_SR_num_mue'+truth_op] #TMP HACK
        #region_ops += ['preselection']
        #region_ops += ['baseline']
        #region_ops += ['ggh_SR']
        #region_ops += ['ztt_CR']
        #region_ops += ['top_CR']
        #region_ops += ['ztt_CR_emu', 'ztt_CR_mue']
        #region_ops += ['top_CR_emu', 'top_CR_mue']
        #region_ops += ['symmetric_emu','symmetric_mue']
        #region_ops += ['symmetric']
        #region_ops += ['zll_CR']

##For Region Compare Looper
# NOTE: Only activate the sample of interest when making region comparisons
region_tuples = []
#region_tuples.append(('top_CR_emu','ggh_SR_emu'))
#region_tuples.append(('top_CR_mue','ggh_SR_mue'))
#region_tuples.append(('ztt_CR_emu','ggh_SR_emu'))
#region_tuples.append(('ztt_CR_mue','ggh_SR_mue'))
#region_tuples.append(('wjets_FF_VR_num_emu','ggh_SR_emu'))
#region_tuples.append(('wjets_FF_VR_num_mue','ggh_SR_mue'))

#replace region names with region objects
REGION_TUPLES = []
for reg_tuple in region_tuples:
    tup = []
    for reg_name in reg_tuple:
        tup.append(next(x for x in REGIONS if x.name == reg_name))
    REGION_TUPLES.append(tuple(tup))

REGIONS = [x for x in REGIONS if x.name in region_ops]
################################################################################
# Toggle options for execution
# - Options expected to change often between different plottings
################################################################################

#######################################
# Yield Table
from PlotTools.YieldTable import YieldTable, YieldTbl
YieldTbl.write_to_latex = True
YLD_TABLE = YieldTbl()
YLD_TABLE.add_row_formula(name="mc_total",displayname="MC Total", formula="MC")
YLD_TABLE.add_row_formula(name="data_mc_ratio",displayname="Data/MC", formula="data/MC")
YLD_TABLE.add_row_formula(name="signal_mc_ratio",displayname="Signal/MC", formula="SIGNAL/MC")
if not (run_fakes or add_truth_den or add_truth_num):
    YLD_TABLE.add_column_formula(name="other",displayname="Other", formula="BASE_REG - prompt - fake")
if run_zjets and not run_fakes:
    pass
if run_base and not run_fakes:
    #YLD_TABLE.add_row_formula(name="nonwjets_data_ratio",displayname="(MC-WJets)/Data", formula="(MC-wjets)/data")
    pass
if run_fakes:
    YLD_TABLE.add_row_formula(name="fakes_mc_ratio",displayname="Fakes/MC", formula="fakes/MC")
    YLD_TABLE.add_row_formula(name="fakes_data_ratio",displayname="Fakes/Data", formula="fakes/data")


YIELD_TBL = YieldTable()
# Add formulas to yield table
# Key values will be the labels on the table
# Formulas should use sample names and these symbols: +, -, *, /, (, ), [0-9]
# 'MC' is also a recognized label for all backgrounds
YIELD_TBL.formulas['Data/MC'] = "data/MC"


#######################################
# What variables to plot
# Strings for plotting
from example_plot_conf import * #TODO: Dont use import *
PlotBase.save_dir = g.plots_dir
Plot1D.auto_set_ylimits = True
Plot1D.doLogY = False
Plot2D.doLogZ = False
Plot2D.auto_set_zlimits = False

vars_to_plot_all = []
#vars_to_plot_all += ['ptll']
#vars_to_plot_all += ['lep_d0sigBSCorr[0]','lep_z0SinTheta[0]','lep_d0sigBSCorr[1]','lep_z0SinTheta[1]']


################################################################################
# Create plots
################################################################################
PLOTS = []
for region_name in region_ops:
    vars_to_plot = list(vars_to_plot_all)
    region = next(x for x in REGIONS if x.name == region_name)
    region.yield_table = deepcopy(YLD_TABLE)
    if "wjets" in region_name:
        #vars_to_plot += ['el1pT_trackclus_ratio', 'taulep1_pT_ratio']
        #vars_to_plot += ['l_pt[0]', 'l_pt[1]', 'MLL', 'nBJets', 'MET']
        #vars_to_plot += ['RelMET', 'nLJets', 'drll', 'DphiLep0MET','DphiLep1MET', 'l_mT[0]', 'l_mT[1]']
        #vars_to_plot += ['dpt_ll:drll','drll','dpt_ll']
        #vars_to_plot += ['l_pt[0]:drll','l_pt[0]:DphiLep1MET','dpt_ll:drll','dpt_ll:DphiLep1MET']
        #vars_to_plot += ['DphiLep1MET:drll',"DphiLep0MET:drll",'DphiLep1MET:l_mT[1]','DphiLep0MET:l_mT[0]']
        #vars_to_plot += ['drll', 'DphiLep1MET',"DphiLep0MET",'l_mT[1]','l_mT[0]']
        vars_to_plot += ['l_truthClass[1]']

        #vars_to_plot += ['DphiLep0MET:MET', 'DphiLep1MET:MET', 'DphiLep0MET:RelMET', 'DphiLep1MET:RelMET']
        #vars_to_plot += ['DphiLep0MET', 'DphiLep1MET', 'MET', 'RelMET']
        
        region.yield_table.add_row_formula(name="top_mc_ratio", displayname="Top/MC", latexname="Top/MC", formula="top/MC")
        region.yield_table.add_row_formula(name="ztt_mc_ratio", displayname="Ztautau/MC", latexname="$Z\\rightarrow\\tau\\tau$/MC", formula="ztt/MC")
        region.yield_table.add_row_formula(name="vv_mc_ratio",displayname="VV/MC", formula="vv/MC")
        if not run_fakes:
            region.yield_table.add_row_formula(name="wjets_mc_ratio",displayname="W+Jets/MC", formula="wjets/MC")
    
    if "zjets" in region_name:
        #vars_to_plot += ['dR_Z_Fake:(ptll - l_pt[2])', 'ptll:(ptll - l_pt[2])', 'dR_Z_Fake:(ptll - l_pt[2])/ptll', '(ptll - l_pt[2])/ptll', '(ptll - l_pt[2] - MET)/ptll', 'lep_met_pT[2]', 'DphiLep2MET']
        #vars_to_plot += ['l_pt[0]']
        #vars_to_plot += ['l_pt[0]','l_pt[1]', 'l_pt[2]', 'ptll', 'MET', 'MLL', 'l_eta[2]', 'nBJets', 'nLJets']
        #vars_to_plot += ['drl2l[0]','drl2l[1]','drll','DphiLep0MET','DphiLep1MET','DphiLep2MET']
        vars_to_plot += ['l_truthClass[2]']
        region.yield_table.add_row_formula(name="zjets_mc_ratio",displayname="Zll/MC", formula="zll/MC")
        region.yield_table.add_row_formula(name="vv_mc_ratio",displayname="VV/MC", formula="vv/MC")

    if "wzCR" in region_name:
        vars_to_plot += ['MLLL']
        #vars_to_plot += ['l_pt[0]', 'l_pt[1]', 'l_pt[2]', 'l_mT[2]', 'MET', 'MLL']
        #vars_to_plot += ['dR_ZLep0_Fake','dR_ZLep1_Fake','dR_Z_Fake']
        region.yield_table.add_row_formula(name="vv_mc_ratio",displayname="VV/MC", formula="vv/MC")
        region.yield_table.add_row_formula(name="norm_factor",displayname="norm factor", formula="(data-(MC-vv))/vv")
    
    if "symmetric" in region_name:
        #vars_to_plot += ['MCollASym', 'l_pt[0]', 'l_pt[1]','l_eta[0]','l_eta[1]','nLJets','nBJets']
        vars_to_plot += ['MCollASym']

    if 'top_CR' in region_name:
        region.yield_table.add_row_formula(name="top_mc_ratio",displayname="Top/MC", formula="top/MC")
        region.yield_table.add_row_formula(name="norm_factor",displayname="norm factor", formula="(data-(MC-top))/top")

    if 'ztt_CR' in region_name:
        region.yield_table.add_row_formula(name="ztt_mc_ratio", displayname="Ztautau/MC", latexname="$Z\\rightarrow\\tau\\tau$/MC", formula="ztt/MC")
        region.yield_table.add_row_formula(name="norm_factor",displayname="norm factor", formula="(data-(MC-ztt))/ztt")
    
    if any(x in region_name for x in ['baseline','preselection','ggh_SR','ztt_CR','top_CR']):
        vars_to_plot = ['l_truthClass[1]']
        #vars_to_plot += ['taulep1_pT_ratio', 'DphiLep1MET', 'l_mT[0]', 'l_mT[1]']
        #vars_to_plot += ['l_pt[0]', 'l_pt[1]', 'MLL', 'nBJets', 'el1pT_trackclus_ratio', 'MET']
        #vars_to_plot += ['RelMET', 'nLJets', 'drll', 'DphiLep0MET','DphiLep1MET', 'l_mT[0]', 'l_mT[1]', 'taulep1_pT_ratio']
    #if any(x in region_name for x in ['zjets','wzCR']):
    #    vars_to_plot += ['l_pt[0]', 'l_pt[1]', 'l_pt[2]', 'MLL', 'nLJets', 'MET']
    #    vars_to_plot += ['drll', 'DphiLep2MET', 'l_mT[2]', 'dR_ZLep0_Fake','dR_ZLep1_Fake','dR_Z_Fake']
    
    # Remove duplicate names
    vars_to_plot = list(set(vars_to_plot))
    #assert vars_to_plot, ("No plots requested")
    
    for var in vars_to_plot:
        key = var
        # Grab the default plot unless a region specific one is defined
        if region_name in region_plots and key in region_plots[region_name]:
            p = deepcopy(region_plots[region_name][key])
        elif key in plot_defaults:
            p = deepcopy(plot_defaults[key])
        else:
            assert False, ("ERROR :: requested plot not defined:", var)

        # Defing the region and variables
        if p.is2D:
            varx = var.split(':')[1]
            vary = var.split(':')[0]
            p.update(region_name, varx, vary)
        else:
            p.update(region_name, var)

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
