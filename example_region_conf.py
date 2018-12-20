from PlotTools.region import Region

################################################################################
# Regions
# - Create all the regions that will be applied to samples
################################################################################

########################################
# Define important selections

# General
#TODO: Clean up lepton selections in ntuple maker and here
emu = 'dilep_flav == 0'
mue = 'dilep_flav == 1'
ee = 'dilep_flav == 2'
mumu = 'dilep_flav == 3'
DF_OS = 'dilep_flav <= 1 && ((l_q[0]*l_q[1])<0)'
SF_OS = 'dilep_flav > 1 && ((l_q[0]*l_q[1])<0)'
Z_emu = 'Z_dilep_flav == 0'
Z_mue = 'Z_dilep_flav == 1'
Z_ee = 'Z_dilep_flav == 2'
Z_mumu = 'Z_dilep_flav == 3'
Z_DF_OS = 'Z_dilep_flav <= 1 && ((l_q[0]*l_q[1])<0)'
Z_SF_OS = 'Z_dilep_flav > 1 && ((l_q[0]*l_q[1])<0)'
lep2_e =  'l_flav[2] == 0'
lep2_m =  'l_flav[2] == 1'

# Same flavor lepton triggers
ee15_trig   = 'HLT_2e12_lhloose_L12EM10VH'
mumu15_trig = 'HLT_mu18_mu8noL1'
ee16_trig   = 'HLT_2e17_lhvloose_nod0'
mumu16_trig = 'HLT_mu22_mu8noL1'
# Different flavor lepton triggs
emu15_trig  = 'HLT_e17_lhloose_mu14'
emu16_trig  = 'HLT_e17_lhloose_nod0_mu14'
# Single lepton triggers
e15_trig    = '(HLT_e24_lhmedium_L1EM20VH || HLT_e60_lhmedium || HLT_e120_lhloose)'
mu15_trig   = '(HLT_mu20_iloose_L1MU15 || HLT_mu40)'
e16_trig    = '(HLT_e26_lhtight_nod0_ivarloose || HLT_e60_lhmedium_nod0 || HLT_e140_lhloose_nod0)'
mu16_trig   = '(HLT_mu26_ivarmedium || HLT_mu50)'

# Triggers with added pT requirements
e15_trig_emu_pT   = '(%s && dilep_flav == 0 && l_pt[0] >= 25)'%e15_trig
mu15_trig_emu_pT  = '(%s && dilep_flav == 0 && l_pt[0] < 25 && l_pt[1] >= 21)'%mu15_trig
e16_trig_emu_pT   = '(%s && dilep_flav == 0 && l_pt[0] >= 27)'%e16_trig
mu16_trig_emu_pT  = '(%s && dilep_flav == 0 && l_pt[0] < 27 && l_pt[1] >= 28)'%mu16_trig
emu15_trig_emu_pT = '(%s && dilep_flav == 0 && 18 <= l_pt[0] && l_pt[0] < 25 && 15 <= l_pt[1] && l_pt[1] < 21)'%emu15_trig
emu16_trig_emu_pT = '(%s && dilep_flav == 0 && 18 <= l_pt[0] && l_pt[0] < 27 && 15 <= l_pt[1] && l_pt[1] < 28)'%emu16_trig
e15_trig_mue_pT   = '(%s && dilep_flav == 1 && l_pt[1] >= 25)'%e15_trig
mu15_trig_mue_pT  = '(%s && dilep_flav == 1 && l_pt[1] < 25 && l_pt[0] >= 21)'%mu15_trig
e16_trig_mue_pT   = '(%s && dilep_flav == 1 && l_pt[1] >= 27)'%e16_trig
mu16_trig_mue_pT  = '(%s && dilep_flav == 1 && l_pt[1] < 27 && l_pt[0] >= 28)'%mu16_trig
emu15_trig_mue_pT = '(%s && dilep_flav == 1 && 18 <= l_pt[1] && l_pt[1] < 25 && 15 <= l_pt[0] && l_pt[0] < 21)'%emu15_trig
emu16_trig_mue_pT = '(%s && dilep_flav == 1 && 18 <= l_pt[1] && l_pt[1] < 27 && 15 <= l_pt[0] && l_pt[0] < 28)'%emu16_trig
e15_trig_ee_pT   = '(%s && Z_dilep_flav == 2 && l_pt[0] >= 25)'%e15_trig
e16_trig_ee_pT   = '(%s && Z_dilep_flav == 2 && l_pt[0] >= 27)'%e16_trig
mu15_trig_mumu_pT  = '(%s && Z_dilep_flav == 3 && l_pt[0] >= 21)'%mu15_trig
mu16_trig_mumu_pT  = '(%s && Z_dilep_flav == 3 && l_pt[0] >= 27)'%mu16_trig

# Combined triggers
e15_trig_pT = '(%s || %s || %s)'%(e15_trig_emu_pT, e15_trig_mue_pT, e15_trig_ee_pT)
mu15_trig_pT = '(%s || %s || %s)'%(mu15_trig_emu_pT, mu15_trig_mue_pT,mu15_trig_mumu_pT)
e16_trig_pT = '(%s || %s || %s)'%(e16_trig_emu_pT, e16_trig_mue_pT, e16_trig_ee_pT)
mu16_trig_pT = '(%s || %s || %s)'%(mu16_trig_emu_pT, mu16_trig_mue_pT, mu16_trig_mumu_pT)
emu15_trig_pT = '(%s || %s)'%(emu15_trig_emu_pT, emu15_trig_mue_pT)
emu16_trig_pT = '(%s || %s)'%(emu16_trig_emu_pT, emu16_trig_mue_pT)

dilep15_trig      = '(treatAsYear==2015 && %s)'%emu15_trig
dilep16_trig     = '(treatAsYear==2016 && %s)'%emu16_trig
singlelep15_trig = '(treatAsYear==2015 && (%s || %s))'%(e15_trig,mu15_trig)
singlelep16_trig = '(treatAsYear==2016 && (%s || %s))'%(e16_trig,mu16_trig)
dilep_trig = '(%s || %s)'%(dilep15_trig, dilep16_trig)
singlelep_trig = '(%s || %s)'%(singlelep15_trig, singlelep16_trig)
lepton_trig = '(%s || %s)'%(singlelep_trig, dilep_trig)

dilep15_trig_pT      = '(treatAsYear==2015 && %s)'%emu15_trig_pT
dilep16_trig_pT     = '(treatAsYear==2016 && %s)'%emu16_trig_pT
singlelep15_trig_pT = '(treatAsYear==2015 && (%s || %s))'%(e15_trig_pT,mu15_trig_pT)
singlelep16_trig_pT = '(treatAsYear==2016 && (%s || %s))'%(e16_trig_pT,mu16_trig_pT)
dilep_trig_pT = '(%s || %s)'%(dilep15_trig_pT, dilep16_trig_pT)
singlelep_trig_pT = '(%s || %s)'%(singlelep15_trig_pT, singlelep16_trig_pT)
lepton_trig_pT = '(%s || %s)'%(singlelep_trig_pT, dilep_trig_pT)

################################################################################
# Create regions
################################################################################
REGIONS = []

REGIONS.append(Region(name = "no_sel", displayname = "No Selections"))
REGIONS[-1].tcut = '1' #DF_OS

REGIONS.append(Region("trig_only", "Lepton Triggers"))
REGIONS[-1].tcut = lepton_trig_pT + "&&" + DF_OS


################################################################################
# Preselection regions

REGIONS.append(Region("preselection", "Preselection"))
preselection_sel = singlelep_trig_pT
#preselection_sel += ' && (fabs(lep_d0sigBSCorr[0]) < 15 && fabs(lep_d0sigBSCorr[1]) < 15)'
#preselection_sel += ' && (fabs(lep_z0SinTheta[0]) < 15 && fabs(lep_z0SinTheta[1]) < 15)'

REGIONS[-1].tcut = preselection_sel
preselection = REGIONS[-1]
REGIONS.append(preselection.build_channel('emu', 'e#mu', '$e\mu$', emu))
preselection.compare_regions.append(REGIONS[-1])
REGIONS.append(preselection.build_channel('mue', '#mue', '$\mu e$', mue))
preselection.compare_regions.append(REGIONS[-1])

# Baseline regions
REGIONS.append(Region("baseline", "Baseline"))
baseline_sel = preselection_sel
baseline_sel += ' && ' + emu
baseline_sel += ' && l_pt[0] >= 45'
baseline_sel += ' && l_pt[1] >= 15'
baseline_sel += ' && (30 < MLL && MLL < 150)'
baseline_sel += ' && nBJets == 0'
baseline_sel += ' && ( !('+mue+') || el1pT_trackclus_ratio < 1.2)'

REGIONS[-1].tcut = baseline_sel
baseline = REGIONS[-1]
REGIONS.append(baseline.build_channel('emu', 'e#mu', '$e\mu$', emu))
baseline.compare_regions.append(REGIONS[-1])
REGIONS.append(baseline.build_channel('mue', '#mue', '$\mu e$', mue))
baseline.compare_regions.append(REGIONS[-1])

REGIONS.append(Region("symmetric","Symmetric"))
symm_sel = preselection_sel
symm_sel += ' && l_pt[0] >= 20 && l_pt[1] >= 20'
symm_sel += ' && fabs(l_eta[0]) <= 2.4 && fabs(l_eta[1]) <= 2.4'
symm_sel += ' && nBJets == 0'
symm_sel += ' && 30 < MLL && MLL < 150'
REGIONS[-1].tcut = symm_sel 
base_sym = REGIONS[-1]
REGIONS.append(base_sym.build_channel('emu', 'e#mu', '$e\mu$', emu))
base_sym.compare_regions.append(REGIONS[-1])
REGIONS.append(base_sym.build_channel('mue', '#mue', '$\mu e$', mue))
base_sym.compare_regions.append(REGIONS[-1])
# Validation Regions

################################################################################
# Signal Regions
VBF_stripped = "JetN_g30 >= 2 && j_pt[0] > 40 && Mjj > 400 && DEtaJJ > 3"
REGIONS.append(Region("vbf", "VBF"))
REGIONS[-1].tcut = "(%s) && (%s)"%(baseline_sel, VBF_stripped)

ggH_sel = baseline_sel
ggH_sel += ' && DphiLep1MET < 1'
ggH_sel += ' && l_mT[0] > 50'
ggH_sel += ' && l_mT[1] < 40'
ggH_sel += ' && taulep1_pT_ratio > 0.5'

promptlep0 = '(0 < l_truthClass[0] && l_truthClass[0] <= 2)' 
promptlep1 = '(0 < l_truthClass[1] && l_truthClass[1] <= 2)'
fakelep1 = '(l_truthClass[1] <= 0 || 2 < l_truthClass[1])'
prompt_truth = " && ".join([promptlep0, promptlep1])
fake_truth = " && ".join([promptlep0, fakelep1])
den_fake_sel = "(!isMC || %s)" % fake_truth
den_prompt_sel = "(!isMC || %s)" % prompt_truth
num_fake_sel = "(nLepID == 1 || !isMC || %s)" % fake_truth
num_prompt_sel = "(nLepID == 1 || !isMC || %s)" % prompt_truth

for num_den, id_aid, fake_sel, prompt_sel in zip(["num","den"],['ID','anti-ID'],[num_fake_sel, den_fake_sel],[num_prompt_sel, den_prompt_sel]):
    REGIONS.append(Region('ggh_SR_'+num_den, 'gg#rightarrowH SR','$gg\\rightarrow H SR$'))
    base_reg = REGIONS[-1]
    base_reg.tcut = ggH_sel
    base_reg.isSR = True

    REGIONS.append(base_reg.build_channel('emu', 'ID e + %s m'%id_aid, 'ID $e$ + %s $\mu$'%id_aid, emu))
    base_reg_emu = REGIONS[-1]
    base_reg.compare_regions.append(base_reg_emu)
    REGIONS.append(base_reg_emu.build_channel('fake', 'Fake', 'Fake', fake_sel))
    REGIONS.append(base_reg_emu.build_channel('real', 'Real', 'Real', prompt_sel))

    REGIONS.append(base_reg.build_channel('mue', 'ID m + %s e'%id_aid, 'ID $\mu$ + %s $e$'%id_aid, mue))
    base_reg_mue = REGIONS[-1]
    base_reg.compare_regions.append(base_reg_mue)
    REGIONS.append(base_reg_mue.build_channel('fake', 'Fake', 'Fake', fake_sel))
    REGIONS.append(base_reg_mue.build_channel('real', 'Real', 'Real', prompt_sel))

################################################################################
# Control Regions
################################################################################

# ttbar 
top_sel = ggH_sel.replace("nBJets == 0","nBJets > 0")
REGIONS.append(Region("top_CR", "Top CR"))
REGIONS[-1].tcut = top_sel
base = REGIONS[-1]
REGIONS.append(base.build_channel('emu', 'e#mu', '$e\mu$', emu))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('mue', '#mue', '$\mu e$', mue))
base.compare_regions.append(REGIONS[-1])

# Z (->tautau) + jets
ztt_sel = ggH_sel.replace("l_pt[0] >= 45","l_pt[0] < 45")
ztt_sel += ' && (drll > 2.5 && dpt_ll < 25)'

REGIONS.append(Region("ztt_CR", "Z->tautau CR","$Z(\\rightarrow\\tau\\tau)$ + jets"))
REGIONS[-1].tcut = ztt_sel
base = REGIONS[-1]
REGIONS.append(base.build_channel('emu', 'e#mu', '$e\mu$', emu))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('mue', '#mue', '$\mu e$', mue))
base.compare_regions.append(REGIONS[-1])

#W + jets (moved to main conf)
#wjets_sel = ztt_sel.replace('(drll > 2.5 && dpt_ll < 25)','!(drll > 2.5 && dpt_ll < 25)')
#REGIONS.append(Region("wjets_CR", "W+jets CR","$W$ + jets"))
#REGIONS[-1].tcut = wjets_sel
#base = REGIONS[-1]
#REGIONS.append(base.build_channel('emu', 'e#mu', '$e\mu$', emu))
#base.compare_regions.append(REGIONS[-1])
#REGIONS.append(base.build_channel('mue', '#mue', '$\mu e$', mue))
#base.compare_regions.append(REGIONS[-1])

# WZ
REGIONS.append(Region("wzCR", "WZ CR"))
wz_cr_cut = preselection_sel
wz_cr_cut += " && 75 < Z_MLL && Z_MLL < 105"
wz_cr_cut += ' && nBJets == 0'
wz_cr_cut += ' && l_mT[2] > 50'
REGIONS[-1].tcut = wz_cr_cut
base = REGIONS[-1]
REGIONS.append(base.build_channel('e','e+ll','$e+\ell\ell$', SF_OS + " && " + lep2_e))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('m','m+ll','$m+\ell\ell$', SF_OS + " && " + lep2_m))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('eee','e+ee','$e+ee$', Z_ee + " && " + lep2_e))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('emm','e+mm','$e+\mu\mu $', Z_mumu + " && " + lep2_e))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('mee','m+ee','$\mu +ee$', Z_ee + " && " + lep2_m))
base.compare_regions.append(REGIONS[-1])
REGIONS.append(base.build_channel('mmm','m+mm','$\mu +\mu\mu $', Z_mumu + " && " + lep2_m))
base.compare_regions.append(REGIONS[-1])

# Z(->ll) + jets
zll_cr_base = "1"
#zll_cr_base = "75 < Z_MLL && Z_MLL < 105 && " + Z_SF_OS
zll_cr_add  = '1'
#zll_cr_add  += " && %s"%singlelep_trig_pT
#zll_cr_add += " && nLepID == 2"
#zll_cr_add += " && nLepAntiID >= 1"
#zll_cr_add += " && nLepID == 3"
#zll_cr_add += " && l_flav[2]==0" # electron
zll_cr_add += " && l_flav[2]==1" # muon

REGIONS.append(Region("zCR", "Z CR"))
REGIONS[-1].tcut = zll_cr_base + " && " + zll_cr_add

REGIONS.append(Region("zCR_ee", "Z CR (Channel: El-El)"))
REGIONS[-1].tcut = zll_cr_base + " && " + zll_cr_add + " && " + Z_ee

REGIONS.append(Region("zCR_mumu", "Z CR (Channel: Mu-Mu)"))
REGIONS[-1].tcut = zll_cr_base + " && " + zll_cr_add + " && " + Z_mumu

