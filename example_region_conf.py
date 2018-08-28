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

# Region building blocks
Baseline_Sel = ('l_pt[0] >= 45 && l_pt[1] >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS
              #+ '&& dilep_flav <= 1 '
              + '&&' + lepton_trig_pT)
VBF_stripped = "JetN_g30 >= 2 && j_pt[0] > 40 && Mjj > 400 && DEtaJJ > 3"


########################################
# Create regions
REGIONS = []

REGIONS.append(Region(name = "no_sel", displayname = "No Selections"))
REGIONS[-1].tcut = '1' #DF_OS

REGIONS.append(Region("trig_only", "Lepton Triggers"))
REGIONS[-1].tcut = lepton_trig_pT + "&&" + DF_OS

# Control Regions
REGIONS.append(Region("topCR", "Top CR"))
REGIONS[-1].tcut = "nBJets >= 1 && MET > 40 &&" + DF_OS + " &&" + lepton_trig_pT

REGIONS.append(Region("wzCR", "WZ CR"))
wz_cr_cut = '1'
wz_cr_cut += ' && %s' % singlelep_trig_pT
wz_cr_cut += ' && fabs(lep_d0sigBSCorr[0]) < 15 && fabs(lep_d0sigBSCorr[1]) < 15 && fabs(lep_d0sigBSCorr[2]) < 15'
wz_cr_cut += ' && fabs(lep_z0SinTheta[0]) < 15 && fabs(lep_z0SinTheta[1]) < 15 && fabs(lep_z0SinTheta[2]) < 15'
wz_cr_cut += " && 75 < Z_MLL && Z_MLL < 105"
wz_cr_cut += ' && nBJets == 0'
wz_cr_cut += ' && l_mT[2] > 50'
REGIONS[-1].tcut = wz_cr_cut
REGIONS.append(Region("wzCR_e", "WZ CR (e+Z)"))
REGIONS[-1].tcut = wz_cr_cut + " && " + SF_OS + " && " + lep2_e
REGIONS.append(Region("wzCR_m", "WZ CR (m+Z)"))
REGIONS[-1].tcut = wz_cr_cut + " && " + SF_OS + " && " + lep2_m
REGIONS.append(Region("wzCR_eee", "WZ CR (e+ee)"))
REGIONS[-1].tcut = wz_cr_cut + " && " + Z_ee + " && " + lep2_e
REGIONS.append(Region("wzCR_emm", "WZ CR (e+mm)"))
REGIONS[-1].tcut = wz_cr_cut + " && " + Z_mumu + " && " + lep2_e
REGIONS.append(Region("wzCR_mee", "WZ CR (m+ee)"))
REGIONS[-1].tcut = wz_cr_cut + " && " + Z_ee + " && " + lep2_m
REGIONS.append(Region("wzCR_mmm", "WZ CR (m+mm)"))
REGIONS[-1].tcut = wz_cr_cut + " && " + Z_mumu + " && " + lep2_m

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

ZTauTau_CR =  ('l_pt[0] >= 30 && l_pt[0] < 45 && l_pt[1] >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS)
REGIONS.append(Region("ztautauCR", "Ztautau CR"))
REGIONS[-1].tcut = ZTauTau_CR

ZTauTau_CR =  ('l_pt[0] >= 30 && l_pt[0] < 45 && l_pt[1] >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS)
REGIONS.append(Region("ztautauCR", "Ztautau CR"))
REGIONS[-1].tcut = ZTauTau_CR

# Baseline regions
REGIONS.append(Region("baseline", "Baseline"))
REGIONS[-1].tcut = Baseline_Sel

REGIONS.append(Region("baseline_emu", "Baseline (Channel: El-Mu)"))
REGIONS[-1].tcut = Baseline_Sel + " && " + emu

REGIONS.append(Region("baseline_mue", "Baseline (Channel: Mu-El)"))
REGIONS[-1].tcut = Baseline_Sel + " && " + mue

REGIONS.append(Region("symmetric","Symmetric"))
symm_sel = singlelep_trig
symm_sel += ' && l_pt[0] >= 20 && l_pt[1] >= 20'
symm_sel += ' && fabs(l_eta[0]) <= 2.4 && fabs(l_eta[1]) <= 2.4'
symm_sel += ' && nBJets == 0'
symm_sel += ' && 30 < MLL && MLL < 150'
REGIONS[-1].tcut = symm_sel 
base_sym = REGIONS[-1]
REGIONS.append(base_sym.build_channel('emu', 'e#mu', '$e\mu$', emu))
base_sym.compare_regions.append(REGIONS[-1])
REGIONS.append(base_sym.build_channel('mue', '#mu e', '$\mu e$', mue))
base_sym.compare_regions.append(REGIONS[-1])
# Validation Regions

# Signal Regions
REGIONS.append(Region("vbf", "VBF"))
REGIONS[-1].tcut = "(%s) && (%s)"%(Baseline_Sel, VBF_stripped)

REGIONS.append(Region("optimized", "Optimized"))
REGIONS[-1].tcut = "(%s) && !(%s) && DphiLep1MET < 1 && MtLep0 > 50 && MtLep1 < 40 && ((MET+l_pt[1])/l_pt[1]) > 0.5"%(Baseline_Sel, VBF_stripped)
