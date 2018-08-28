from PlotTools.sample import Sample, MCsample, Data, Background, Signal
import ROOT
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

color_palette = {
    # Preferred ordering for colors when overlaying/stacking multiple plots
    # Colors should all be visually distinguishable and look good :)
    'gray'    : [ROOT.kGray    +0, ROOT.kGray    +1, ROOT.kGray    +2, ROOT.kGray    +3, ROOT.kBlack   +0],
    'red'     : [ROOT.kRed     -7, ROOT.kRed     -3, ROOT.kRed     +2, ROOT.kRed     +3, ROOT.kRed     +4],
    'orange'  : [ROOT.kOrange  +1, ROOT.kOrange  +2, ROOT.kOrange  +8, ROOT.kOrange  +5, ROOT.kOrange  +0],
    'yellow'  : [ROOT.kYellow  +2, ROOT.kYellow  -1, ROOT.kYellow  -2, ROOT.kYellow  -7, ROOT.kYellow  +4],
    'spring'  : [ROOT.kSpring  -8, ROOT.kSpring  -7, ROOT.kSpring  +2, ROOT.kSpring  +3, ROOT.kSpring  +4],
    'green'   : [ROOT.kGreen   -2, ROOT.kGreen   -8, ROOT.kGreen   +2, ROOT.kGreen   +3, ROOT.kGreen   +4],
    'teal'    : [ROOT.kTeal    +0, ROOT.kTeal    +1, ROOT.kTeal    +2, ROOT.kTeal    +3, ROOT.kTeal    +4],
    'cyan'    : [ROOT.kCyan    +0, ROOT.kCyan    +1, ROOT.kCyan    +2, ROOT.kCyan    +3, ROOT.kCyan    +4],
    'azure'   : [ROOT.kAzure   -7, ROOT.kAzure   -9, ROOT.kAzure   +8, ROOT.kAzure   +3, ROOT.kAzure   +4],
    'blue'    : [ROOT.kBlue    -7, ROOT.kBlue    +2, ROOT.kBlue    -8, ROOT.kBlue    -1, ROOT.kBlue    +4],
    'violet'  : [ROOT.kViolet  +0, ROOT.kViolet  +1, ROOT.kViolet  +2, ROOT.kViolet  +3, ROOT.kViolet  +4],
    'magenta' : [ROOT.kMagenta +0, ROOT.kMagenta +1, ROOT.kMagenta +2, ROOT.kMagenta +3, ROOT.kMagenta +4],
    'pink'    : [ROOT.kPink    +0, ROOT.kPink    +1, ROOT.kPink    +2, ROOT.kPink    +3, ROOT.kPink    +4],
}

MCsample.weight_str = 'eventweight'
MCsample.scale_factor = lumi
Sample.input_file_treename = 'superNt'

# Data
data = Data(); SAMPLES.append(data)

################################################################################
# Fakes
## Fakes
fakes = Background("fakes", "Fakes"); SAMPLES.append(fakes)
fakes.scale_factor = 1 # Derived from data so no need for scaling
fakes.color = color_palette['gray'][0] 

################################################################################
# Signal
signal_branching_ratio = 0.01
signal_SF = 10
signal_label = "Higgs LFV" if signal_SF == 1 else "Higgs LFV (%dX)"%signal_SF

signal = Background("higgs_lfv", signal_label); SAMPLES.append(signal)
signal.scale_factor *= signal_branching_ratio * signal_SF
signal.color = color_palette['green'][0] 

signal_taum = Background("higgs_lfv_taum", "H#rightarrow#tau#mu"); SAMPLES.append(signal_taum)
signal_taum.scale_factor *= signal_branching_ratio * signal_SF
signal_taum.color = color_palette['green'][0] 
signal_taue = Background("higgs_lfv_taue", "H#rightarrow#taue"); SAMPLES.append(signal_taue)
signal_taue.scale_factor *= signal_branching_ratio * signal_SF
signal_taue.color = color_palette['green'][1] 
################################################################################
# Backgrounds
#######################################
# Initialize all backgrounds
# Higgs (Htt + HWW + HV + ttH)
higgs = Background("higgs", "Higgs"); SAMPLES.append(higgs)
higgs.color = color_palette['red'][0] 

ggF_Htt = Background("ggF_Htt", "H#tau#tau (ggF)", "$H\\rightarrow\\tau\\tau (ggF)$)"); SAMPLES.append(ggF_Htt)
ggF_Htt.color = color_palette['red'][0] 
VBF_Htt = Background("VBF_Htt", "H#tau#tau (VBF)"); SAMPLES.append(VBF_Htt)
VBF_Htt.color = color_palette['red'][1] 
HWW = Background("HWW", "HWW"); SAMPLES.append(HWW)
HWW.color = color_palette['red'][2] 
HV = Background("HV", "HV"); SAMPLES.append(HV)
HV.color = color_palette['red'][3] 
ttH = Background("ttH", "ttbar + H"); SAMPLES.append(ttH)
ttH.color = color_palette['red'][4] 

### Top (ttbar + singletop + Wt)
top = Background("top", "Top"); SAMPLES.append(top)
top.color = color_palette['orange'][0] 

ttbar = Background("ttbar", "t#bar{t}"); SAMPLES.append(ttbar)
ttbar.color = color_palette['orange'][0] 
stop = Background("st", "Single top"); SAMPLES.append(stop)
stop.color = color_palette['orange'][1] 
wtop = Background("wt", "Wt"); SAMPLES.append(wtop)
wtop.color = color_palette['orange'][2] 
ttbarX = Background("ttbarX", "t#bar{t}+X","$t\\bar{t}+X$"); SAMPLES.append(ttbarX)
ttbarX.color = color_palette['orange'][3] 

# VV combined
VV = Background("vv", "VV"); SAMPLES.append(VV)
VV.color = color_palette['spring'][0]
VV.scale_factor *= 1.05

# VVV combined
VVV = Background("vvv", "VVV"); SAMPLES.append(VVV)
VVV.color = color_palette['spring'][1] 

# Z(->ll) + jets
zll = Background("zll", "Z(->ll)+jets", "$Z(\\rightarrow\ell\ell)+$jets"); SAMPLES.append(zll)
zll.color = color_palette['azure'][0] 

zee = Background("zee", "Zee"); SAMPLES.append(zee)
zee.color = color_palette['azure'][0] 
zmumu = Background("zmumu", "Zmumu"); SAMPLES.append(zmumu)
zmumu.color = color_palette['azure'][1] 

# Ztt
ztt = Background("ztt", "Z#tau#tau","$Z(\\rightarrow\\tau\\tau)+$jets"); SAMPLES.append(ztt)
ztt.color = color_palette['blue'][0] 
    
ztt_sherpa = Background("ztt_sherpa", "Z#tau#tau (Sherpa 2.2.2)"); SAMPLES.append(ztt_sherpa)
ztt_sherpa.color = color_palette['blue'][0]
ztt_l13l7 = Background("ztt_l13l7", "Z#tau#tau (+2 jets)"); SAMPLES.append(ztt_l13l7)
ztt_l13l7.color = color_palette['blue'][0]
ztt_lowMLL = Background("ztt_lowMLL", "Z#tau#tau (Low MLL)"); SAMPLES.append(ztt_lowMLL)
ztt_lowMLL.color = color_palette['blue'][1]
ztt_2jets = Background("ztt_2jets", "Z#tau#tau (+2 jets)"); SAMPLES.append(ztt_2jets)
ztt_2jets.color = color_palette['blue'][2]
    

# Wjets
wjets = Background("wjets", "W+jets"); SAMPLES.append(wjets)
wjets.color = color_palette['yellow'][0]

# V+gamma
vgamma = Background("vgamma", "V+#gamma"); SAMPLES.append(vgamma)
vgamma.color = color_palette['yellow'][1]

# W+gamma
wgamma = Background("wgamma", "W+#gamma","$W+\gamma$"); SAMPLES.append(wgamma)
wgamma.color = color_palette['yellow'][1]

# Z+gamma; 
zgamma = Background("zgamma", "Z+gamma","$Z+\\gamma$"); SAMPLES.append(zgamma)
zgamma.color = color_palette['yellow'][2]
