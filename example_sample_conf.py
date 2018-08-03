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

MCsample.weight_str = 'eventweight'
MCsample.scale_factor = lumi
Sample.input_file_treename = 'superNt'

# Data
data = Data()
SAMPLES.append(data)

################################################################################
# Fakes
## Fakes
fakes = Background("fakes", "Fakes")
fakes.scale_factor = 1 # Derived from data so no need for scaling
fakes.color = ROOT.kGray
SAMPLES.append(fakes)

################################################################################
# Signal
signal_branching_ratio = 0.01
signal_SF = 1
signal_label = "Higgs LFV" if signal_SF == 1 else "Higgs LFV (%dX)"%signal_SF

signal = Signal("higgs_lfv", signal_label)
signal.scale_factor *= signal_branching_ratio * signal_SF
signal.color = ROOT.kGreen
SAMPLES.append(signal)

################################################################################
# Backgrounds

#######################################
# Initialize all backgrounds
# Higgs (Htt + HWW + HV + ttH)
higgs = Background("higgs", "Higgs")
higgs.color = ROOT.kRed - 7
SAMPLES.append(higgs)

# ttbar
ttbar = Background("ttbar", "t#bar{t}")
ttbar.color = ROOT.kOrange+2
SAMPLES.append(ttbar)

# singletop
stop = Background("st", "Single top")
stop.color = ROOT.kOrange+1
SAMPLES.append(stop)

# W+top
wtop = Background("wt", "Wt")
wtop.color = ROOT.kOrange+8
SAMPLES.append(wtop)

# ttbare + lep
ttbarlep = Background("ttbarlep", "ttbarlep")
ttbarlep.color = ROOT.kOrange+5
SAMPLES.append(ttbarlep)

# Top (ttbar + singletop + Wt)
top = Background("top", "Top")
top.color = ROOT.kOrange+1
SAMPLES.append(top)

# VV combined
VV = Background("vv", "Diboson")
VV.color = ROOT.kSpring-8
SAMPLES.append(VV)

# VVV combined
VVV = Background("vvv", "Triboson")
VVV.color = ROOT.kSpring-7
SAMPLES.append(VVV)

# Zee
zee = Background("zee", "Zee")
zee.color = ROOT.kAzure-7
SAMPLES.append(zee)

# Zmumu
zmumu = Background("zmumu", "Zmumu")
zmumu.color = ROOT.kAzure-9
SAMPLES.append(zmumu)

# Ztt
ztt = Background("ztt", "Z#tau#tau")
ztt.color = ROOT.kAzure+8
SAMPLES.append(ztt)

# Wjets
wjets = Background("wjets", "W+jets")
wjets.color = ROOT.kYellow +2
SAMPLES.append(wjets)

# V+gamma
vgamma = Background("vgamma", "V+gamma")
vgamma.color = ROOT.kYellow - 1
SAMPLES.append(vgamma)

# W+gamma
wgamma = Background("wgamma", "W+gamma")
wgamma.color = ROOT.kYellow - 2
SAMPLES.append(wgamma)

# Z+gamma
zgamma = Background("zgamma", "Z+gamma")
zgamma.color = ROOT.kYellow - 3
SAMPLES.append(zgamma)
