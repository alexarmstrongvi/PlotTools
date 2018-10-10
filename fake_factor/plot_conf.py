from PlotTools.plot import PlotBase, Plot1D, Plot2D, Plot3D, Types
from copy import deepcopy

################################################################################
# Variables
# - Create all the variables that one might like to plot
################################################################################

PlotBase.output_format = 'pdf'
Plot1D.auto_set_ylimits = False
Plot1D.doLogY = False
Plot1D.ymax = 4e3

# Improved plot defining setup
# These plots will be copied into new plots for each region being run
for_root_file = True
pt_max = 1000.0 if for_root_file else 100.0

plot_defaults = {
    'l_pt[2]'            : Plot1D( bin_range=[0.0, pt_max],  bin_width=0.5, ptype=Types.stack, doLogY=False, add_overflow = False, xunits='GeV', xlabel='Fake candidate lepton p_{T}'),
    'l_eta[2]'           : Plot1D( bin_range=[-3.0, 3.0, 0, 4000],   bin_width=0.01, ptype=Types.stack, doLogY= False, add_underflow = True, xlabel='Fake candidate lepton #eta'),
    'l_eta[2]'           : Plot1D( bin_range=[-3.0, 3.0, 0, 4000],   bin_width=0.01, ptype=Types.stack, doLogY= False, add_underflow = True, xlabel='Fake candidate lepton #eta'),
    'fabs(l_eta[2]):l_pt[2]'   : Plot2D( bin_range=[0, pt_max, 0, 3.0], xbin_width = 1, ybin_width = 0.01, ylabel='Fake probe lepton |#eta|', xunits='GeV', xlabel='Leading lepton p_{T}'),
    'n_jets:fabs(l_eta[2]):l_pt[2]' : Plot3D ( bin_range=[0, pt_max, 0, 3.0, -0.5, 10.5], xbin_width=1, ybin_width=0.01, zbin_width=1)
}

plot_defaults['l_pt[2]'].rebin_bins = [0,10,11,15,20,25,1000]
plot_defaults['l_eta[2]'].rebin_bins = [-3.0, -2.5, -1.45, 0, 1.45, 2.5, 3.0]
plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_xbins = plot_defaults['l_pt[2]'].rebin_bins 
plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_ybins = [0, 1.45, 2.5, 3.0] 
plot_defaults['n_jets:fabs(l_eta[2]):l_pt[2]'].rebin_xbins = plot_defaults['l_pt[2]'].rebin_bins
plot_defaults['n_jets:fabs(l_eta[2]):l_pt[2]'].rebin_ybins = plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_ybins
plot_defaults['n_jets:fabs(l_eta[2]):l_pt[2]'].rebin_zbins = [-0.5,0.5,1.5,2.5,10.5]
