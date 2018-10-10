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
pt_bins = [0,10,11,15,20,25,pt_max]
eta_bins = [0, 1.45, 2.5, 3.0]
njet_bins = [-0.5,0.5,1.5,2.5,10.5]

plot_defaults = {
    'l_pt[2]'            : Plot1D( bin_edges=pt_bins, ptype=Types.stack, doLogY=False, add_overflow = False, xunits='GeV', xlabel='Fake candidate lepton p_{T}'),
    'l_eta[2]'           : Plot1D( bin_edges=eta_bins, ptype=Types.stack, doLogY= False, add_underflow = True, xlabel='Fake candidate lepton #eta'),
    'fabs(l_eta[2]):l_pt[2]'   : Plot2D( xbin_edges=pt_bins, ybin_edges=eta_bins, ylabel='Fake probe lepton |#eta|', xunits='GeV', xlabel='Leading lepton p_{T}'),
    'n_jets:fabs(l_eta[2]):l_pt[2]' : Plot3D ( xbin_edges=pt_bins, ybin_edges=eta_bins, zbin_edges=njet_bins)
}

#plot_defaults['l_pt[2]'].rebin_bins = [0,10,11,15,20,25,100]
#plot_defaults['l_eta[2]'].rebin_bins = [-3.0, -2.5, -1.45, 0, 1.45, 2.5, 3.0]
#plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_xbins = plot_defaults['l_pt[2]'].rebin_bins 
#plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_ybins = [0, 1.45, 2.5, 3.0] 
#plot_defaults['n_jets:fabs(l_eta[2]):l_pt[2]'].rebin_xbins = plot_defaults['l_pt[2]'].rebin_bins
#plot_defaults['n_jets:fabs(l_eta[2]):l_pt[2]'].rebin_ybins = plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_ybins
#plot_defaults['n_jets:fabs(l_eta[2]):l_pt[2]'].rebin_zbins = [-0.5,0.5,1.5,2.5,10.5]
