"""
================================================================================
Plot classes for ROOT plotting

Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with almost the entire code borrowed exactly from Danny Antrim
        <dantrim@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""
# General python
import os, sys
import glob
import re
from math import floor
from enum import Enum
from collections import OrderedDict
from array import array

# Root data analysis framework
import ROOT as r
r.TCanvas.__init__._creates = False

import PlotTools.plot_utils as pu
#import AtlasStyle #TODO: Can't load it correctly
#SetAtlasStyle()

################################################################################
# Enumerating plot types
################################################################################
class Types(Enum):
    default = 0
    stack = 1
    ratio = 2
    double_ratio = 3
    comparison = 4
    two_dim = 5
    three_dim = 6
    profile = 7
    undefined = 8
################################################################################
# Plot Classes
################################################################################
# TODO add some generic functions to plot utils

class PlotBase(object):
    save_dir = './'
    output_format = 'pdf'
    def __init__(self):
        self.suffix = ""
    def save_plot(self, can, suffix = ""):
        if suffix:
            suffix = "_" + suffix
        if self.suffix:
            suffix += "_" + self.suffix
        save_name = self.name + suffix + "." + self.output_format
        save_path = os.path.join(self.save_dir, save_name)
        save_path = os.path.normpath(save_path)
        can.SaveAs(save_path)

class Plot1D(PlotBase) :
    # Class defaults for axis range depending on log and normalization settings
    xmin = 0
    xmax = 50.0
    ymin = 0
    ymax =1e4

    logy_min = 1e-1
    logy_max = 1e10

    norm_ymin = 0
    norm_ymax = 1.5

    logy_norm_min = 1e-7
    logy_norm_max = 1e4

    auto_set_ylimits = True
    doLogY = True


    def __init__(self,
        region = "",
        name = "",
        variable = "",
        xlabel = "x-Label",
        xunits = "",
        ylabel = "",
        yunits = "",
        suffix = "",
        bin_range = [], # [x0, x1, y0, y1]
        bin_width = None, # Can specify to override nbins
        nbins = 20,
        bin_edges = [],
        doLogY = None,
        doNorm = False,
        add_overflow = True,
        add_underflow = False,
        leg_is_left = False,
        leg_is_bottom_right = False,
        leg_is_bottom_left = False,
        xcut_is_max = True, #toggle for CutScan1D hists
        ptype = Types.default,
        ) :
        '''
        Constructor

        Order is important so be careful when rearranging.
        '''
        # Descriptors
        self.region = region
        self.variable = variable
        self.name = name if name else determine_name(region, variable)
        self.suffix = suffix


        # Flags
        self.is2D = False
        self.is3D = False
        if doLogY: self.doLogY = doLogY
        self.doNorm = doNorm

        self.add_overflow = add_overflow
        self.add_underflow = add_underflow
        self.leg_is_left = leg_is_left
        self.leg_is_bottom_right = leg_is_bottom_right
        self.leg_is_bottom_left = leg_is_bottom_left
        self.xcut_is_max = xcut_is_max
        self.ptype = ptype
        
        # Objects
        self.pads = None

        assert not bool(bin_edges) or not bool(bin_width), (
                "ERROR Plot1D :: Only provide one bin value option")
        if bin_edges:
            bin_range = [bin_edges[0], bin_edges[-1]]
        self.xmin, self.xmax, self.ymin, self.ymax = self.determine_range(bin_range)
        self.nbins, self.bin_edges = determine_bins(bin_edges, bin_width, nbins, self.xmin, self.xmax, variable)

        self.xunits = xunits
        self.yunits = yunits
        self.xlabel, self.ylabel = self.determine_labels(xlabel, ylabel)

        self.bin_labels = []
        self.rebin_bins = []

    def update(self,
        region = None,
        variable = None,
        name = None,
        bin_range = None,
        bin_width = None,
        bin_edges = None,
        nbins = None,
        xmin = None,
        xmax = None,
        ymin = None,
        ymax = None,
        logy_min = None,
        logy_max = None,
        norm_ymin = None,
        norm_ymax = None,
        logy_norm_min = None,
        logy_norm_max = None,
        doLogY = None,
        doNorm = None,
        add_overflow = None,
        add_underflow = None):
        '''
        Update plot properties

        This is mainly for updating plot properties of a cloned plot. Therefore,
        only a minimum set of properties that one expects could be updated are
        included.
        '''
        if region: self.region = region
        if variable: self.variable = variable
        if name:
            self.name = name
        elif not name and (region or variable):
            self.name = determine_name(region, variable)


        # Flags
        if doLogY: self.doLogY = doLogY
        if doNorm: self.doNorm = doNorm
        if add_overflow: self.add_overflow = add_overflow
        if add_underflow: self.add_underflow = add_underflow

        # Properties
        if xmin: self.xmin = xmin
        if xmax: self.xmax = xmax
        if ymin: self.ymin = ymin
        if ymax: self.ymax = ymax
        if logy_min: self.logy_min = logy_min
        if logy_max: self.logy_max = logy_max
        if norm_ymin: self.norm_ymin = norm_ymin
        if norm_ymax: self.norm_ymax = norm_ymax
        if logy_norm_min: self.logy_norm_min = logy_norm_min
        if logy_norm_max: self.logy_norm_max = logy_norm_max

        self.xmin, self.xmax, self.ymin, self.ymax = self.determine_range(bin_range)
        if nbins or bin_width:
            self.nbins = determine_nbins(self.xmax, self.xmin, bin_width, self.variable) if bin_width else nbins
        
        assert bool(bin_edges) + bool(nbins) + bool(bin_width) <= 1, (
                "ERROR Plot1D :: Only provide one bin value option")
        if bin_edges:
            bin_range = [bin_edges[0], bin_edges[-1]]
        self.xmin, self.xmax, self.ymin, self.ymax = self.determine_range(bin_range)
        if bin_range or nbins or bin_width:
            self.nbins, self.bin_edges = determine_bins(bin_edges, bin_width, nbins, self.xmin, self.xmax, variable)

    def bin_width(self):
        return (self.xmax - self.xmin) / self.nbins

    def determine_range(self, bin_range):
        '''
        Inteligently determine bin range

        The user can provide 2 or 4 values in a list. The default range values
        will be modified. It is assumed providing two values indicates one only
        wants to set the x-axis range.

        args:
            bin_range (list(int or float)) - list of bin range values

        returns:
            tuple (4 int or float) - range values for both x- and y-axis

        '''
        # default values

        assert not bin_range or len(bin_range) in [2,4],(
            'ERROR :: Unrecognized bin range format:', bin_range)

        if self.doLogY and self.doNorm:
            ymin = self.logy_norm_min
            ymax = self.logy_norm_max
        elif self.doLogY:
            ymin = self.logy_min
            ymax = self.logy_max
        elif self.doNorm:
            ymin = self.norm_ymin
            ymax = self.norm_ymax
        else:
            ymin = self.ymin
            ymax = self.ymax

        xmin = self.xmin
        xmax = self.xmax

        if not bin_range:
            return  xmin, xmax, ymin, ymax
        elif len(bin_range) == 2:
            return bin_range[0], bin_range[1], ymin, ymax
        elif len(bin_range) == 4:
            return bin_range[0], bin_range[1], bin_range[2], bin_range[3]

    def determine_labels(self, xlabel, ylabel):

        # Set x-axis label
        if self.xunits:
            xlabel = "%s [%s]"%(xlabel, self.xunits)

        # set y-axis label defaults
        if ylabel:
            pass
        elif self.doNorm:
            ylabel = "a.u."
        else:
            ylabel = "Events"

        # set y-axis label
        width_label = str(round(self.bin_width(), 2))
        if not self.xunits and width_label == '1.0':
            pass
        elif self.xunits and width_label == '1.0':
            ylabel = "%s / %s"%(ylabel, self.xunits)
        else:
            ylabel = "%s / %s %s"%(ylabel,width_label,self.xunits)

        return xlabel, ylabel

    def set_bin_labels(self, axis_hist):
        '''
        Add arbitrary labels to TH1 bins

        Set labels for a TH1 histogram after defining bin_labels. It is
        required that the number of labels be the same as the number of bins. If
        one wants only certain bins to be labeled, then put empty strings for
        those indices and they will be skipped. Put a blank string to create
        empty labels

        args:
            axis_hist (TH1): histogram that includes plot axis

        '''
        n_labels = len(self.bin_labels)
        n_bins = axis_hist.GetNbinsX()
        assert n_labels == n_bins, (
            "ERROR :: nbins (%d) != nlabels (%d)"%(n_bins, n_labels))

        x_axis = axis_hist.GetXaxis()
        for ibin, label in zip(range(n_labels), self.bin_labels):
            if not label: continue
            x_axis.SetBinLabel(ibin+1, label)
        x_axis.SetLabelOffset(0.005)
        x_axis.SetLabelSize(1.5*x_axis.GetLabelSize())

    def setDefaultPads(self, name) :
        self.pads = Pads(name)
        self.ptype = Types.default

    def setStackPads(self, name):
        self.pads = StackPads(name)
        self.ptype = Types.stack

    def setRatioPads(self, name) :
        self.pads = RatioPads(name)
        self.ptype = Types.ratio

    def setDoubleRatioPads(self, name) :
        self.pads = DoubleRatioPads(name)
        self.ptype = Types.double_ratio

    def setCutScanPads1D(self, name):
        self.pads = CutScanPads1D(name)
        #self.ptype = Types.cut_scan_1d #Do I need this?

    def make_data_mc_stack_plot(self, reg_name, hists, suffix=''):
        # Get Canvas
        can = self.pads.canvas
        can.cd()
        if self.doLogY : can.SetLogy(True)

        # Checks - Move on to next plot in case of failure
        if not hists.mc_stack:
            print "WARNING :: Stack plot has either no MC. Skipping."
            return

        self.draw_data_mc_stack_plot(reg_name, hists)

        # Reset axis
        can.RedrawAxis()
        can.SetTickx()
        can.SetTicky()
        can.Update()

        # Save the histogram

        self.save_plot(can, suffix)

        # Clean up
        #plot.pads.canvas.Clear() #TODO: Figure out when to clear canvas (do I need to)
    def draw_data_mc_stack_plot(self, reg_name, hists):
        ''' In a separate function so it can be used when making stack + ratio plots '''
        # Draw the histograms
        hists.axis.Draw()
        hists.mc_stack.Draw("HIST SAME")
        hists.mc_errors.Draw("E2 same")
        hists.mc_total.Draw('hist same')
        for hsig in hists.signals: hsig.Draw("hist same")
        if hists.data: hists.data.Draw("option same pz 0")
        hists.leg.Draw()
        hists.leg_sig.Draw()
        pu.draw_atlas_label('Internal','Higgs LFV', reg_name)

    def make_data_mc_stack_with_ratio_plot(self, reg_name, stack_hists, ratio_hists):
        # Pads
        rcan = self.pads #remove relableing

        ############################################################################
        # Top stack plot
        rcan.upper_pad.cd()
        if self.doLogY : rcan.upper_pad.SetLogy(True)
        rcan.upper_pad.Update()


        # Draw the histograms
        self.draw_data_mc_stack_plot(reg_name, stack_hists)

        # Reset axis
        rcan.upper_pad.RedrawAxis()
        rcan.upper_pad.SetTicks()
        rcan.upper_pad.Update()

        ############################################################################
        # Bottom ratio plot
        rcan.lower_pad.cd()

        ratio_hists.axis.Draw("AXIS")
        ratio_hists.errors.Draw("E2")
        ratio_hists.ratio.Draw("option same pz 0")

        pu.draw_line(self.xmin, 1.5, self.xmax, 1.5, style = 3, width = 1)
        pu.draw_line(self.xmin, 1.0, self.xmax, 1.0, style = 2, width = 1, color = r.kBlack)
        pu.draw_line(self.xmin, 0.5, self.xmax, 0.5, style = 3, width = 1)

        rcan.lower_pad.SetTicks()
        rcan.lower_pad.Update()

        ############################################################################
        # Save the histogram
        self.save_plot(rcan.canvas)

        # Clean up
        #self.pads.canvas.Clear() #TODO: check if this can be uncommented

    def make_region_compare_plot(self, sample_name, hists, suffix=""):
        # Get Canvas
        can = self.pads.canvas
        can.cd()
        if self.doLogY : can.SetLogy(True)

        hists.axis.Draw()
        for hist in hists.hists:
            hist.SetLineWidth(3)
            hist.Draw("HIST SAME")
        hists.leg.Draw()
        pu.draw_atlas_label('Internal','Higgs LFV', sample_name)

        # Reset axis
        can.RedrawAxis()
        can.SetTickx()
        can.SetTicky()
        can.Update()

        # Save the histogram

        self.save_plot(can, suffix)

    def make_cutscan1d_plot(self, hists, reg_name, roc_curve=False):
        # Get Canvas
        pads = self.pads
        pads.hist_pad.cd()
        if self.doLogY : pads.hist_pad.SetLogy(True)

        # Make signal background overlay
        hists.axis.Draw()
        hists.signal.SetFillColor(0)
        hists.signal.SetLineWidth(3)
        hists.signal.Draw("HIST SAME")
        hists.bkgd.SetFillColor(0)
        hists.bkgd.SetLineWidth(3)
        hists.bkgd.Draw("HIST SAME")
        hists.hist_leg.Draw()
        pu.draw_atlas_label('Internal','Higgs LFV', reg_name)
        
        # Reset axis
        pads.hist_pad.RedrawAxis()
        pads.hist_pad.SetTickx()
        pads.hist_pad.SetTicky()
        pads.hist_pad.Update()

        # Add Eff/Rej plot
        pads.eff_pad.cd()
        hists.signal_eff.SetFillColor(0)
        hists.signal_eff.SetLineWidth(3)
        hists.signal_eff.SetLineWidth(3)
        hists.signal_eff.GetYaxis().SetLabelOffset(0.85*hists.signal_eff.GetYaxis().GetLabelOffset())
        hists.signal_eff.Draw("SAME")
        hists.bkgd_rej.SetFillColor(0)
        hists.bkgd_rej.SetLineWidth(3)
        hists.bkgd_rej.Draw("SAME")
        hists.eff_leg.Draw()

        # Add S/B plot
        pads.s_over_b_pad.cd()
        hists.s_over_b.SetFillColor(hists.s_over_b.GetLineColor())
        hists.s_over_b.GetYaxis().SetLabelOffset()
        hists.s_over_b.Draw("SAME")

        # Add ROC Curve
        pads.roc_pad.cd()
        #pads.roc_pad.SetLogy(True)
        hists.roc_graph.Draw("AC*")

        # Save the histogram
        self.save_plot(pads.canvas, suffix="cutScan")


    def Print(self) :
        print "Plot1D    plot: %s  (region: %s  var: %s)"%(
            self.name, self.region, self.variable)
 
class Plot2D(PlotBase) :
    xmin = 0
    xmax = 50.0
    ymin = 0
    ymax = 50.0
    zmin = 0
    zmax = 1e4

    logz_min = 1e-1
    logz_max = 1e10

    norm_zmin = 0
    norm_zmax = 1.5

    logz_norm_min = 1e-7
    logz_norm_max = 1e4

    auto_set_zlimits = True
    doLogZ = True
    def __init__(self,
        region = "",
        name = "",
        xvariable = "",
        yvariable = "",
        xlabel = "x-Label",
        xunits = "",
        ylabel = "y-Label",
        yunits = "",
        zlabel = "",
        zunits = "",
        suffix = "",
        bin_range = [], # [x0, x1, y0, y1, z0, z1]
        ybin_width = None, # Can specify to override nbins
        xbin_width = None, # Can specify to override nbins
        nxbins = 20,
        nybins = 20,
        xbin_edges = [],
        ybin_edges = [],
        doLogZ = None,
        doNorm = False,
        add_overflow = False,
        add_underflow = False,
        ptype = Types.two_dim,
        style = 'colz'
        ) :
        # Descriptors
        self.region = region
        self.xvariable = xvariable
        self.yvariable = yvariable
        self.name = name if name else determine_name(region, xvariable, yvariable)
        self.suffix = suffix

        # Objects
        self.pads = None

        #Flags
        self.is2D = True
        self.is3D = False
        if doLogZ: self.doLogZ = doLogZ
        self.doNorm = doNorm
        self.add_overflow = add_overflow
        self.add_underflow = add_underflow
        self.ptype = ptype

        # Properties
        if xbin_edges and ybin_edges:
            bin_range = [xbin_edges[0], xbin_edges[-1], 
                         ybin_edges[0], ybin_edges[-1]]
        bin_ranges = self.determine_range(bin_range)
        self.xmin, self.xmax = bin_ranges[0:2]
        self.ymin, self.ymax = bin_ranges[2:4]
        self.zmin, self.zmax = bin_ranges[4:6]

        self.nxbins, self.xbin_edges = determine_bins(xbin_edges, xbin_width, nxbins, self.xmin, self.xmax, xvariable)
        self.nybins, self.ybin_edges = determine_bins(ybin_edges, ybin_width, nybins, self.ymin, self.ymax, yvariable)

        self.xunits = xunits
        self.yunits = yunits
        self.zunits = zunits
        self.xlabel, self.ylabel, self.zlabel = self.determine_labels(xlabel, ylabel, zlabel)
        self.style = style

        self.rebin_xbins = []
        self.rebin_ybins = []

    def update(self, region, xvar, yvar):
        self.region = region
        self.xvariable = xvar
        self.yvariable = yvar
        self.name = determine_name(region, xvar, yvar)


    def xbin_width(self):
        return (self.xmax - self.xmin) / self.nxbins

    def ybin_width(self):
        return (self.ymax - self.ymin) / self.nybins

    def determine_range(self, bin_range):
        '''
        Inteligently determine bin range

        The user can provide 4 or 6 values in a list. The default range values
        will be modified. It is assumed providing 4 values indicates one only
        wants to set the x- and y-axes.

        args:
            bin_range (list(int or float)) - list of bin range values

        returns:
            tuple (6 int or float) - range values for both x- and y-axis

        '''
        # default values

        assert len(bin_range) in [4,6],(
            'ERROR :: Unrecognized bin range format:', bin_range)

        if len(bin_range) == 6: return tuple(bin_range)

        if self.doLogZ and self.doNorm:
            zmin = self.logz_norm_min
            zmax = self.logz_norm_max
        elif self.doLogZ:
            zmin = self.logz_min
            zmax = self.logz_max
        elif self.doNorm:
            zmin = self.norm_zmin
            zmax = self.norm_zmax
        else:
            zmin = self.zmin
            zmax = self.zmax

        return bin_range[0], bin_range[1], bin_range[2], bin_range[3], zmin, zmax

    def determine_labels(self, xlabel, ylabel, zlabel):

        # Set x-axis label
        if self.xunits:
            xlabel = "%s [%s]"%(xlabel, self.xunits)
        if self.yunits:
            ylabel = "%s [%s]"%(ylabel, self.yunits)

        # set y-axis label defaults
        if zlabel:
            pass
        elif self.doNorm:
            zlabel = "a.u."
        else:
            zlabel = "Events"

        # set z-axis label
        xwidth_label = str(round(self.xbin_width(), 2))
        ywidth_label = str(round(self.ybin_width(), 2))
        if not self.xunits and xwidth_label == '1.0':
            pass
        elif self.xunits and xwidth_label == '1.0':
            zlabel = "%s / %s"%(zlabel, self.xunits)
        else:
            zlabel = "%s / %s %s"%(zlabel, xwidth_label, self.xunits)

        if not self.yunits and ywidth_label == '1.0':
            pass
        elif self.yunits and ywidth_label == '1.0':
            zlabel = "%s / %s"%(zlabel, self.yunits)
        else:
            zlabel = "%s / %s %s"%(zlabel, ywidth_label, self.yunits)

        return xlabel, ylabel, zlabel

    def setDefaultPads(self, name) :
        self.pads = Pads(name)
        self.ptype = Types.default

    def set2DPads(self, name) :
        self.pads = Pads(name)
        self.ptype = Types.two_dim
    
    def setCutScanPads2D(self, name):
        self.pads = CutScanPads2D(name)

    def setTProfilePads(self, name) :
        self.pads = Pads(name)
        self.ptype = Types.profile

    def make_2d_hist(self, hists, suffix = ""):
        if not hists.hist: return
        # Get canvas
        can = self.pads.canvas
        can.cd()
        if self.doLogZ : can.SetLogz(True)
        can.SetLeftMargin(0.15)
        can.SetBottomMargin(0.15)
        can.SetRightMargin(0.2)

        # Formatting
        if self.doNorm:
            pu.normalize_plot(self.hist)
        if self.auto_set_zlimits:
            reformat_zaxis(self.hist)
        hists.axis.Draw()
        hists.hist.Draw("%s SAME" % self.style)

        can.RedrawAxis()
        can.SetTickx()
        can.SetTicky()
        can.Update()

        self.save_plot(can, suffix)
        # self.pads.canvas.Clear()

    def make_cutscan2d_plot(self, hists, reg_name, roc_curve=False):
        # Get Canvas
        pads = self.pads
        pads.sig_hist_pad.cd()
        if self.doLogZ : pads.sig_hist_pad.SetLogz(True)
        hists.axis.Draw()
        hists.signal.SetFillColor(0)
        hists.signal.SetLineWidth(3)
        hists.signal.Draw("COLZ SAME")
        pu.draw_text(x=0.2, y=0.9, text="Signal: %s" % hists.signal.GetTitle())

        pads.bkgd_hist_pad.cd()
        if self.doLogZ : pads.bkgd_hist_pad.SetLogz(True)
        hists.axis.Draw()
        hists.bkgd.SetFillColor(0)
        hists.bkgd.SetLineWidth(3)
        hists.bkgd.Draw("COLZ SAME")
        pu.draw_text(x=0.2, y=0.9, text="Bkgd: %s" % hists.bkgd.GetTitle())
        
        # Reset axis
        hists.axis.SetTitle("")
        pads.sig_hist_pad.RedrawAxis()
        pads.sig_hist_pad.SetTickx()
        pads.sig_hist_pad.SetTicky()
        pads.sig_hist_pad.Update()
        pads.bkgd_hist_pad.RedrawAxis()
        pads.bkgd_hist_pad.SetTickx()
        pads.bkgd_hist_pad.SetTicky()
        pads.bkgd_hist_pad.Update()

        # Add Eff/Rej plot
        pads.sig_eff_pad.cd()
        hists.axis.Draw()
        hists.signal_eff.SetFillColor(0)
        hists.signal_eff.SetLineWidth(3)
        hists.signal_eff.SetLineWidth(3)
        hists.signal_eff.GetZaxis().SetTitle("% Acceptence")
        hists.signal_eff.Draw("COLZ SAME")
        pads.sig_eff_pad.RedrawAxis()
        pads.sig_eff_pad.Update()

        pads.bkgd_eff_pad.cd()
        hists.axis.Draw()
        hists.bkgd_rej.SetFillColor(0)
        hists.bkgd_rej.SetLineWidth(3)
        hists.bkgd_rej.GetZaxis().SetTitle("% Rejection")
        hists.bkgd_rej.Draw("COLZ SAME")
        pads.bkgd_eff_pad.RedrawAxis()
        pads.bkgd_eff_pad.Update()

        # Add S/B plot
        pads.s_over_b_pad.cd()
        hists.axis.Draw()
        hists.s_over_b.SetFillColor(hists.s_over_b.GetLineColor())
        hists.s_over_b.Draw("COLZ SAME")
        pads.s_over_b_pad.RedrawAxis()
        pads.s_over_b_pad.Update()

        # Add ROC Curve
        pads.roc_pad.cd()
        #pads.roc_pad.SetLogy(True)
        hists.roc_graph.Draw("A*")

        # Save the histogram
        self.save_plot(pads.canvas, suffix="cutScan")

    def reformat_zaxis(self):
        # Get maximum histogram z-value
        maxz = self.hist.GetMaximum()
        minz = self.hist.GetMinimum()
        assert maxz > 0

        # Get default z-axis max and min limits
        if self.doLogZ:
            maxz = 10**(pu.get_order_of_mag(maxz))
            if minz > 0:
                minz = 10**(pu.get_order_of_mag(minz))
            else:
                minz = 10**(pu.get_order_of_mag(maxz) - 7)
        else:
            minz = 0

        # Get z-axis max multiplier to fit labels

        # reformat the axis
        self.hist.SetMaximum(maxz)
        self.hist.SetMinimum(minz)


    def Print(self) :
        print "Plot2D    plot: %s  (region: %s  xVar: %s  yVar: %s)"%(
            self.name, self.region, self.xVariable, self.yVariable)

class Plot3D(PlotBase) :
    def __init__(self,
        region = "",
        name = "",
        xvariable = "",
        yvariable = "",
        zvariable = "",
        xlabel = "x-Label",
        xunits = "",
        ylabel = "y-Label",
        yunits = "",
        zlabel = "z-Label",
        zunits = "",
        suffix = "",
        bin_range = [], # [x0, x1, y0, y1, z0, z1]
        xbin_width = None, # Can specify to override nbins
        ybin_width = None, # Can specify to override nbins
        zbin_width = None, # Can specify to override nbins
        xbin_edges = None, # Can specify to override nbins
        ybin_edges = None, # Can specify to override nbins
        zbin_edges = None, # Can specify to override nbins
        nxbins = 20,
        nybins = 20,
        nzbins = 20,
        add_overflow = False,
        add_underflow = False,
        ptype = Types.three_dim,
        ) :
        # Descriptors
        self.region = region
        self.xvariable = xvariable
        self.yvariable = yvariable
        self.zvariable = zvariable
        self.name = name if name else determine_name(region, xvariable, yvariable, zvariable)
        self.suffix = suffix

        #Flags
        self.is2D = False
        self.is3D = True
        self.add_overflow = add_overflow
        self.add_underflow = add_underflow
        self.ptype = ptype

        # Properties

        if xbin_edges and ybin_edges and zbin_edges:
            bin_range = [xbin_edges[0], xbin_edges[-1], 
                         ybin_edges[0], ybin_edges[-1],
                         zbin_edges[0], zbin_edges[-1]]
        self.xmin, self.xmax = bin_range[0:2]
        self.ymin, self.ymax = bin_range[2:4]
        self.zmin, self.zmax = bin_range[4:6]

        self.nxbins, self.xbin_edges = determine_bins(xbin_edges, xbin_width, nxbins, self.xmin, self.xmax, xvariable)
        self.nybins, self.ybin_edges = determine_bins(ybin_edges, ybin_width, nybins, self.ymin, self.ymax, yvariable)
        self.nzbins, self.zbin_edges = determine_bins(zbin_edges, zbin_width, nzbins, self.zmin, self.zmax, zvariable)

        self.xunits = xunits
        self.yunits = yunits
        self.zunits = zunits
        self.xlabel, self.ylabel, self.zlabel = self.determine_labels(xlabel, ylabel, zlabel)

        self.rebin_xbins = []
        self.rebin_ybins = []
        self.rebin_zbins = []

    def update(self, region, xvar, yvar, zvar):
        self.region = region
        self.xvariable = xvar
        self.yvariable = yvar
        self.zvariable = zvar
        self.name = determine_name(region, xvar, yvar, zvar)

    def determine_labels(self, xlabel, ylabel, zlabel):
        if self.xunits:
            xlabel = "%s [%s]"%(xlabel, self.xunits)
        if self.yunits:
            ylabel = "%s [%s]"%(ylabel, self.yunits)
        if self.zunits:
            zlabel = "%s [%s]"%(zlabel, self.zunits)

        return xlabel, ylabel, zlabel

def determine_name(region, *argv):
    ''' Determine default plot name from variables and region'''
    vars_stripped = [pu.strip_for_root_name(var) for var in argv]
    return "_".join([region] + vars_stripped)

def determine_nbins(ax_max, ax_min, bin_width, variable = "?", update_range = True):
    '''
    Intelligently determine the number of bins.

    Given a desired bin width and axis range, the correct number of bins
    is determined for use with TH1. The axis range or bin width will
    likely be adjusted so that the latter divides the former. If the
    boundaries are int type values, these will be maintained.

    args:
        ax_min (float/int) - minimum axis value
        ax_max (float/int) - maximum axis value
        bin_width (float) - desired width of axis bins
        variable (str) - variable of plot. Used for printing a warning
        update_range - option to update axis range to be divisible by
            bin width as opposed to modying bin width to be a divisor

    returns:
        (int) - the number axis bins
    '''
    assert ax_min != ax_max, ("ERROR :: axis range not set")
    bin_width = float(bin_width)
    ax_range = float(ax_max - ax_min)
    remainder = ax_range % bin_width
    cutoff_range = bin_width - remainder if remainder else 0
    int_bins = (isinstance(ax_min, (int, long))
            and isinstance(ax_max, (int, long)))

    eps = 0.0000001
    notify = cutoff_range > eps and abs(cutoff_range - bin_width) > eps
    if update_range and notify:
        print "WARNING :: bin_width is not a divisor of axis range.",
        print "Variable = %s, Range = [%f,%f], bin_width = %f"%(
                variable, ax_min, ax_max, bin_width)
        print "Expanding range by %f to fit"%cutoff_range
    if update_range and ax_min == 0:
        ax_max += cutoff_range
    elif update_range and ax_min != 0:
        ax_min -= 0.5 * cutoff_range
        ax_max += 0.5 * cutoff_range

    if int_bins:
        ax_min = int(round(ax_min))
        ax_max = int(round(ax_max))

    ax_range = (ax_max - ax_min)

    nbins = int(round( ax_range / bin_width ))

    return nbins

def determine_bins(edges, width, nbins, lo, hi, var):
    assert edges or width or nbins
    # If edges are not provided calculate them
    if not edges:
        # If a bin width is requested, calculate the required number of bins
        if width:
            nbins = determine_nbins(hi, lo, width, var)
        
        edges = pu.determine_bin_edges(lo, hi, nbins)
    else:
        nbins = len(edges) - 1
    return nbins, array('d', edges)
################################################################################
# TPad handler classes
################################################################################
class Pads :
    def __init__(self, name, width=800, height=600):
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, width, height)
        self.set_pad_dimensions()

    def set_pad_dimensions(self):
        pass

class StackPads(Pads):
    def __init__(self, name):
        Pads.__init__(self, name)

    def set_pad_dimensions(self):
        can = self.canvas
        can.cd()

        # Color
        can.SetFrameFillColor(0)
        can.SetFillColor(0)

        # Margins
        can.SetRightMargin(0.05)
        can.SetLeftMargin(0.14)
        can.SetBottomMargin(1.3*can.GetBottomMargin())

        can.Update()
        self.canvas = can

class RatioPads :
    ylabel = 'Data / MC'
    ymax = 2
    def __init__(self, name) :
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, 768, 768)
        self.upper_pad = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
        self.lower_pad = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dimensions()

    def set_pad_dimensions(self) :
        can = self.canvas
        up  = self.upper_pad
        dn  = self.lower_pad

        can.cd()
        up_height = 0.75
        dn_height = 0.30
        up.SetPad(0.0, 1.0-up_height, 1.0, 1.0)
        dn.SetPad(0.0, 0.0, 1.0, dn_height)

        up.SetTickx(0)
        dn.SetGrid(0)
        dn.SetTicky(0)

        up.SetFrameFillColor(0)
        up.SetFillColor(0)

        # set margins
        up.SetRightMargin(0.05)
        up.SetLeftMargin(0.14)
        up.SetTopMargin(0.7 * up.GetTopMargin())
        up.SetBottomMargin(0.09)

        dn.SetRightMargin(up.GetRightMargin())
        dn.SetLeftMargin(up.GetLeftMargin())
        dn.SetBottomMargin(0.4)

        up.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can
        self.upper_pad = up
        self.lower_pad = dn

class DoubleRatioPads :
    def __init__(self, name) :
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, 300, 350)
        self.upper_pad = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
        self.middle_pad = r.TPad("middle", "middle", 0.0, 0.0, 1.0, 1.0)
        self.lower_pad = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dimensions()

    def set_pad_dimensions(self) :
        can = self.canvas
        up  = self.upper_pad
        mid = self.middle_pad
        dn = self.lower_pad

        can.cd()
        up_height = 0.90
        mid_height_low = 0.25
        mid_height_high = 0.40
        dn_height = 0.25

        up.SetPad(0.0, mid_height_high, 1.0, 1.0)
        mid.SetPad(0.0, mid_height_low, 1.0, mid_height_high)
        dn.SetPad(0.0, 0.0, 1.0, mid_height_low)

        up.SetTickx(0)
        mid.SetGrid(0)
        mid.SetTicky(0)
        dn.SetGrid(0)
        dn.SetTicky(0)

        up.SetFrameFillColor(0)
        up.SetFillColor(0)

        # set right margins
        right_margin = 0.05
        up .SetRightMargin(right_margin)
        mid.SetRightMargin(right_margin)
        dn .SetRightMargin(right_margin)

        # set left margins
        left_margin = 0.14
        up .SetLeftMargin(left_margin)
        mid.SetLeftMargin(left_margin)
        dn .SetLeftMargin(left_margin)

        # bottom margins
        up.SetBottomMargin(0.04)
        mid.SetBottomMargin(0.15)
        dn.SetBottomMargin(0.47)

        # set top margins
        up.SetTopMargin(0.09)
        mid.SetTopMargin(0.05)
        dn.SetTopMargin(0.02)

        up.Draw()
        mid.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can
        self.upper_pad = up
        self.middle_pad = mid
        self.lower_pad = dn

class CutScanPads1D(Pads):
    def __init__(self, name, width=600, height=1200):
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, width, height)
        self.hist_pad = r.TPad("hist_pad", "hist_pad", 0.0, 0.0, 1.0, 1.0)
        self.eff_pad = r.TPad("eff_pad", "eff_pad", 0.0, 0.0, 1.0, 1.0)
        self.s_over_b_pad = r.TPad("s_over_b_pad", "s_over_b_pad", 0.0, 0.0, 1.0, 1.0)
        self.roc_pad = r.TPad("roc_pad", "roc_pad", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dimensions()

    def set_pad_dimensions(self):
        can = self.canvas
        up  = self.hist_pad
        mid1 = self.eff_pad
        mid2 = self.s_over_b_pad
        dn = self.roc_pad

        # Define percentage heigh of each pad
        can.cd()
        up_h= 0.30
        mid1_h= 0.20
        mid2_h= 0.20
        dn_h= 0.30

        top = 1.0
        div1 = top - up_h
        div2 = top - up_h - mid1_h
        div3 = top - up_h - mid1_h - mid2_h
        bottom = 0 

        up.SetPad(0.0, div1, 1.0, top)
        mid1.SetPad(0.0, div2, 1.0, div1)
        mid2.SetPad(0.0, div3, 1.0, div2)
        dn.SetPad(0.0, bottom, 1.0, dn_h)

        #up.SetGrid(0)
        #up.SetTickx(0)
        mid1.SetGrid(1)
        #mid1.SetTickx(0)
        mid2.SetGrid(1)
        #mid2.SetTicky(0)
        dn.SetGrid(1)
        #dn.SetTicky(0)

        up.SetFrameFillColor(0)
        up.SetFillColor(0)

        # set right margins
        right_margin = 0.05
        up .SetRightMargin(right_margin)
        mid1.SetRightMargin(right_margin)
        mid2.SetRightMargin(right_margin)
        dn .SetRightMargin(right_margin)

        # set left margins
        up .SetLeftMargin(0.14)
        mid1.SetLeftMargin(0.14)
        mid2.SetLeftMargin(0.14)
        dn .SetLeftMargin(0.14)

        # bottom margins
        up.SetBottomMargin(0.01)
        mid1.SetBottomMargin(0.01)
        mid2.SetBottomMargin(0.2)
        dn.SetBottomMargin(0.15)

        # set top margins
        up.SetTopMargin(0.09)
        mid1.SetTopMargin(0.01)
        mid2.SetTopMargin(0.01)
        dn.SetTopMargin(0.1)
       
        # Add pads to main canvas
        up.Draw()
        mid1.Draw()
        mid2.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can 
        self.hist_pad = up  
        self.eff_pad = mid1 
        self.s_over_b_pad = mid2 
        self.roc_pad = dn 

class CutScanPads2D(Pads):
    def __init__(self, name, width=600, height=1200):
        self.name = "c_" + name
        self.canvas = r.TCanvas(self.name, self.name, width, height)
        self.sig_hist_pad = r.TPad("sig_hist_pad", "sig_hist_pad", 0.0, 0.0, 1.0, 1.0)
        self.bkgd_hist_pad = r.TPad("bkgd_hist_pad", "bkgd_hist_pad", 0.0, 0.0, 1.0, 1.0)
        self.sig_eff_pad = r.TPad("sig_eff_pad", "sig_eff_pad", 0.0, 0.0, 1.0, 1.0)
        self.bkgd_eff_pad = r.TPad("bkgd_eff_pad", "bkgd_eff_pad", 0.0, 0.0, 1.0, 1.0)
        self.s_over_b_pad = r.TPad("s_over_b_pad", "s_over_b_pad", 0.0, 0.0, 1.0, 1.0)
        self.roc_pad = r.TPad("roc_pad", "roc_pad", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dimensions()

    def set_pad_dimensions(self):
        can = self.canvas
        up_l  = self.sig_hist_pad
        up_r  = self.bkgd_hist_pad
        mid1_l = self.sig_eff_pad
        mid1_r = self.bkgd_eff_pad
        mid2 = self.s_over_b_pad
        dn = self.roc_pad

        # Define percentage heigh of each pad
        can.cd()
        up_h= 0.25
        mid1_h= 0.25
        mid2_h= 0.25
        dn_h= 0.25

        top = 1.0
        div1 = top - up_h
        div2 = top - up_h - mid1_h
        div3 = top - up_h - mid1_h - mid2_h
        bottom = 0 

        up_l.SetPad(0.0, div1, 0.5, top)
        up_r.SetPad(0.5, div1, 1.0, top)
        mid1_l.SetPad(0.0, div2, 0.5, div1)
        mid1_r.SetPad(0.5, div2, 1.0, div1)
        mid2.SetPad(0.0, div3, 1.0, div2)
        dn.SetPad(0.0, bottom, 1.0, dn_h)

        #up.SetGrid(0)
        #up.SetTickx(0)
        mid1_l.SetGrid(1)
        mid1_r.SetGrid(1)
        #mid1.SetTickx(0)
        mid2.SetGrid(1)
        #mid2.SetTicky(0)
        dn.SetGrid(1)
        #dn.SetTicky(0)

        up_l.SetFrameFillColor(0)
        up_l.SetFillColor(0)
        up_r.SetFrameFillColor(0)
        up_r.SetFillColor(0)

        # set right margins
        right_margin = 0.2
        up_l.SetRightMargin(right_margin)
        up_r.SetRightMargin(right_margin)
        mid1_l.SetRightMargin(right_margin)
        mid1_r.SetRightMargin(right_margin)
        mid2.SetRightMargin(right_margin)
        dn .SetRightMargin(0.1)

        # set left margins
        left_margin = 0.2
        up_l.SetLeftMargin(left_margin)
        up_r.SetLeftMargin(left_margin)
        mid1_l.SetLeftMargin(left_margin)
        mid1_r.SetLeftMargin(left_margin)
        mid2.SetLeftMargin(0.15)
        dn.SetLeftMargin(0.15)

        # bottom margins
        up_l.SetBottomMargin(0.01)
        up_r.SetBottomMargin(0.01)
        mid1_l.SetBottomMargin(0.15)
        mid1_r.SetBottomMargin(0.15)
        mid2.SetBottomMargin(0.15)
        dn.SetBottomMargin(0.2)

        # set top margins
        up_l.SetTopMargin(0.15)
        up_r.SetTopMargin(0.15)
        mid1_l.SetTopMargin(0.01)
        mid1_r.SetTopMargin(0.01)
        mid2.SetTopMargin(0.01)
        dn.SetTopMargin(0.1)
       
        # Add pads to main canvas
        up_l.Draw()
        up_r.Draw()
        mid1_l.Draw()
        mid1_r.Draw()
        mid2.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can 
        self.sig_hist_pad = up_l  
        self.bkgd_hist_pad = up_r  
        self.sig_eff_pad = mid1_l 
        self.bkgd_eff_pad = mid1_r 
        self.s_over_b_pad = mid2 
        self.roc_pad = dn 
