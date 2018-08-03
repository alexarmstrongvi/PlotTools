#!/usr/bin/env python
"""
================================================================================
As comprehensive a plotting script as there is this side of the Mississipi

Examples
    ./plotter.py -c config_file.py

Description:
    Creates various histogram plots (e.g. stack, ratio plot, overlay, etc.)
    composed of several samples (i.e. data, background, and signal). The script
    requires a detailed configuration file be provided in that defines
    backgrounds, data, systematics, regions, and plots.


Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with lots of ideas borrowed from Danny Antrim <dantrim@cern.ch>

License:
    Copyright: (C) <May 16th, 2018>; University of California, Irvine
================================================================================
"""
# General python
import sys, os, traceback, argparse
import time
import importlib
import re
from array import array
from collections import OrderedDict
from copy import deepcopy

# Root data analysis framework
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

# Turn off root ownership when creating classes
r.TEventList.__init__._creates = False
r.TH1F.__init__._creates = False
r.TGraphErrors.__init__._creates = False
r.TGraphAsymmErrors.__init__._creates = False
r.TLegend.__init__._creates = False
r.THStack.__init__._creates = False

# Local classes for plotting
import PlotTools.plot_utils as pu
from PlotTools.plot import Types
from PlotTools.YieldTable import UncFloat
from global_variables import event_list_dir, plots_dir

################################################################################
def main ():
    """ Main Function """

    global args
    # Slim plots to make if requested
    rReg = args.requestRegion
    rPlot = args.requestPlot
    if not rReg and not rPlot:
        pass # Use all plots
    elif rReg and not rPlot:
        PLOTS = [p for p in PLOTS if p.region == rReg]
    elif not rReg and rPlot:
        PLOTS = [p for p in PLOTS if p.name == rPlot]
    elif rReg and rPlot:
        PLOTS = [p for p in PLOTS if p.name == rPlot and p.region == rReg]

    make_plots()

    print "Yields found during plotting...\n"
    for yld_tbl in YIELD_TABLES:
        yld_tbl.Print()
        print '\n'

################################################################################
# PLOTTING FUNCTIONS
def make_plots() :
    ''' '''
    for reg in REGIONS:
        # Get plots for this region
        plots_with_region = [p for p in PLOTS if p.region == reg.name]
        if not len(plots_with_region): continue


        print '\n', 20*'-', "Plots for %s region"%reg.displayname, 20*'-', '\n'

        ########################################################################
        print "Setting EventLists for %s"%reg.name
        cut = r.TCut(reg.tcut)
        for sample in SAMPLES :
            list_name = "list_" + reg.name + "_" + sample.name
            sample.set_event_list(cut, list_name, event_list_dir)
        ########################################################################
        # Make 1D and 2D plots
        n_total = len(plots_with_region)
        for n_current, plot in enumerate(plots_with_region, 1):
            print "[%d/%d] Plotting %s"%(n_current, n_total, plot.name), 40*'-'

            # Clear the yield table
            YIELD_TBL.reset()
            YIELD_TBL.region = reg.displayname

            # Determine the correct plot
            if plot.is2D :
                YIELD_TBL.variable = plot.xvariable+":"+plot.yvariable
                make_plots2D(plot, reg)
            else:
                YIELD_TBL.variable = plot.variable
                make_plots1D(plot, reg)

            # Print yield table
            for yld_tbl in YIELD_TABLES:
                if YIELD_TBL == yld_tbl:
                    yld_tbl.variable = ', '.join([YIELD_TBL.variable, yld_tbl.variable])
                    break
            else:
                YIELD_TABLES.append(deepcopy(YIELD_TBL))

    print 30*'-', 'PLOTS COMPLETED', 30*'-','\n'

################################################################################
def make_plots1D(plot, reg) :
    ''' '''
    if plot.ptype == Types.ratio :
        make_plotsRatio(plot, reg)
    elif plot.ptype == Types.comparison :
        make_plotsComparison(plot, reg)
    elif plot.ptype == Types.stack :
        make_plotsStack(plot, reg)
    else:
        print "ERROR :: Unknown plot type,", plot.ptype.name

def make_plots2D(plot, reg):
    # Get canvas
    can = plot.pads.canvas
    can.cd()
    if plot.doLogZ : can.SetLogz(True)
    can.SetLeftMargin(0.15)
    can.SetBottomMargin(0.15)
    can.SetRightMargin(0.2)

    # Get plot primitive
    axis = make_2D_axis(plot)

    mc_samples = [s for s in SAMPLES if s.isMC and not s.isSignal]
    data_samples = [s for s in SAMPLES if not s.isMC]
    signal_samples = [s for s in SAMPLES if s.isMC and s.isSignal]

    mc_plot = make_2D_plot(plot, reg, mc_samples) if mc_samples else None
    data_plot = make_2D_plot(plot, reg, data_samples) if data_samples else None
    sig_plot = make_2D_plot(plot, reg, signal_samples) if signal_samples else None

    # Formatting
    for h, suffix in [(mc_plot,'mc'), (data_plot,'data'), (sig_plot,'signal')]:
        if not h: continue
        if plot.doNorm:
            normalize_plot(h)
        if plot.auto_set_zlimits:
            reformat_zaxis(plot, h)
        axis.Draw()
        h.Draw("%s SAME" % plot.style)

        can.RedrawAxis()
        can.SetTickx()
        can.SetTicky()
        can.Update()

        save_plot(can, plot.name+"_"+suffix+".pdf")
        plot.pads.canvas.Clear()
    root_delete([axis, mc_plot, data_plot, sig_plot])

def make_2D_plot(plot, reg, samples):
    h = r.TH2D(plot.name, "", plot.nxbins, plot.xmin, plot.xmax, plot.nybins, plot.ymin, plot.ymax)
    for sample in samples:

        x_var = re.sub(r'[(){}[\]]+','',plot.xvariable)
        y_var = re.sub(r'[(){}[\]]+','',plot.yvariable)
        h_name_tmp = "h_"+reg.name+'_'+sample.name+"_"+x_var+"_"+y_var
        h_tmp = r.TH2D(h_name_tmp, "", plot.nxbins, plot.xmin, plot.xmax, plot.nybins, plot.ymin, plot.ymax)
        # Draw final histogram (i.e. selections and weights applied)


        if not sample.isMC:
            weight_str = '1'
        elif plot.xvariable != sample.weight_str and plot.yvariable != sample.weight_str:
            weight_str = "%s * %s"%(sample.weight_str, str(sample.scale_factor))
        else:
            weight_str = '1'
        draw_cmd = "%s>>%s"%(plot.yvariable+":"+plot.xvariable, h_tmp.GetName())
        sample.tree.Draw(draw_cmd, weight_str, "goff")

        # Yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h_tmp.IntegralAndError(0,-1,0,-1,stat_err)
        if sample.isMC and sample.isSignal:
            YIELD_TBL.signals[sample.name] = UncFloat(integral, stat_err)
        elif sample.isMC and not sample.isSignal:
            YIELD_TBL.mc[sample.name] = UncFloat(integral, stat_err)
        elif not sample.isMC:
            YIELD_TBL.data[sample.name] = UncFloat(integral, stat_err)

        h.Add(h_tmp)

    zax = h.GetZaxis()
    zax.SetTitle(plot.zlabel)
    zax.SetTitleFont(42)
    zax.SetLabelFont(42)
    zax.SetTitleOffset(1.5)
    zax.SetLabelOffset(0.013)
    zax.SetLabelSize(1.2 * 0.035)
    return h

def normalize_plot(hist):
    if hist and hist.Integral():
        hist.Scale(1.0/hist.Integral())

def reformat_zaxis(plot, h):
    # Get maximum histogram z-value
    maxz = h.GetMaximum()
    minz = h.GetMinimum()
    assert maxz > 0

    # Get default z-axis max and min limits
    if plot.doLogZ:
        maxz = 10**(pu.get_order_of_mag(maxz))
        if minz > 0:
            minz = 10**(pu.get_order_of_mag(minz))
        else:
            minz = 10**(pu.get_order_of_mag(maxz) - 7)
    else:
        minz = 0

    # Get z-axis max multiplier to fit labels

    # reformat the axis
    h.SetMaximum(maxz)
    h.SetMinimum(minz)

def make_2D_axis(plot):
    hax = r.TH2D("axes", "", plot.nxbins, plot.xmin, plot.xmax, plot.nybins, plot.ymin, plot.ymax)
    hax.SetMinimum(plot.zmin)
    hax.SetMaximum(plot.zmax)

    xax = hax.GetXaxis()
    xax.SetTitle(plot.xlabel)
    xax.SetTitleFont(42)
    xax.SetLabelFont(42)
    xax.SetLabelSize(0.035)
    xax.SetTitleSize(0.048 * 0.85)
    xax.SetLabelOffset(1.15 * 0.02)
    xax.SetTitleOffset(1.5 * xax.GetTitleOffset())

    #if plot.bin_labels:
    #    plot.set_bin_labels(hax)

    yax = hax.GetYaxis()
    yax.SetTitle(plot.ylabel)
    yax.SetTitleFont(42)
    yax.SetLabelFont(42)
    yax.SetTitleOffset(1.4)
    yax.SetLabelOffset(0.013)
    yax.SetLabelSize(1.2 * 0.035)
    yax.SetTitleSize(0.055 * 0.85)

    zax = hax.GetZaxis()
    zax.SetTitle(plot.zlabel)
    #zax.SetTitleFont(42)
    #zax.SetLabelFont(42)
    #zax.SetTitleOffset(1.4)
    #zax.SetLabelOffset(0.013)
    #zax.SetLabelSize(1.2 * 0.035)
    #zax.SetTitleSize(0.055 * 0.85)

    #if plot.bin_labels:
    #    plot.set_ybin_labels(hax)

    #if plot.rebin_xbins:
    #    new_bins = array('d', plot.rebin_xbins)
    #    hax = hax.RebinX(len(new_bins)-1, 'axes', new_bins)
    #if plot.rebin_ybins:
    #    new_bins = array('d', plot.rebin_ybins)
    #    hax = hax.RebinY(len(new_bins)-1, 'axes', new_bins)

    return hax
################################################################################
def make_plotsStack(plot, reg):
    ''' '''
    # Get Canvas
    can = plot.pads.canvas
    can.cd()
    if plot.doLogY : can.SetLogy(True)

    # Get plotting primitives
    # Not all primitives are for drawing. Some are for preserving pointers
    legend = make_stack_legend(plot)
    axis = make_stack_axis(plot)
    mc_stack, mc_total, signals, hists = add_stack_backgrounds(plot, reg)
    data, data_hist = add_stack_data(plot, legend, reg)
    error_leg, mc_errors = add_stack_mc_errors(plot, legend, hists, mc_stack)
    if plot.doNorm:
        normalize_stack(mc_total, signals, data_hist, data, mc_stack, mc_errors)
    if plot.auto_set_ylimits:
        reformat_axis(plot, legend, mc_stack, data, axis, signals)

    # Checks - Move on to next plot in case of failure
    if not mc_stack:
        print "WARNING :: Stack plot has either no MC. Skipping."
        return


    # Draw the histograms
    draw_stack(axis, mc_stack, mc_errors, mc_total, signals, data, legend, reg.displayname, plot.doNorm)

    # Reset axis
    can.RedrawAxis()
    can.SetTickx()
    can.SetTicky()
    can.Update()

    # Save the histogram
    save_plot(can, plot.name+".pdf")

    # Clean up
    root_delete([axis, mc_stack, mc_errors, mc_total, error_leg, signals, data, data_hist, legend, hists])
    plot.pads.canvas.Clear()

def make_plotsRatio(plot, reg) :
    ''' '''
    # Pads
    rcan = plot.pads

    ############################################################################
    # Top stack plot
    rcan.upper_pad.cd()
    if plot.doLogY : rcan.upper_pad.SetLogy(True)
    rcan.upper_pad.Update()

    # Get plotting primitives
    # Not all primitives are for drawing. Some are for preserving pointers
    legend = make_stack_legend(plot)
    axis = make_stack_axis(plot)
    mc_stack, mc_total, signals, hists = add_stack_backgrounds(plot, reg)
    data, data_hist = add_stack_data(plot, legend, reg)
    error_leg, mc_errors = add_stack_mc_errors(plot, legend, hists, mc_stack)
    if plot.doNorm:
        normalize_stack(mc_total, signals, data_hist, data, mc_stack, mc_errors)
    if plot.auto_set_ylimits:
        reformat_axis(plot, legend, mc_stack, data, axis, signals)

    # Checks - Move on to next plot in case of failure
    if not mc_stack and data:
        print "WARNING :: Ratio plot has either no MC or not data. Skipping."
        return


    # Draw the histograms
    draw_stack(axis, mc_stack, mc_errors, mc_total, signals, data, legend, reg.displayname, plot.doNorm)

    # Reset axis
    rcan.upper_pad.RedrawAxis()
    rcan.upper_pad.SetTicks()
    rcan.upper_pad.Update()

    ############################################################################
    # Bottom ratio plot
    rcan.lower_pad.cd()

    ratio_axis = get_ratio_axis(plot, mc_stack, rcan.ylabel, rcan.ymax)
    ratio_errors = get_ratio_errors(mc_errors)
    ratio = get_ratio_graph(data_hist, mc_errors)

    draw_ratio(plot, ratio_axis, ratio_errors, ratio)

    rcan.lower_pad.SetTicks()
    rcan.lower_pad.Update()

    ############################################################################
    # Save the histogram
    save_plot(rcan.canvas, plot.name+".pdf")

    # Clean up
    root_delete([axis, mc_stack, mc_errors, mc_total, error_leg, signals, data,
                 data_hist, legend, hists, ratio_axis, ratio_errors, ratio])
    plot.pads.canvas.Clear()

def save_plot(can, outname):
    save_path = os.path.join(plots_dir, args.outdir, outname)
    save_path = os.path.normpath(save_path)
    can.SaveAs(save_path)
    #OFILE.cd()
    #can.Write()

def root_delete(root_objects):
    for ro in root_objects:
        if not ro or isinstance(ro, r.THStack):
            continue
        elif isinstance(ro, list):
            root_delete(ro)
        elif issubclass(type(ro), r.TObject):
            ro.Delete()
        else:
            print "ERROR :: Unknown primitive type", type(ro)

################################################################################
#Stack Plot Functions
################################################################################
#TODO Clean up and optimize the methods for stack
def make_stack_legend(plot):
    if plot.leg_is_left :
        leg = pu.default_legend(xl=0.2,yl=0.7,xh=0.47, yh=0.87)
    elif plot.leg_is_bottom_right :
        leg = pu.default_legend(xl=0.7, yl=0.17,xh=0.97,yh=0.41)
    elif plot.leg_is_bottom_left :
        leg = pu.default_legend(xl=0.2,yl=0.2,xh=0.47,yh=0.37)
    else :
        leg = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)

    leg.SetNColumns(2)
    # TODO: Incorporate signal legend
    leg_sig = pu.default_legend(xl=0.55, yl=0.6, xh=0.91, yh=0.71)
    leg_sig.SetNColumns(1)

    return leg

def make_stack_axis(plot):
    hax = r.TH1F("axes", "", int(plot.nbins), plot.xmin, plot.xmax)
    hax.SetMinimum(plot.ymin)
    hax.SetMaximum(plot.ymax)
    xax = hax.GetXaxis()
    xax.SetTitle(plot.xlabel)
    xax.SetTitleFont(42)
    xax.SetLabelFont(42)
    xax.SetLabelSize(0.035)
    xax.SetTitleSize(0.048 * 0.85)
    if plot.ptype == Types.ratio:
        hax.GetXaxis().SetTitleOffset(-999)
        hax.GetXaxis().SetLabelOffset(-999)
    else:
        xax.SetLabelOffset(1.15 * 0.02)
        xax.SetTitleOffset(1.5 * xax.GetTitleOffset())

    yax = hax.GetYaxis()
    yax.SetTitle(plot.ylabel)
    yax.SetTitleFont(42)
    yax.SetLabelFont(42)
    yax.SetTitleOffset(1.4)
    yax.SetLabelOffset(0.013)
    yax.SetLabelSize(1.2 * 0.035)
    yax.SetTitleSize(0.055 * 0.85)

    if plot.bin_labels and plot.ptype == Types.stack:
        plot.set_bin_labels(hax)
    if plot.rebin_bins:
        new_bins = array('d', plot.rebin_bins)
        hax = hax.Rebin(len(new_bins)-1, 'axes', new_bins)

    return hax

def add_stack_backgrounds(plot, reg):
    stack = r.THStack("stack_"+plot.name, "")

    # Initilize lists and defaults
    mc_samples = [s for s in SAMPLES if s.isMC]
    histos = []
    all_histos = []
    sig_histos = []
    avoid_bkg = []

    hists_to_clear = []

    # Make MC sample hists
    for mc_sample in mc_samples :
        # Initilize histogram
        h_name_tmp = re.sub(r'[(){}[\]]+','',plot.variable)
        h_name = "h_"+reg.name+'_'+mc_sample.name+"_"+h_name_tmp
        hists_to_clear.append(h_name)
        h = pu.th1d(h_name, "", int(plot.nbins),
                    plot.xmin, plot.xmax,
                    plot.xlabel, plot.ylabel)

        h.SetLineColor(mc_sample.color)
        h.GetXaxis().SetLabelOffset(-999)
        if mc_sample.isSignal:
            h.SetLineWidth(2)
            h.SetLineStyle(2)
            h.SetFillStyle(0)
        else:
            h.SetFillColor(mc_sample.color)
            h.SetFillStyle(1001)
        h.Sumw2

        # Draw final histogram (i.e. selections and weights applied)
        if plot.variable != mc_sample.weight_str:
            weight_str = "%s * %s"%(mc_sample.weight_str, str(mc_sample.scale_factor))
        else:
            weight_str = 1
        cut = "(%s) * %s"%(reg.tcut, weight_str)
        cut = r.TCut(cut)
        sel = r.TCut("1")
        draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
        mc_sample.tree.Draw(draw_cmd, cut * sel, "goff")


        # Yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1,stat_err)

        # Rebin
        if plot.rebin_bins:
            new_bins = array('d', plot.rebin_bins)
            h = h.Rebin(len(new_bins)-1, h_name, new_bins)

        h.leg_name = mc_sample.displayname #dynamic class members...ooo yeah!

        # Add overflow
        if plot.add_overflow:
            pu.add_overflow_to_lastbin(h)
        if plot.add_underflow:
            pu.add_underflow_to_firstbin(h)

        # Record all and non-empty histograms
        if mc_sample.isSignal:
            leg_sig.AddEntry(h, mc_sample.displayname, "l")
            sig_histos.append(h)
            YIELD_TBL.signals[mc_sample.name] = UncFloat(integral, stat_err)
        else:
            all_histos.append(h)
            histos.append(h) if integral > 0 else avoid_bkg.append(mc_sample.name)
            YIELD_TBL.mc[mc_sample.name] = UncFloat(integral, stat_err)

    if not len(histos):
        print "ERROR (make_stack_background) :: All SM hists are empty. Skipping"
        return None, None, sig_histos, all_histos

    # Order the hists by total events
    histos = sorted(histos, key=lambda h: h.Integral())

    for h in histos :
        stack.Add(h)

    h_leg = sorted(all_histos, key=lambda h: h.Integral(), reverse=True)
    histos_for_leg = histos_for_legend(h_leg)

    # draw the total bkg line
    hist_sm = stack.GetStack().Last().Clone("hist_sm")
    hist_sm.SetLineColor(r.kBlack)
    hist_sm.SetLineWidth(3)
    hist_sm.SetLineStyle(1)
    hist_sm.SetFillStyle(0)
    hist_sm.SetLineWidth(3)

    return stack, hist_sm, sig_histos, histos_for_leg

def add_stack_data(plot, leg, reg):
    #TODO: Look for a way to combine with backgrounds
    data = [s for s in SAMPLES if not s.isMC]
    assert len(data) <= 1, "ERROR :: Multiple data samples"
    if len(data) == 0: return None, None
    data = data[0]

    hd_name = "h_"+reg.name+'_data_'+plot.variable
    hd = pu.th1d(hd_name, "", int(plot.nbins),
                              plot.xmin, plot.xmax,
                              plot.xlabel, plot.ylabel)
    hd.Sumw2

    cut = "(" + reg.tcut + ")"
    cut = r.TCut(cut)
    sel = r.TCut("1")
    draw_cmd = "%s>>%s"%(plot.variable, hd.GetName())
    data.tree.Draw(draw_cmd, cut * sel, "goff")
    hd.GetXaxis().SetLabelOffset(-999)

    # print the yield +/- stat error
    stat_err = r.Double(0.0)
    integral = hd.IntegralAndError(0,-1,stat_err)
    YIELD_TBL.data[data.name] = UncFloat(integral, stat_err)

    # Rebin
    if plot.rebin_bins:
        new_bins = array('d', plot.rebin_bins)
        hd = hd.Rebin(len(new_bins)-1, hd_name, new_bins)

    # Add overflow
    if plot.add_overflow:
        pu.add_overflow_to_lastbin(hd)
    if plot.add_underflow:
        pu.add_underflow_to_firstbin(hd)

    gdata = pu.convert_errors_to_poisson(hd)
    #gdata.SetLineWidth(2)
    #uglify
    gdata.SetLineWidth(1)
    gdata.SetMarkerStyle(20)
    gdata.SetMarkerSize(1.5)
    gdata.SetLineColor(1)
    leg.AddEntry(gdata, "Data", "p")

    return gdata, hd

def add_stack_mc_errors(plot, leg, hists, stack):
    if not stack:
        return None, None
    r.gStyle.SetHatchesSpacing(0.9)

    mcError = r.TH1F("mcError", "mcError", 2,0,2)
    mcError.SetLineWidth(3)
    mcError.SetFillStyle(3345)
    mcError.SetFillColor(r.kBlue)
    mcError.SetLineColor(r.kBlack)
    leg.AddEntry(mcError, "Standard Model", "fl")

    # now add backgrounds to legend
    for h in hists :
        if h.Integral(0, -1) <= 0: continue
        leg.AddEntry(h, h.leg_name, "f")

    totalSM = stack.GetStack().Last().Clone("totalSM")
    nominalAsymErrors = pu.th1_to_tgraph(totalSM)
    totalSM.Delete()
    nominalAsymErrors.SetMarkerSize(0)
    nominalAsymErrors.SetLineWidth(0)
    nominalAsymErrors.SetFillStyle(3345)
    nominalAsymErrors.SetFillColor(r.kBlue)

        # symmetrize the errors
    for i in xrange(nominalAsymErrors.GetN()) :
        ehigh = nominalAsymErrors.GetErrorYhigh(i)
        elow  = nominalAsymErrors.GetErrorYlow(i)


        error_sym = r.Double(0.0)
        error_sym = (ehigh + elow) / 2.0

        if ehigh != error_sym:
            print "initial error (+%.2f,-%.2f), symmetrized = (+%.2f,-%.2f)"%(ehigh,elow, error_sym, error_sym)


        nominalAsymErrors.SetPointEYhigh(i,0.0)
        nominalAsymErrors.SetPointEYhigh(i, error_sym)
        nominalAsymErrors.SetPointEYlow(i,0.0)
        nominalAsymErrors.SetPointEYlow(i,error_sym)

    return mcError, nominalAsymErrors

def normalize_stack(mc_total, signals, data_hist, data_graph, mc_stack, mc_errors):
    if mc_total and mc_total.Integral():
        mc_norm_factor = 1.0/mc_total.Integral()
        pu.scale_thstack(mc_stack, mc_norm_factor)
        mc_total.Scale(mc_norm_factor)
        pu.scale_tgraph(mc_errors, mc_norm_factor)
    for s in signals:
        sig_norm_factor = 1.0/s.Integral() if s.Integral() else 1
        s.Scale(sig_norm_factor)
    if data_hist and data_hist.Integral():
        data_norm_factor = 1.0/data_hist.Integral()
        pu.scale_tgraph(data_graph, data_norm_factor)

def reformat_axis(plot, leg, stack, data, hax, signals):
    ''' Reformat axis to fit content and labels'''
    if not stack:
        return
    # Get maximum histogram y-value
    if data:
        maxy = max(pu.get_tgraph_max(data), stack.GetMaximum())
        miny = min(pu.get_tgraph_min(data), stack.GetMinimum())
    else:
        maxy = stack.GetMaximum()
        miny = stack.GetMinimum()
    if maxy <= 0:
        print "WARNING :: Max value of plot is <= 0"
        return

    # Get default y-axis max and min limits
    logy = plot.doLogY
    if logy:
        maxy = 10**(pu.get_order_of_mag(maxy))
        if miny > 0:
            miny = 10**(pu.get_order_of_mag(miny))
        else:
            miny = 10**(pu.get_order_of_mag(maxy) - 4)
    else:
        maxy = maxy
        miny = 0

    # Get y-axis max multiplier to fit labels
    if logy:
        max_mult = 1e6 if signals else 1e5
    else:
        max_mult = 2.0 if signals else 1.8

    # reformat the axis
    stack.SetMaximum(max_mult*maxy)
    stack.SetMinimum(miny)
    hax.SetMaximum(max_mult*maxy)
    hax.SetMinimum(miny)

def draw_stack(axis, mc_stack, mc_errors, mc_total, signals, data, legend, reg_name, do_norm) :
    axis.Draw()
    mc_stack.Draw("HIST SAME")
    mc_errors.Draw("E2 same")
    mc_total.Draw('hist same')
    for hsig in signals: hsig.Draw("hist same")
    if data: data.Draw("option same pz 0")
    legend.Draw()
    pu.draw_atlas_label('Internal','Higgs LFV', reg_name)

################################################################################
#Stack Plot Functions
################################################################################
def get_ratio_axis(plot, stack, ylabel, ymax):
    # yaxis
    h_sm = stack.GetStack().Last().Clone("h_sm")
    yax = h_sm.GetYaxis()
    yax.SetRangeUser(0,ymax)
    yax.SetTitle(ylabel)
    yax.SetTitleSize(0.14 * 0.83)
    yax.SetLabelSize(0.13 * 0.81)
    yax.SetLabelOffset(0.98 * 0.013 * 1.08)
    yax.SetTitleOffset(0.45 * 1.2)
    yax.SetLabelFont(42)
    yax.SetTitleFont(42)
    yax.SetNdivisions(5)

    # xaxis
    xax = h_sm.GetXaxis()
    xax.SetTitleSize(1.1 * 0.14)
    xax.SetLabelSize(yax.GetLabelSize())
    xax.SetLabelOffset(1.15*0.02)
    xax.SetTitleOffset(0.85 * xax.GetTitleOffset())
    xax.SetLabelFont(42)
    xax.SetTitleFont(42)

    h_sm.SetTickLength(0.06)

    if plot.bin_labels and plot.ptype == Types.ratio:
        plot.set_bin_labels(h_sm)

    #if plot.rebin_bins:
    #    print "WARNING :: rebinning not yet implemented"
    #    #TODO: Implement rebinning

    return h_sm

def get_ratio_errors(mc_errors):
    ratio_errors = r.TGraphAsymmErrors(mc_errors)
    pu.buildRatioErrorBand(mc_errors, ratio_errors)

    return ratio_errors

def get_ratio_graph(data_hist, mc_errors):
    g_data = pu.convert_errors_to_poisson(data_hist)
    #g_sm = pu.th1_to_tgraph(h_sm)
    #g_ratio = pu.tgraphAsymmErrors_divide(g_data, g_sm)

    # For Data/MC only use the statistical error for data
    # since we explicity draw the MC error band
    nominalAsymErrorsNoSys = r.TGraphAsymmErrors(mc_errors)
    for i in xrange(nominalAsymErrorsNoSys.GetN()) :
        nominalAsymErrorsNoSys.SetPointError(i-1,0,0,0,0)
    ratio_raw = pu.tgraphAsymmErrors_divide(g_data, nominalAsymErrorsNoSys)
    #ratio_raw = pu.tgraphAsymmErrors_divide(nominalAsymErrorsNoSys,g_data)

    # Make final ratio plot
    ratio = r.TGraphAsymmErrors()
    x1, y1 = r.Double(0.0), r.Double(0.0)
    index = 0
    for i in xrange(ratio_raw.GetN()) :
        ratio_raw.GetPoint(i, x1, y1)
        if y1 <= 0. : y1 = r.Double(-2.0)

        ratio.SetPoint(index, x1, y1)
        xlo, xhi = ratio_raw.GetErrorXlow(i), ratio_raw.GetErrorXhigh(i)
        ylo, yhi = ratio_raw.GetErrorYlow(i), ratio_raw.GetErrorYhigh(i)
        ratio.SetPointError(index, xlo, xhi, ylo, yhi)
        index+=1

    # Format
    #ratio.SetLineWidth(2)
    #uglify
    ratio.SetLineWidth(1)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.5)
    ratio.SetLineColor(1)
    ratio.SetMarkerSize(1.5)

    # Clean up
    ratio_raw.Delete()
    nominalAsymErrorsNoSys.Delete()
    g_data.Delete()

    return ratio

def draw_ratio(plot, axis, ratio_errors, ratio):
    axis.Draw("AXIS")
    ratio_errors.Draw("E2")
    ratio.Draw("option same pz 0")

    xmin, xmax = plot.xmin, plot.xmax
    pu.draw_line(xmin, 1.5, xmax, 1.5, style = 3, width = 1)
    pu.draw_line(xmin, 1.0, xmax, 1.0, style = 2, width = 1, color = r.kBlack)
    pu.draw_line(xmin, 0.5, xmax, 0.5, style = 3, width = 1)

################################################################################
def histos_for_legend(histos) :
    '''
    rearrange histogram list for legend

    param:
        histos : list(TH1F)
            histograms in plot ordered by total events (assumes len >= 4)

    returns:
        list(TH1F)
    '''


    if len(histos) == 7 :
        indices = [0, 4, 1, 5, 2, 6, 3]
    elif len(histos) == 6 :
        indices = [0, 3, 1, 4, 2, 5]
    elif len(histos) == 5 :
        indices = [0, 3, 1, 4, 2]
    elif len(histos) == 4 :
        indices = [0, 2, 1, 3]
    else:
        return histos

    return [histos[idx] for idx in indices]

################################################################################
# SETUP FUNCTIONS
def check_args(args):
    """ Check the input arguments are as expected """
    configuration_file = os.path.normpath(args.plotConfig)
    if not os.path.exists(configuration_file):
        print "ERROR :: Cannot find config file:", configuration_file
        sys.exit()

def check_environment():
    """ Check if the shell environment is setup as expected """
    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    if python_ver < 2.7:
        print "ERROR :: Running old version of python\n", sys.version
        sys.exit()

def print_inputs(args):
    """ Print the program inputs """
    full_path = os.path.abspath(__file__)
    prog_name = os.path.basename(full_path)
    prog_dir = os.path.dirname(full_path)

    print " ==================================================================\n"
    print " Program : %s "%prog_name
    print " Run from: %s "%prog_dir
    print ""
    print " Options:"
    print "     plot config      :  %s "%args.plotConfig
    print "     specific plot    :  %s "%args.requestPlot
    print "     specific region  :  %s "%args.requestRegion
    print "     output directory :  %s "%args.outdir
    print "     verbose          :  %s "%args.verbose
    print ""
    print "===================================================================\n"

    # print out the loaded samples and plots
    print " ============================"
    if SAMPLES :
        print "Loaded samples:    "
        for sample in SAMPLES :
            print '\t',
            sample.Print()
    if args.verbose:
        print "Loaded plots:"
        for plot in PLOTS :
            plot.Print()
    print " ============================"

def check_for_consistency() :
    '''
    Make sure that the plots are not asking for undefined region

    param:
        plots : list(plot class)
            plots defined in config file
        regions : list(region class)
            regions defined in config file
    '''
    region_names = [r.name for r in REGIONS]
    bad_regions = set([p.region for p in PLOTS if p.region not in region_names])
    if len(bad_regions) > 0 :
        print 'check_for_consistency ERROR    '\
        'You have configured a plot for a region that is not defined. '\
        'Here is the list of "bad regions":'
        print bad_regions
        print 'check_for_consistency ERROR    The regions that are defined in the configuration ("%s") are:'%args.plotConfig
        print region_names
        print "check_for_consistency ERROR    Exiting."
        sys.exit()

################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-c", "--plotConfig",
                                default="",
                                help='name of the config file')
        parser.add_argument("-r", "--requestRegion",
                                default="",
                                help='request a region to plot -- will make all plots in the config that are in this region')
        parser.add_argument("-p", "--requestPlot",
                                default="",
                                help='request a specific plot -- provide the name of the plot')
        parser.add_argument("-o", "--outdir",
                                default="./",
                                help='name of the output directory to save plots.')
        parser.add_argument("-v", "--verbose",
                                action="store_true",
                                help='set verbosity mode')
        args = parser.parse_args()

        if args.verbose:
            print '>'*40
            print 'Running {}...'.format(os.path.basename(__file__))
            print time.asctime()

        check_args(args)
        check_environment()

        # Import configuration file
        import_conf = args.plotConfig.replace(".py","")
        conf = importlib.import_module(import_conf)

        SAMPLES = conf.SAMPLES
        REGIONS = conf.REGIONS
        PLOTS = conf.PLOTS
        YIELD_TBL = conf.YIELD_TBL
        YIELD_TABLES = []

        check_for_consistency()
        print_inputs(args)

        #tfile_path = os.path.join(plots_dir, args.outdir, 'plotter_hists.root')
        #OFILE = r.TFile(tfile_path,'RECREATE')
        main()
        #OFILE.Close()

        if args.verbose:
            print time.asctime()
            time = (time.time() - start_time)
            print 'TOTAL TIME: %fs'%time,
            print ''
            print '<'*40
    except KeyboardInterrupt, e: # Ctrl-C
        print 'Program ended by keyboard interruption'
        raise e
    except SystemExit, e: # sys.exit()
        print 'Program ended by system exit'
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)

