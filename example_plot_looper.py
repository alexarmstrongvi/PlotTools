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
import PlotTools.hist as Hist
from PlotTools.YieldTable import UncFloat
from global_variables import event_list_dir, plots_dir

################################################################################
def main ():
    """ Main Function """

    global args
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
            if args.suffix:
                if plot.suffix:
                    plot.suffix += "_" + args.suffix
                else:
                    plot.suffix = args.suffix

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
    #TODO Move combined samples to sample conf file
    mc_samples = [s for s in SAMPLES if s.isMC and not s.isSignal]
    data_samples = [s for s in SAMPLES if not s.isMC]
    signal_samples = [s for s in SAMPLES if s.isMC and s.isSignal]

    with Hist.Hist2D(plot, reg, YIELD_TBL, data_samples) as main_hist:
        plot.make_2d_hist(main_hist, suffix = "_data")
    with Hist.Hist2D(plot, reg, YIELD_TBL, mc_samples) as main_hist:
        plot.make_2d_hist(main_hist, suffix = "_mc")
    with Hist.Hist2D(plot, reg, YIELD_TBL, signal_samples) as main_hist:
        plot.make_2d_hist(main_hist, suffix = "_signal")

################################################################################
def make_plotsStack(plot, reg):
    ''' '''
    backgrounds = [s for s in SAMPLES if s.isMC and not s.isSignal]
    data = [s for s in SAMPLES if not s.isMC][0]
    signals = [s for s in SAMPLES if s.isMC and s.isSignal]
    with Hist.DataMCStackHist1D(plot, reg, YIELD_TBL, data=data, bkgds=backgrounds, sig=signals) as main_hist:
        plot.make_data_mc_stack_plot(reg.displayname, main_hist)

    ## For fake factor histograms
    #with Hist.DataCorrMCRealStackHist1D(plot, reg, data=data, data_corr=bkgds, bkgds=backgrounds) as main_hist:
    #    plot.make_fake_datacorr_mctruth_stack_plot(main_hist)

    #with Hist.FakeFactor1D(plot, reg, ata_num, data_den, mc_num, mc_den) as main_hist:
    #    plot.make_fake_factor_plot(main_hist)

def make_plotsRatio(plot, reg) :
    ''' '''
    backgrounds = [s for s in SAMPLES if s.isMC and not s.isSignal]
    data = [s for s in SAMPLES if not s.isMC][0]
    signals = [s for s in SAMPLES if s.isMC and s.isSignal]
    with Hist.DataMCStackHist1D(plot, reg, YIELD_TBL, data=data, bkgds=backgrounds, sig=signals) as main_hist:
        with Hist.DataMCRatioHist1D(plot, reg, main_hist) as ratio_hist:
            plot.make_data_mc_stack_with_ratio_plot(reg.displayname, main_hist, ratio_hist)

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
        parser.add_argument("-s", "--suffix",
                                default="",
                                help='Suffix to append to output plot names')
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

