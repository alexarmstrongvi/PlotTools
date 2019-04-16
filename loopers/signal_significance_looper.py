#!/usr/bin/env python
import sys, os, traceback, argparse
import time
import importlib
from YieldTable import YieldTbl, UncFloat
from math import sqrt
from copy import deepcopy
from collections import defaultdict

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

def main():
    global args
    yld_dict = defaultdict(dict)
    for reg in REGIONS:
        print '\n', 20*'-', "Yields for %s region"%reg.displayname, 20*'-', '\n'

        ########################################################################
        print "Setting EventLists for %s"%reg.name,
        yld_dict[reg.name]['mc_bkgd_yld'] = UncFloat() 
        yld_dict[reg.name]['data_yld'] = UncFloat() 
        yld_dict[reg.name]['mc_sig_yld'] = defaultdict(UncFloat) 
        for sample in SAMPLES :
            if sample.isMC:
                weight_var = sample.weight_str
            elif not sample.isMC and sample.blinded and reg.isSR:
                weight_var = "0"
            else:
                weight_var = ""

            scale_factor = sample.scale_factor if sample.isMC else 1
            list_name = "list_" + reg.name + "_" + sample.name
            sample.set_event_list(reg.tcut, list_name, EVENT_LIST_DIR)
            yld, error = get_yield_and_error(sample.tree, weight_var, scale_factor)
            if sample.name == 'data': error = 0
            
            if sample.isMC and sample.isSignal:
                mass_x = sample.mass_x
                mass_y = sample.mass_y
                yld_dict[reg.name]['mc_sig_yld'][(mass_x, mass_y)] += UncFloat(yld, error)
            elif sample.isMC and not sample.isSignal:
                yld_dict[reg.name]['mc_bkgd_yld'] += UncFloat(yld, error)
            elif not sample.isMC:
                yld_dict[reg.name]['data_yld'] += UncFloat(yld, error)
           
    ofile_str = 'region_name,mc_bkgd_yld,mc_bkgd_yld_unc,data_yld,mass_x,mass_y,sig_yld,sig_unc\n'
    for reg_name, d in yld_dict.iteritems():
        by = d['mc_bkgd_yld'].value
        byu = d['mc_bkgd_yld'].uncertainty
        dy = d['data_yld'].value
        for (mass_x, mass_y), yld in d['mc_sig_yld'].iteritems():
            sy = yld.value
            syu = yld.uncertainty
            row = [reg_name, by, byu, dy, mass_x, mass_y, sy, syu]
            row = [str(x) for x in row]
            ofile_str += "%s\n" % (",".join(row)) 

    ofile_path = "%s/signal_sensativity.txt" % YIELD_TBL_DIR
    with open(ofile_path, 'w') as ofile:
        ofile.write(ofile_str)
    print "INFO :: Output file written to", ofile_path




def get_yield_and_error(ttree, weight_var="", scale=1, dummy_var="isMC"):
    error = r.Double(0.0)
    weight_str = "%s * %f" % (weight_var, scale) if weight_var else "1"
    h_tmp = r.TH1D('h_get_yields_tmp','',1,0,-1) #Beware of saturation errors when using TH1F instead of TH1D. Sample be getting huge these days
    draw_cmd = "%s >> %s" % (dummy_var, h_tmp.GetName())
    ttree.Draw(draw_cmd, weight_str)
    yld = h_tmp.IntegralAndError(0,-1, error)
    h_tmp.Delete()

    return yld, error
        
        
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
    print "     suffix           :  %s "%args.suffix
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
    print " ============================"

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
                                help='Suffix to append to output files')
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
        EVENT_LIST_DIR = conf.EVENT_LIST_DIR
        YIELD_TBL_DIR = conf.YIELD_TBL_DIR

        print_inputs(args)

        main()

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

