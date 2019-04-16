#!/usr/bin/env python
import sys, os, traceback, argparse
import time
import importlib
from YieldTable import YieldTbl, UncFloat
from math import sqrt
from global_variables import event_list_dir, plots_dir, yield_tbl_dir
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

def main():
    global args

    for reg in REGIONS:
        print '\n', 20*'-', "Yields for %s region"%reg.displayname, 20*'-', '\n'

        yld_table = reg.yield_table if reg.yield_table else deepcopy(YLD_TABLE)
        ########################################################################
        for sample in SAMPLES :
            if sample.isMC:
                weight_var = sample.weight_str
            elif not sample.isMC and sample.blinded and reg.isSR:
                weight_var = "0"
            else:
                weight_var = ""
            
            scale_factor = sample.scale_factor if sample.isMC else 1
            list_name = "list_" + reg.name + "_" + sample.name
            sample.tree.SetEventList(0)
            trig_cut = reg.cutflow[1].cut_str            
            for ii, cut in enumerate(reg.cutflow):
                print "Setting eventlist for Cut %2d of %2d: %50s" % (ii+1, len(reg.cutflow), cut.cut_str.replace(trig_cut,"Trigger"))
                sample.set_event_list(cut.cut_str, list_name, event_list_dir, reset=True)
                yld, error = get_yield_and_error(sample.tree, weight_var, scale_factor)
                if sample.name == 'data': error = 0
                #print "TESTING :: Yield = ", yld
                yld_table.add_entry(row_name = cut.name,
                                    col_name = sample.name,
                                    val = yld,
                                    error = error,
                                    row_displayname = cut.displayname,
                                    row_latexname = cut.latexname, 
                                    col_displayname = sample.displayname,
                                    col_latexname = sample.latexname,
                                    ) 

        #TODO
        #yld_table.apply_column_formulas()
        #yld_table.apply_row_formulas()

        if yld_table.write_to_latex:
            name = reg.name + "_cutflow"
            if args.suffix: name += "_" + args.suffix
            save_path = os.path.join(yield_tbl_dir, name + ".tex")
            print "Saving cutflow yield table to", save_path
            yld_table.save_table(save_path, latex=True)

        yld_table.Print()

def get_yield_and_error(ttree, weight_var="", scale=1, dummy_var="isMC"):
    error = r.Double(0.0)
    weight_str = "%s * %f" % (weight_var, scale) if weight_var else "1"
    draw_cmd = "%s >> h_get_yields_tmp" % dummy_var
    #print "TESTING :: weight_str = ", weight_str, "; draw_cmd =", draw_cmd
    ttree.Draw(draw_cmd, weight_str)
    h_tmp = r.gROOT.FindObject('h_get_yields_tmp')
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
                                help='Suffix to append to output yield tables')
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
        YLD_TABLE = conf.YLD_TABLE

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

