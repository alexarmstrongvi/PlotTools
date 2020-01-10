#!/usr/bin/env python

################################################################################
# For each region, determine the composition from different fake categories
# such as double fake, leading fake, etc.
# For each category, determine the sample and fake composition
################################################################################

import sys, os, traceback, argparse
import time
import importlib
from PlotTools.YieldTable import YieldTbl, UncFloat
UncFloat.precision = 2
from math import sqrt
from copy import deepcopy
from collections import defaultdict, namedtuple
from IFFTruthClassifierDefs import IFF_Type

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
    if args.type == '2lep':
        fake_event = '(2 <= truthLepOrderType && truthLepOrderType <= 4)'
        leadfake_event  = 'truthLepOrderType == 2'
        sleadfake_event = 'truthLepOrderType == 3'
        multifake_event = 'truthLepOrderType == 4'
        fake_sels = [
            ('TotalYld' , ''),
            ('FakeEvents' , fake_event),
            ('LeadingFakeLep' , leadfake_event),
            ('SubleadingFakeLep' , sleadfake_event),
            ('MultipleFakeLep' , multifake_event),
        ]
    elif args.type == '3lep':
        fake_event = '(6 <= truthLepOrderType && truthLepOrderType <= 12)'
        multifake_event = '(9 <= truthLepOrderType && truthLepOrderType <= 12)'
        Zfake_event = '!(%s) && (lepIsFNP[ZLepIdx[0]] || lepIsFNP[ZLepIdx[1]])' % multifake_event
        probefake_event = '!(%s) && lepIsFNP[probeLepIdx[0]]' % multifake_event
        fake_sels = [
            ('TotalYld' , ''),
            ('FakeEvents' , fake_event),
            ('MultipleFakeLep' , multifake_event),
            ('ZFakeLep' , Zfake_event),
            ('ProbeFakeLep' , probefake_event),
        ]
    else:
        print "ERROR :: Unknown event type:", args.type
        sys.exit()
    total_total_yld = 0
    for reg in REGIONS:
        print '\n', 20*'-', "Yields for %s region"%reg.displayname, 20*'-', '\n'
        final_print_str = ''
        for name, truth_cut in fake_sels:

            ########################################################################
            print "Setting EventLists for %s + %s"% (reg.name, name)
            total_yld = UncFloat()
            sample_yld = {}
            faketype_yld = defaultdict(UncFloat)
            for sample in SAMPLES :
                if not sample.isMC: continue
                weight_var = sample.weight_str

                scale_factor = sample.scale_factor if sample.isMC else 1
                list_name = "list_" + reg.name + "_" + sample.name
                cut = reg.tcut + " && " + truth_cut if truth_cut else reg.tcut
                sample.set_event_list(cut, list_name, EVENT_LIST_DIR)
                yld, error = get_yield_and_error(sample.tree, weight_var, scale_factor)
                
                result = UncFloat(yld, error)
                total_yld += result
                sample_yld[sample.name] = result
                h_truth = get_truth_hist(sample.tree, weight_var, scale_factor)
                for faketype, yld in get_faketype_ylds(h_truth).items():
                    faketype_yld[faketype] += yld
                h_truth.Delete()
                
                # Done looping over samples
            print_str = "Breakdown for %s\n" % name
            if name == "TotalYld":
                print_str += "Total Yield : %s\n" % str(total_yld)
                total_total_yld = total_yld
            else:
                print_str += "Total Yield : %s [%s %%]\n" % (str(total_yld), str((total_yld/total_total_yld) * 100))
            print_str += "\tFake processes: \n%s" % rank_ylds_str(sample_yld, tabs='\t')
            print_str += "\tFake lepton types: \n%s" % rank_ylds_str(faketype_yld, tabs='\t')

            final_print_str += print_str
        print final_print_str

def get_yield_and_error(ttree, weight_var="", scale=1, dummy_var="isMC"):
    error = r.Double(0.0)
    weight_str = "%s * %f" % (weight_var, scale) if weight_var else "1"
    h_tmp = r.TH1D('h_get_yields_tmp','',1,0,-1) #Beware of saturation errors when using TH1F instead of TH1D. Samples be getting huge these days
    draw_cmd = "%s >> %s" % (dummy_var, h_tmp.GetName())
    ttree.Draw(draw_cmd, weight_str)
    yld = h_tmp.IntegralAndError(0,-1, error)
    h_tmp.Delete()

    return yld, error
        
def get_truth_hist(ttree, weight_var="", scale=1):
    weight_str = "%s * %f" % (weight_var, scale) if weight_var else "1"
    nbins = len(IFF_Type)
    xlow = min(IFF_Type.values()) - 0.5
    xhigh = max(IFF_Type.values()) + 0.5
    h_name = ttree.GetName()+"_hist"
    h_tmp = r.TH1D(h_name,'', nbins, xlow, xhigh)
    draw_cmd = "fnpLepTruthIFFClass >> %s" % (h_tmp.GetName())
    ttree.Draw(draw_cmd, weight_str)
    return h_tmp

def get_faketype_ylds(h_fake):
    fake_dict = {}
    for faketype, enum in IFF_Type.items():
        ibin = enum + 1 # enum starts at 0 but bins start at 1
        yld = h_fake.GetBinContent(ibin)
        err = h_fake.GetBinError(ibin)
        result = UncFloat(yld, err)
        fake_dict[faketype] = result
    return fake_dict

def rank_ylds_str(yld_dict, tabs='', min_yld_perc = 0.0):
    rank_str = ''
    total = sum(yld_dict.values())
    if total == UncFloat():
        return ''
    other = UncFloat()
    for k, yld in sorted(yld_dict.items(), key=lambda (k,v): v, reverse=True):
        perc = yld / total
        if perc.value > min_yld_perc:
            rank_str += '%s%20s (%-15s)[%-15s%%]\n' % (tabs, k, yld, perc * 100)
        else:
            other += yld
    if other > UncFloat():
        rank_str += '%s%20s (%-15s)[%-15s%%]; ' % (tabs, 'Remainder', other, (other/total) * 100)
    return rank_str
    

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
        parser.add_argument("--type",
                                default='2lep',
                                help='event type (2lep, 3lep)')
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

