#!/usr/bin/env python
"""
================================================================================
Compare two flat ntuples to check for differences

Examples
    python compare_flat_ntuples.py file1.root file2.root
    python compare_flat_ntuples.py file1.root file2.root -t tree_name
    python compare_flat_ntuples.py file1.root file2.root -s -d ./Plots -p Test
    python compare_flat_ntuples.py file1.root file2.root -r -d ./Files -o name.root

Compares the TH1 histograms or flat (i.e. fully split) TTree branches in two 
TFile objects. All other objects found with TFile::GetListOfKeys() are ignored.
The outputs are printouts of the number of shared and unique keys as well as the
number of shared keys for TH1/TBranch objects that are identical and different.

Default is to run on top level TH1 objects but providing a ttree name will
cause the script to compare branches in the ttree instead.

For compared objects that differ, an option can be set to save a comparison 
plot as pdfs or in a root file. 

Set verbose mode to print out names of histograms/branches falling into the 
various categories.

Author:
    Alex Armstrong <alarmstr@cern.ch>
    December 9, 2018
================================================================================
"""

# TODO: Add ability to make multiple plots for vector branches.
#       Keep making plots for each index until plots are empty

# Imports
import sys, os, traceback, argparse
import time
import subprocess

from ROOT import TH1F, TH1, TFile
from ROOT import kBlue, kRed, kWarning
import ROOT  # needed for ROOT globals
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = kWarning
ROOT.TPad.__init__._creates = False

from plot_utils import print_hist, get_branch_max, get_branch_min

# User Argument defaults and help information
_help_ifile_name1 = 'First input file name'
_help_ifile_name2 = 'Second input file name'
_help_save_cmp_plots = 'Save comparison plots of differing histograms as separate pdf images'
_help_tree_name = 'Compare branches of a TTree instead of top level histograms'
_df_output_dir = './'
_help_output_dir = 'Directory for outputs if requested [default: %s]' % _df_output_dir 
_help_ofile_name = 'Root file name for storing comparison plots instead of saving individual plots [default: flatNt_cmp_file1_file2.root]'
_help_plot_name_suffix = 'Suffix for all comparison plot save names'
_df_verbose = False #TODO: Test if this is needed. If not update skeleton
_help_verbose = 'verbose output'


################################################################################
def main ():
    """ Main Function """

    global args
    check_environment()
    check_inputs(args)

    ########################################################################### 
    # Do the magic
    f1 = TFile(args.ifile_name1, "READ") 
    f2 = TFile(args.ifile_name2, "READ") 

    if args.tree_name:
        # Check TTree exists in both files before accessing
        found_key = f1.GetListOfKeys().Contains(args.tree_name)
        found_key = found_key and f2.GetListOfKeys().Contains(args.tree_name)
        if not found_key:
            print "ERROR :: TTree (%s) name not found in both TFiles: %s and %s" % (
                    tree_name, file1.GetName(), file2.GetName())
            sys.exit()
        t1 = f1.Get(args.tree_name)
        t2 = f2.Get(args.tree_name)

        # Do main comparison
        shared_keys, unique1_keys, unique2_keys = get_ttree_branch_sets(t1, t2)
        identical_keys, diff_keys = compare_flat_ttrees(t1, t2, shared_keys)
    else:
        shared_keys, unique1_keys, unique2_keys = get_tfile_key_sets(f1, f2, args.tree_name)
        identical_keys, diff_keys = compare_flat_ntuples(f1, f2, shared_keys)
    

    ########################################################################### 
    # Print results
    n_shared = len(shared_keys)
    n_identical = len(identical_keys)
    n_diff = len(diff_keys)
    n_uniq1 = len(unique1_keys)
    n_uniq2 = len(unique2_keys)
    print "\n===== Results ====="
    print "Total keys: %d" % (n_shared + n_uniq1 + n_uniq2)
    print "\tShared keys: %d" % n_shared
    print "\t\tDifferent: %d" % n_diff
    print "\t\tIdentical: %d" % n_identical
    print "\tUnique to %s: %d" % (args.ifile_name1, n_uniq1)
    print "\tUnique to %s: %d" % (args.ifile_name2, n_uniq2)
    
    if args.verbose:
        print "\n===== Printout of keys ===== \n"
        print "\nDifferent: "
        for key in sorted(diff_keys): print "\t", key
        print "\nIdentical: "
        for key in sorted(identical_keys): print "\t", key
        print "\nUnique to %s: " % args.ifile_name1
        for key in sorted(unique1_keys): print "\t", key
        print "\nUnique to %s: " % args.ifile_name2
        for key in sorted(unique2_keys): print "\t", key

    ########################################################################### 
    # Make plots if requested
    if args.tree_name and (args.save_plots or args.save_root):
        flat_ttree_compare_plots(args.ifile_name1, args.ifile_name2, 
                                 args.tree_name,
                                 diff_keys,
                                 args.out_dir,
                                 args.save_root,
                                 args.ofile_name,
                                 args.plot_suffix)   
    elif not args.tree_name and (args.save_plots or args.save_root):
        flatnt_compare_plots(f1, f2, diff_keys, args.out_dir, args.ofile_name, args.plot_suffix)   


################################################################################
# FUNCTIONS
def check_inputs(args):
    """ Check the input arguments are as expected """

    if args.ifile_name1 and not os.path.exists(args.ifile_name1):
        print "ERROR :: Cannot find input file:", args.ifile_name1
        sys.exit()
    if args.ifile_name2 and not os.path.exists(args.ifile_name2):
        print "ERROR :: Cannot find input file:", args.ifile_name2
        sys.exit()
    
    # Check for inconsistent options
    if args.plot_suffix and not args.save_plots:
        print "ERROR :: Plot suffix provided but saving plot option not requested"
        sys.exit()
    if args.ofile_name and not args.save_root:
        print "ERROR :: Output root name provided but saving root option not requested"
        sys.exit()
    if args.out_dir and not (args.save_plots or args.save_root):
        print "ERROR :: Output directory provided but no save option requested (i.e. root or plots)"
        sys.exit()

    # Set default directory
    if not args.out_dir:
        args.out_dir = _df_output_dir

    # Build default output file name if needed
    if not args.ofile_name:
        #TODO: Test this
        name1 = args.ifile_name1.replace(".root","")
        name2 = args.ifile_name2.replace(".root","")
        args.ofile_name = 'flatNt_cmp_%s_%s.root' % (name1, name2 )

    # Check if output file exists
    # If so, check with user if they want to overwrite it
    if os.path.exists(args.ofile_name):
        #TODO: Test this then add to skeleton
        usr_msg =  "Output file already exists: %s\n" % args.ofile_name
        usr_msg += "Would you like to overwrite it? [Y/N] "
        overwrite_op = raw_input(usr_msg)
        
        # Only accept Y or N 
        while overwrite_op not in ["Y","N"]:
            usr_msg = "Unacceptable answer: %s" % overwrite_op
            usr_msg += "Would you like to overwrite it? [Y/N] "
            overwrite_op = raw_input(usr_msg)

        if overwrite_op == "N":
            sys.exit()
    
def check_environment():
    """ Check if the shell environment is setup as expected """
    assert os.environ['USER'], "USER variable not set"

    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    assert python_ver >= 2.7, ("Running old version of python\n", sys.version)



def get_tfile_key_sets(file1, file2, tree_name=''):
    """
    Get the shared and unique key sets for two TFiles
    args:
        file1 (ROOT.TFile) - First input root file
        file2 (ROOT.TFile) - Second input root file
    returns:
        (tuple of 3 sets):
        set 1: Keys shared by both TFiles
        set 2: Keys unique to the first TFile
        set 3: Keys unique to the second TFile
    """
    f1_keys = set([k.GetName() for k in file1.GetListOfKeys()])
    f2_keys = set([k.GetName() for k in file2.GetListOfKeys()])
    
    shared_keys = f1_keys & f2_keys
    unique1_keys = f1_keys - f2_keys
    unique2_keys = f2_keys - f1_keys

    return shared_keys, unique1_keys, unique2_keys

def get_ttree_branch_sets(ttree1, ttree2):
    """
    Get the shared and unique key sets for two TFiles
    args:
        ttree1 (ROOT.TTree) - First input root TTree
        ttree2 (ROOT.TTree) - Second input root TTree
    returns:
        (tuple of sets(str)):
        set 1: Keys shared by both TTrees
        set 2: Keys unique to the first TTree
        set 3: Keys unique to the second TTree
    """
    t1_branches = set([k.GetName() for k in ttree1.GetListOfBranches()])
    t2_branches = set([k.GetName() for k in ttree2.GetListOfBranches()])
    
    shared_branches = t1_branches & t2_branches
    unique1_branches = t1_branches - t2_branches
    unique2_branches = t2_branches - t1_branches

    return shared_branches, unique1_branches, unique2_branches

def compare_flat_ntuples(file1, file2, keys, tree_name=''):
    """
    Determine identical and different elements between two files 
    args:
        file1 (ROOT.TFile) - First input root file
        file2 (ROOT.TFile) - Second input root file
        keys (list(str)) - list of histogram keys to compare
    returns:
        (tuple of sets(str)):
        set 1: Keys for identical histograms
        set 2: Keys for differing histograms
    """
    not_found_keys = []
    skipped_keys = []
    identical_keys = set()
    diff_keys = set()
    for key in keys:
        # Check that key exists in both files
        found_key = file1.GetListOfKeys().Contains(key)
        found_key = found_key and file2.GetListOfKeys().Contains(key)
        if not found_key:
            not_found_keys.append(key_check)
            continue
        
        # Make sure key is for TH1 
        hist1 = file1.Get(key)
        hist2 = file2.Get(key)
        are_hists = isinstance(hist1, TH1) and isinstance(hist2, TH1)
        if not are_hists: 
            skipped_keys.append(key)
            continue
        
        # Determine if there is a match
        if hists_are_diff(hist1, hist2): 
            diff_keys.add(key)
        else: 
            identical_keys.add(key)

    if not_found_keys:
        print "WARNING :: %d requested keys were not found:\n\t" % len(not_found_keys), not_found_keys
    if skipped_keys:
        print "WARNING :: %d requested keys were not histograms:\n\t" % len(skipped_keys), skipped_keys

    return identical_keys, diff_keys

def compare_flat_ttrees(tree1, tree2, branch_names):
    """
    Determine identical and different branches between two flat TTrees 
    args:
        tree1 (ROOT.TTree) - First input root tree
        tree2 (ROOT.TTree) - Second input root tree
        branch_names (list(str)) - list of branch names to compare
    returns:
        (tuple of sets(str)):
        set 1: Keys for identical branches
        set 2: Keys for differing branches
    """
    not_found_br = []
    identical_br = set()
    diff_br = set()
    for br in branch_names:
        are_branches = tree1.GetListOfBranches().Contains(br)
        are_branches = are_branches and tree2.GetListOfBranches().Contains(br)
        if not are_branches:
            not_found_br.append(br)
            continue

        # Determine if there is a match
        hist1, hist2 = make_ttree_compare_hists(tree1, tree2, br)
        if not hist1 or not hist2: continue
    
        if hists_are_diff(hist1, hist2): 
            diff_br.add(br)
        else: 
            identical_br.add(br)

    if not_found_br:
        print "WARNING :: %d requested branches were not found:\n\t" % len(not_found_br), not_found_br
    return identical_br, diff_br 

def make_ttree_compare_hists(tree1, tree2, branch_name):
    """
    Make comparable (i.e. same range and nbins) histograms of a branch 
    found in two trees. 
    In cases where it is not possible to draw a histogram, None is 
    returned for both hists
    args:
        tree1 (ROOT.TTree) - first Root tree
        tree2 (ROOT.TTree) - second Root tree
        branch_name (str) - branch to be plotted
    return:
        (ROOT.TH1F) - histogram from first Root tree
        (ROOT.TH1F) - histogram from second Root tree
    """
    nbins = 100
    min1 = get_branch_min(tree1, branch_name)
    min2 = get_branch_min(tree2, branch_name)
    final_min = min(min1, min2)

    max1 = get_branch_max(tree1, branch_name)
    max2 = get_branch_max(tree2, branch_name)
    final_max = max(max1, max2)

    val_range = final_max - final_min
    if val_range > 10e12:
        print "WARNING :: Branch %s has very large range [%s, %s]" % (
                branch_name, final_min, final_max)
    if val_range == 0:
        final_min -= 1
        final_max += 1
    if ROOT.gDirectory.FindObject("h1"):
        ROOT.gDirectory.Get("h1").Delete()
    h1 = TH1F("h1","",nbins,final_min, final_max)
    draw_cmd = "%s >> %s" % (branch_name, h1.GetName())
    tree1.Draw(draw_cmd,"","goff")

    if ROOT.gDirectory.FindObject("h2"):
        ROOT.gDirectory.Get("h2").Delete()
    h2 = TH1F("h2","",nbins,final_min, final_max)
    draw_cmd = "%s >> %s" % (branch_name, h2.GetName())
    tree2.Draw(draw_cmd,"","goff")

    return h1, h2

def flatnt_compare_plots(file1, file2, keys, out_dir="./", ofile_name="", suffix=""):
    """
    Make comparison plots (i.e. overlay with ratio) for two files
    args:
        file1 (ROOT.TFile) - First input root file
        file2 (ROOT.TFile) - Second input root file
        keys (list(str)) - list of keys to plot
        out_dir (str) - directory for saving output
        ofile_name (str) - optional root file name for storing outputs in
        suffix (str) - suffix for each plot name
    """
    print "ERROR :: Not yet setup to make comparison plots from histograms"
    return

def flat_ttree_compare_plots(file1_path, file2_path, tree_name, branch_names, out_dir="./", save_as_root=False, ofile_name="", suffix=""):
    """
    Make comparison plots (i.e. overlay with ratio) for two files
    args:
        file1_path (str) - First input root file name
        file2_path (str) - Second input root file name
        tree_name (str) - name of TTree stored in both root files
        branch_names (list(str)) - list of branch names to plot
        out_dir (str) - directory for saving output
        ofile_name (str) - optional root file name for storing outputs in a root file
        suffix (str) - suffix for each plot name
    """
    if save_as_root:
        print "ERROR :: Not yet setup to produce root files"
        return

    # Set samples to be processed
    from sample import Sample
    Sample.input_file_treename = tree_name
    samples = []
    
    file1_name = file1_path.split('/')[-1]
    file1_dir = file1_path.replace(file1_name,"") 
    sample_name = file1_name.replace(".root","")
    new_or_old = "New" if "new" in sample_name else "Old"
    sample1 = Sample(sample_name, new_or_old)
    sample1.color = kBlue
    sample1.set_chain_from_root_file(file1_name, file1_dir)
    samples.append(sample1)

    file2_name = file2_path.split('/')[-1]
    file2_dir = file2_path.replace(file2_name,"") 
    sample_name = file2_name.replace(".root","")
    new_or_old = "New" if "new" in sample_name else "Old"
    sample2 = Sample(sample_name, new_or_old)
    sample2.color = kRed
    sample2.set_chain_from_root_file(file2_name, file2_dir)
    samples.append(sample2)
    
    # Make region
    from region import Region
    region = Region("no_sel","No Selection")
    region.tcut = "1"

    # Build list of plots
    from plot import PlotBase, Plot1D, Types
    PlotBase.output_format = 'pdf'
    PlotBase.save_dir = out_dir
    Plot1D.doLogY = False
    Plot1D.auto_set_ylimits = True
    nbins = 25

    plots = []
    for br in branch_names:
        flag = False
        max_val = max(get_branch_max(s.tree, br) for s in samples)
        min_val = min(get_branch_min(s.tree, br) for s in samples)
        val_range = max_val - min_val
        if val_range == 0:
            max_val += 1
            min_val -= 1
            val_range = max_val - min_val
        range_min = min_val - 0.05*val_range
        range_max = max_val + 0.05*val_range
        plot = Plot1D(region=region.name, 
                      variable=br, 
                      bin_range=[range_min, range_max], 
                      nbins=nbins, 
                      xlabel=br, 
                      ptype=Types.ratio,
                      suffix=suffix)
        #if plot.variable == "METPhi":
        #    import pdb; pdb.set_trace()
        plot.setRatioPads(plot.name)
        plots.append(plot)
    
    # Loop over each plot and save image
    from hist import SampleCompare1D, RatioHist1D 
    for plot in plots:
        with SampleCompare1D(plot, region, samples) as hists: 
            num = hists.hists[0]
            den = hists.hists[1]
            if num.Integral() == 0 and den.Integral() == 0:
                print "INFO :: Skipping %s, all histograms are empty" % plot.variable
            with RatioHist1D(plot, num, den, ymax = 2, ymin = 0) as ratio_hist:
                ratio_label = "%s / %s" % (samples[0].displayname, samples[1].displayname)
                plot.make_overlay_with_ratio_plot(region.displayname, ratio_label, hists, ratio_hist) 


def hists_are_diff(hist1, hist2):
    ''' Compares string versions of hists to check if equal '''
    return (print_hist(hist1) != print_hist(hist2))
################################################################################

def get_args():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('ifile_name1',
                        help=_help_ifile_name1)
    parser.add_argument('ifile_name2',
                        help=_help_ifile_name2)
    parser.add_argument('-t', '--tree_name',
                        help=_help_tree_name)
    parser.add_argument('-s', '--save_plots',
                        action='store_true',
                        help=_help_save_cmp_plots)
    parser.add_argument('-d', '--out_dir',
                        help=_help_output_dir)
    parser.add_argument('-r', '--save_root',
                        action='store_true',
                        help=_help_save_cmp_plots)
    parser.add_argument('-o', '--ofile_name',
                        help=_help_ofile_name)
    parser.add_argument('-p', '--plot_suffix',
                        help=_help_plot_name_suffix)
    parser.add_argument('-v', '--verbose',
                        action='store_true', default=_df_verbose,
                        help=_help_verbose)
    args = parser.parse_args()
    return args

################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        # TODO: Add ability to check standard input so things can be piped
        args = get_args()
        if args.verbose:
            print '>'*40
            print 'Running {}...'.format(os.path.basename(__file__))
            print time.asctime()
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
