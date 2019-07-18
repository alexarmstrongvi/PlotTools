"""
================================================================================
Class for handling samples for plotting

Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with lots of ideas borrowed from Danny Antrim <dantrim@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""
# General python
import sys, os
import glob
import hashlib
import re
from copy import copy
# Root data analysis framework
import ROOT as r
import time
r.TColor.__init__._creates = False
r.TEventList.__init__._creates = False

################################################################################
# Main Base Class
################################################################################
class Sample :
    ''' '''
    # Static variables for when all inputs have the same property.
    # Can be updated for specific samples if needed
    input_file_treename = 'Unset'
    check_for_duplicates_flag = False
    weight_str = ''
    scale_factor = 1.0

    def __init__(self, name = "", displayname = "", latexname = ""):
        if not displayname: displayname = name
        if not latexname: latexname = displayname
        self.name = name
        self.displayname = displayname
        self.latexname = latexname
        self.file_path = ""
        self.color = r.kRed
        self.tree = None
        self.isMC = None
        self.cut = None

    # Setup TChain/TTree
    def set_chain_from_root_file(self, file_name, flat_ntuple_dir):
        self.file_path = os.path.realpath(flat_ntuple_dir)
        full_name = flat_ntuple_dir+file_name
        self.set_chain_from_list([full_name])

    def set_chain_from_dir(self, flat_ntuple_dir, save_dir='', search_str = "*"):
        '''
        Build chain from directory of ntuples 

        Args:
            flat_ntuple_dir (str): Name of dir containing all .root flat ntuples
            search_str (str): A globbing pattern to use in grabbing files from flat_ntuple_dir
        '''

        # Get list of flat ntuples file names from sample directory
        self.file_path = os.path.realpath(flat_ntuple_dir)
        print "Getting ntuples from", self.file_path
        flat_ntuples = glob.glob(flat_ntuple_dir + search_str)
        assert len(flat_ntuples), "No root files found at %s"%flat_ntuple_dir

        self.set_chain_from_list(flat_ntuples, save_dir)

    def set_chain_from_dsid_list(self, dsid_list, flat_ntuple_dir, checklist=[], search_strs=[], exclude_strs=[]):
        '''
        Build chain of background ntuples

        Args:
            dsid_list (list(int)): list of DSIDs to include in TChain
            flat_ntuple_dir (str): Name of dir containing all .root flat ntuples
            checklist (list(str)): substrings that must be found in the set of filenames
                matched to each DSID. Used for checking that all MC campaigns are found
            search_strs (list(str)): substrings that must be found in each filename
            exclude_strs (list(str)): substrings that must not be found in each filename

        Returns:
            TChain: TChain of flat ntuples from input directory
        '''
        # Get list of flat ntuples file names from sample directory
        self.file_path = os.path.realpath(flat_ntuple_dir)
        print "Getting ntuples from", self.file_path

        flat_ntuples = glob.glob(self.file_path + "/*.root")
        assert len(flat_ntuples), "No root files found at %s"%flat_ntuple_dir

        # Select out flat ntuples found in DSID list
        chosen_ntuples = []
        for dsid in dsid_list:
            matched_files = []
            for fname in flat_ntuples:
                # Apply filters and selections
                if any(s not in fname for s in search_strs if s): continue
                if any(s in fname for s in exclude_strs if s) : continue
                
                # Check for DSID match
                if str(dsid) in fname :
                    matched_files.append(fname)

            # Checks
            missed = filter(lambda x : not any(x in y for y in matched_files), checklist)
            if not matched_files:
                print "WARNING :: Unable to find any files for DSID =", dsid
            elif checklist and missed:
                print "WARNING :: Unable to find all files",
                print "for DSID = %s (missing: %s)" % (dsid, ", ".join(missed))

            # Store matched files
            chosen_ntuples += matched_files;

        # Checks
        if not chosen_ntuples:
            print "WARNING :: No samples found for", self.name

        self.set_chain_from_list(chosen_ntuples)

    def set_chain_from_file_list(self, file_list, flat_ntuple_dir):
        full_file_list = [flat_ntuple_dir+f for f in file_list]
        self.set_chain_from_list(full_file_list)

    def set_chain_from_list(self, files):
        chain = r.TChain(self.input_file_treename)
        for fname in files:
            chain.Add(fname)

            # Optional check for duplicates
            if self.check_for_duplicates_flag:
                f = r.TFile(fname)
                dup_events = self.check_for_duplicates(f)
                if dup_events:
                    print "Duplicate events found in", fname
                    print dup_events
                f.Close()
        self.tree = chain
        print "%10s : Chained %d files"%(self.name, len(files))
        
        # Code below was for a failed attempt at caching TChain as a single
        # file to avoid building TChain every time. Problem was the files get
        # too big
        #chain_final = r.TChain(self.input_file_treename)
        #if len(files) > 1 and save_dir:
        #    # Define variables
        #    n_files = str(len(files))
        #    totmtime = str(int(sum(map(os.path.getmtime, files))))
        #    full_names = map(lambda x : os.path.abspath(os.path.realpath(x)), files)
        #    cache_name = get_cache_name(full_names) + '.root'
        #    save_path = os.path.join(save_dir, cache_name)
        #    save_path = os.path.realpath(os.path.normpath(save_path))
        #    load_file = True
        #    if os.path.isfile(save_path):
        #        print "%10s : Sample file found at %s" % (self.name, save_path)
        #        rfile = r.TFile.Open(save_path) 
        #        try:
        #            stored_save_path = rfile.Get("save_path").GetTitle()
        #            stored_totmtime_path = rfile.Get("totmtime").GetTitle()
        #            stored_n_files = rfile.Get("n_files").GetTitle()
        #            if save_path != stored_save_path:
        #                print "Eventlist path has changed.",
        #                print "Playing it safe and remaking EventList."
        #                load_file = False
        #            elif totmtime != stored_totmtime_path:
        #                print "Sum of file modification times is different.",
        #                print "Playing it safe and remaking EventList."
        #                load_file = False
        #            elif n_files != stored_n_files:
        #                print "Number of entries in tree has changed. Remaking Eventlist"
        #                load_file = False
        #        except ReferenceError:
        #            print "TEventList not formatted as expected. Remaking Eventlist."
        #            load_file = False
        #        rfile.Close()
        #    else:
        #        load_file = False
        #     
        #    if not load_file:
        #        print "Creating combined flat ntuple for", self.name
        #        chain = r.TChain(self.input_file_treename)
        #        for fname in files:
        #            chain.Add(fname)

        #            # Optional check for duplicates
        #            if self.check_for_duplicates_flag:
        #                f = r.TFile(fname)
        #                dup_events = self.check_for_duplicates(f)
        #                if dup_events:
        #                    print "Duplicate events found in", fname
        #                    print dup_events
        #                f.Close()
        #        print "\tAdding %d files ..."%(len(files)),
        #        combined_file = r.TFile.Open(save_path, 'RECREATE')
        #        combined_file.cd()
        #        combined_tree = chain.CloneTree(-1,"fast")
        #        r.TNamed('save_path', save_path).Write()
        #        r.TNamed('totmtime', totmtime).Write()
        #        r.TNamed('n_files',n_files).Write()
        #        combined_file.Write()
        #        combined_file.Close()
        #        print "Done"
        #    
        #    self.file_path = save_path
        #elif len(files) == 1:
        #    self.file_path = files[0]
        #else:
        #    print "ERROR :: No files provided for sample"
        #    return
        #chain_final.Add(self.file_path)
        #self.tree = chain_final

    def set_event_list(self, cut, list_name, save_dir, reset=True):
        # Checks
        if not self.tree:
            print "WARNING :: no tree set for", self.name
            return

        # Reset event list
        if reset:
            self.tree.SetEventList(0)

        # Define useful variables
        # TODO: Remove TCut conversion from string
        cut = str(cut)
        if self.cut:
            cut = "%s && %s" % (cut, self.cut)
        tcut  = r.TCut(cut)
        n_entries = str(self.tree.GetEntries())
        identifiers = [cut, save_dir, self.file_path, str(n_entries), list_name]
        list_file_name = get_cache_name(identifiers) + '.root'
        save_path = os.path.join(save_dir, list_file_name)
        save_path = os.path.realpath(os.path.normpath(save_path))


        # Check if the list already exists
        load_eventlist = True
        if os.path.isfile(save_path) :
            print "%10s : EventList found at %s..."%(self.name, save_path),
            sys.stdout.flush()
            start = time.time()
            rfile = r.TFile.Open(save_path)

            # Check for any changes
            try:
                # Check that the expected variables are stored in the TEventList.
                # Above, TFile::Get will return a nullptr TObject if the variable
                # does not exist. So one must apply some attribute check (GetTitle)
                # to make sure there is no reference error that gets raised. If
                # the variables are correctly grabbed then proceed with other
                # checks
                stored_cut = rfile.Get("cut").GetTitle()
                stored_save_path = rfile.Get("save_path").GetTitle()
                stored_file_path = rfile.Get("file_path").GetTitle()
                stored_n_entries = rfile.Get("n_entries").GetTitle()

                if not rfile.GetListOfKeys().Contains(list_name):
                    print "\nNo TEventList found with name %s. Remaking EventList" % list_name
                    load_eventlist = False
                elif tcut != stored_cut:
                    print "\nEventList cuts have changed. Remaking EventList."
                    load_eventlist = False
                elif save_path != stored_save_path:
                    print "\nEventlist path has changed.",
                    print "Playing it safe and remaking EventList."
                    load_eventlist = False
                elif self.file_path != stored_file_path:
                    print "\nPath to sample files has changed. Remaking Eventlist."
                    load_eventlist = False
                elif n_entries != stored_n_entries:
                    print "\nNumber of entries in tree has changed. Remaking Eventlist"
                    load_eventlist = False
            except ReferenceError:
                print "\nTEventList not formatted as expected. Remaking Eventlist."
                load_eventlist = False
        else:
            load_eventlist = False

        # Load/Create evenlist
        if load_eventlist:
            event_list = rfile.Get(list_name)
            self.tree.SetEventList(event_list)
            end = time.time()
            dur = end - start
            print "LOADED (%.2fsec)" % dur
        else:
            print "Creating TEventList for %s..." % self.name,
            sys.stdout.flush()
            start = time.time()

            rfile = r.TFile(save_path,'recreate')
            draw_list = ">> " + list_name
            self.tree.Draw(draw_list, tcut)
            event_list = r.gROOT.FindObject(list_name).Clone()
            self.tree.SetEventList(event_list)
            event_list.Write(list_name)

            # Append other information
            tcut.Write("cut")
            r.TNamed("save_path",save_path).Write()
            r.TNamed("file_path",self.file_path).Write()
            r.TNamed("n_entries",n_entries).Write()
            
            end = time.time()
            dur = end - start
            print "Done (%dmin %.2fsec)" % (int(dur)/60, dur%60) 

        rfile.Close()

    # Comparison
    def __eq__(self, other) :
        return (isinstance(other, Sample)
            and self.displayname == other.displayname
            and self.name == other.name)

     # Comparison
    def Print(self) :
        print 'Sample "%s" (tree %s)'%(self.displayname, self.name)

    def __gt__(self, other, tcut) :
        '''
        Comparison operator to order background samples by
        their yields in a given region defined by tcut
        '''
        cut = r.TCut(tcut)
        return ( self.tree.Draw("isMC", cut, "goff") > other.tree.Draw("isMC", cut, "goff") )

    # Sanity checks
    def is_setup(self):
        return self.tree and self.tree.GetEntries() and self.isMC != None

    def check_for_duplicates(self, ifile):
        tree = ifile.Get(self.input_file_treename)
        events = [x.event_number for x in tree]
        
        # Get duplicates events
        seen = set()
        dup_evts = set()
        for x in events:
            if x not in seen: seen.add(x)
            else: dup_evts.add(x)

        if dup_evts:
            print "There are %d/%d duplicate events"%(len(dup_evts), len(events))
        return dup_evts

################################################################################
# Data class
################################################################################
class Data(Sample):
    def __init__(self, name = 'data', displayname='Data', latexname='Data') :
        Sample.__init__(self, name, displayname, latexname)
        self.weight_str = '1'
        self.scale_factor = '1'
        self.color = r.kBlack
        self.isMC = False
        self.isDataBkg = False
        self.blinded = True

    def Print(self) :
        print 'Data (tree %s)'%(self.name)

class DataBackground(Data) :
    def __init__(self, name = "", displayname = "", latexname="") :
        Data.__init__(self, name, displayname)
        self.isDataBkg = True

    def Print(self) :
        print 'Data background "%s" (tree %s)'%(self.displayname, self.name)
################################################################################
# MC classes
################################################################################
class MCsample(Sample):

    def __init__(self, name = "", displayname = "", latexname=""):
        Sample.__init__(self, name, displayname, latexname)
        self.dsid = ""
        self.line_style = 1
        self.fill_style = 0
        self.isSignal = None
        self.isMC = True

    def isSignal(self) :
        return self.isSignal
    def Print(self) :
        print 'MC Sample "%s" (tree %s)'%(self.displayname, self.name)

class MCBackground(MCsample) :
    def __init__(self, name = "", displayname = "", latexname="") :
        MCsample.__init__(self, name, displayname, latexname)
        self.isSignal = False

    def Print(self) :
        print 'MC Background "%s" (tree %s)'%(self.displayname, self.name)

class Signal(MCsample) :
    def __init__(self, name = "", displayname ="", latexname="") :
        MCsample.__init__(self, name, displayname, latexname)
        self.isSignal = True

    def Print(self) :
        print 'Signal "%s" (tree %s)'%(self.displayname, self.name)

def get_cache_name(str_identifiers):
    hash_input = '_'.join(str_identifiers)
    hash_object = hashlib.md5(b'%s' % hash_input)
    return hash_object.hexdigest()

color_palette = {
    # Preferred ordering for colors when overlaying/stacking multiple plots
    # Colors should all be visually distinguishable and look good :)
    'gray'    : [r.kGray    +0, r.kGray    +1, r.kGray    +2, r.kGray    +3, r.kBlack   +0],
    'red'     : [r.kRed     +1, r.kRed     -5, r.kRed     +3, r.kRed     -7, r.kRed     -2],
    'orange'  : [r.kOrange  +1, r.kOrange  -5, r.kOrange  +3, r.kOrange  -7, r.kOrange  -2],
    'yellow'  : [r.kYellow  +1, r.kYellow  -5, r.kYellow  +3, r.kYellow  -7, r.kYellow  -2],
    'spring'  : [r.kSpring  +1, r.kSpring  -5, r.kSpring  +3, r.kSpring  -7, r.kSpring  -2],
    'green'   : [r.kGreen   +1, r.kGreen   -5, r.kGreen   +3, r.kGreen   -7, r.kGreen   -2],
    'teal'    : [r.kTeal    +1, r.kTeal    -5, r.kTeal    +3, r.kTeal    -7, r.kTeal    -2],
    'cyan'    : [r.kCyan    +1, r.kCyan    -5, r.kCyan    +3, r.kCyan    -7, r.kCyan    -2],
    'azure'   : [r.kAzure   +1, r.kAzure   -5, r.kAzure   +3, r.kAzure   -7, r.kAzure   -2],
    'blue'    : [r.kBlue    +1, r.kBlue    -5, r.kBlue    +3, r.kBlue    -7, r.kBlue    -2],
    'violet'  : [r.kViolet  +1, r.kViolet  -5, r.kViolet  +3, r.kViolet  -7, r.kViolet  -2],
    'magenta' : [r.kMagenta +1, r.kMagenta -5, r.kMagenta +3, r.kMagenta -7, r.kMagenta -2],
    'pink'    : [r.kPink    +1, r.kPink    -5, r.kPink    +3, r.kPink    -7, r.kPink    -2],
}
