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

# Root data analysis framework
import ROOT as r
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

    # Setup TChain/TTree
    def set_chain_from_root_file(self, file_name, flat_ntuple_dir):
        self.file_path = os.path.realpath(flat_ntuple_dir)
        full_name = flat_ntuple_dir+file_name
        self.set_chain_from_list([full_name])

    def set_chain_from_dsid_list(self, dsid_list, flat_ntuple_dir, search_strs=[], exclude_strs=[]):
        '''
        Build chain of background ntuples

        Args:
            dsid_list (list(int)): list of DSIDs to include in TChain
            flat_ntuple_dir (str): Name of dir containing all .root flat ntuples

        Returns:
            TChain: TChain of flat ntuples from input directory
        '''
        # Get list of flat ntuples file names from sample directory
        self.file_path = os.path.realpath(flat_ntuple_dir)
        print "Getting ntuples from", self.file_path

        flat_ntuples = glob.glob(flat_ntuple_dir + "*.root")
        assert len(flat_ntuples), "No root files found at %s"%flat_ntuple_dir

        # Select out flat ntuples found in DSID list
        chosen_ntuples = []
        for dsid in dsid_list:
            for fname in flat_ntuples:
                if any(s not in fname for s in search_strs if s): continue
                if any(s in fname for s in exclude_strs if s) : continue
                if str(dsid) in fname :
                    chosen_ntuples.append(fname);
                    break
            else:
                print "WARNING :: Unable to find file for DSID =", dsid
        if not chosen_ntuples:
            print "WARNING :: No samples found for", self.name

        self.set_chain_from_list(chosen_ntuples)

    def set_chain_from_file_list(self, file_list, flat_ntuple_dir):
        self.file_path = os.path.realpath(flat_ntuple_dir)
        full_file_list = [flat_ntuple_dir+f for f in file_list]
        self.set_chain_from_list(full_file_list)

    def set_chain_from_list(self, files):
        chain = r.TChain(self.input_file_treename)
        for n_files, fname in enumerate(files):
            chain.Add(fname)

            # Optional check for duplicates
            if self.check_for_duplicates_flag:
                f = r.TFile(fname)
                dup_events = self.check_for_duplicates(f)
                if dup_events:
                    print "Duplicate events found in", fname
                    print dup_events
                f.Close()

        print "%10s : ADDED %d FILES"%(self.name, n_files+1)
        self.tree = chain

    def set_event_list(self, cut, list_name, save_dir):
        # Checks
        if not self.tree:
            print "WARNING :: no tree set for", self.name
            return

        # Reset event list
        self.tree.SetEventList(0)

        # Define useful variables
        cut = str(cut)
        tcut  = r.TCut(cut)
        n_entries = str(self.tree.GetEntries())
        hash_input = '_'.join([cut, save_dir, self.file_path, str(n_entries)])
        hash_object = hashlib.md5(b'%s' % hash_input)
        list_name = hash_object.hexdigest()
        list_file_name = list_name + ".root"
        save_path = os.path.join(save_dir, list_file_name)
        save_path = os.path.realpath(os.path.normpath(save_path))


        # Check if the list already exists
        load_eventlist = True
        if os.path.isfile(save_path) :
            rfile = r.TFile.Open(save_path)
            print "%10s : EventList found at %s"%(self.name, save_path)

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

                if tcut != stored_cut:
                    print "EventList cuts have changed. Remaking EventList."
                    load_eventlist = False
                elif save_path != stored_save_path:
                    print "Eventlist path has changed.",
                    print "Playing it safe and remaking EventList."
                    load_eventlist = False
                elif self.file_path != stored_file_path:
                    print "Path to sample files has changed. Remaking Eventlist."
                    load_eventlist = False
                elif n_entries != stored_n_entries:
                    print "Number of entries in tree has changed. Remaking Eventlist"
                    load_eventlist = False
            except ReferenceError:
                print "TEventList not formatted as expected. Remaking Eventlist."
                load_eventlist = False
        else:
            load_eventlist = False

        # Load/Create evenlist
        if load_eventlist:
            event_list = rfile.Get(list_name)
            self.tree.SetEventList(event_list)
        else:
            print "Creating TEventList for", self.name
            rfile = r.TFile(save_path,'recreate')
            draw_list = ">> " + list_name
            self.tree.Draw(draw_list, tcut)
            event_list = r.gROOT.FindObject(list_name).Clone()
            self.tree.SetEventList(event_list)
            event_list.Write(list_name)

            # Append other information
            tcut.Write("cut")
            r.TNamed("save_path",save_path).Write()
            #new_save_path = r.TNamed("save_path",save_path)
            #new_save_path.Write()
            new_file_path = r.TNamed("file_path",self.file_path)
            new_file_path.Write()
            new_n_entries = r.TNamed("n_entries",n_entries)
            new_n_entries.Write()
        
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
        sel = r.TCut("1")
        return ( self.tree.Draw("isMC", cut * sel, "goff") > other.tree.Draw("isMC", cut * sel, "goff") )

    # Sanity checks
    def is_setup(self):
        return self.tree and self.isMC != None

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
    def __init__(self, name = 'data', displayname='Data') :
        Sample.__init__(self, name, displayname)
        self.color = r.kBlack
        self.isMC = False

    def Print(self) :
        print 'Data (tree %s)'%(self.name)

################################################################################
# MC classes
################################################################################
class MCsample(Sample):

    # Static variables for when all inputs have the same property.
    # Can be updated for specific samples if needed
    weight_str = ''
    scale_factor = 1.0

    def __init__(self, name = "", displayname = "", latexname=""):
        Sample.__init__(self, name, displayname, latexname)
        self.dsid = ""
        self.line_style = 1
        self.fill_style = 0
        self.isSignal = None
        self.isMC = True

    def isSignal(self) :
        return self.is_signal
    def Print(self) :
        print 'MC Sample "%s" (tree %s)'%(self.displayname, self.name)

class Background(MCsample) :
    def __init__(self, name = "", displayname = "", latexname="") :
        MCsample.__init__(self, name, displayname, latexname)
        self.is_signal = False

    def Print(self) :
        print 'Background "%s" (tree %s)'%(self.displayname, self.name)

class Signal(MCsample) :
    def __init__(self, name = "", displayname ="", latexname="") :
        MCsample.__init__(self, name, displayname, latexname)
        self.is_signal = True

    def Print(self) :
        print 'Signal "%s" (tree %s)'%(self.displayname, self.name)


