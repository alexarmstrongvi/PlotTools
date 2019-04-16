import sys
sys.path.append('../..') # Is this needed?

import collections # OrderedDict

class Cut :
    def __init__(self, cut, name = "", displayname = "", latexname = ""):
        if not name: 
            if len(cut) > 100: name = cut[:22]+' ... '+cut[-22:]
            else: name = cut
        if not displayname: displayname = name
        if not latexname: latexname = displayname
        self.name = name
        self.displayname = displayname
        self.latexname = latexname
        self.cut_str = cut

class Region :
    def __init__(self, name="", displayname="", latexname="", cuts=None) :
        
        if not displayname: displayname = name
        if not latexname: latexname = displayname
        
        self.name = name
        self.displayname = displayname
        self.latexname = latexname
        self.tcut = cuts
        self.compare_regions = []
        self.isSR = False

        self.yld_table = None # For linking to a YieldTbl object

        self.truth_fake_sel = None
        self.truth_bkg_sel = None


    def build_channel(self, name, displayname='', latexname='', cuts=''):

        if not displayname: displayname = name
        if not latexname: latexname = displayname
        
        name = "%s_%s" % (self.name, name)
        displayname = "%s (%s)" % (self.displayname, displayname)
        latexname = "%s (%s)" % (self.latexname, latexname)
        cuts = "%s && (%s)" % (self.tcut, cuts)

        new_reg = Region(name, displayname, latexname, cuts)
        new_reg.isSR = self.isSR

        return new_reg
    
    @property
    def cutflow(self):
        if not self._cutflow:
            cuts_list = split_cut_str(self.tcut)
            self._cutflow = build_cutflow(cuts_list)
        return self._cutflow
    
    def Print(self) :
        if self.is_cutflow :
            print 'Region "%s" ("%s") -- CutFlow -- %s'%(self.displayname, self.name, self.cutFlow[-1])
        else :
            print 'Region "%s" ("%s"): %s'%(self.displayname, self.name, self.tcut)

def split_cut_str(cut_str):
    # Split cut string into main AND'd parts
    # Keep together terms inside of parentheses
    cuts_list = []
    
    naive_split = [x.strip() for x in cut_str.split("&&")]
    parentheses_count = 0
    cut_parts = []
    for cut in naive_split:
        parentheses_count += cut.count("(") - cut.count(")")
        cut_parts.append(cut)
        if parentheses_count == 0:
            cut_class = Cut(' && '.join(cut_parts))
            cuts_list.append(cut_class)
            cut_parts = []
    return cuts_list    

def build_cutflow(cuts_list):
    # Build incremented cutflow cut list from list of individual cuts
    # Currently assumes list is filled with Cut class objects
    cutflow = [Cut("1","no_selection","No Selection")]    
    for ii, cut in enumerate(cuts_list):
        sel = ' && '.join([x.cut_str for x in cuts_list[:(ii+1)]])
        cutflow_cut = Cut(sel, cut.name, cut.displayname, cut.latexname)
        cutflow.append(cutflow_cut)
        
    return cutflow

