import sys
sys.path.append('../..')

import collections # OrderedDict

class Region :
    def __init__(self, name="", displayname="") :
        self.name = name
        self.displayname = displayname
        self.tcut = None
        self.truth_fake_sel = None
        self.truth_bkg_sel = None

        self.is_cutflow = False
        self.cutDict = collections.OrderedDict()
        self.cutFlow = []


    def setCutFlow(self) :
        '''
        Toggle whether to use this Region
        to perform a cutflow
        '''
        self.is_cutflow = True

    def isCutFlow(self) :
        '''
        Check whether the Region is configured
        as a cutflow
        '''
        return self.is_cutflow

    def getCutFlowDict(self) :
        return self.cutDict
    def getCutFlowList(self) :
        return self.cutFlow

    def addCut(self, name="", tcut="") :
        '''
        Method to add a cut and increment the cutflow
        '''
        self.cutDict[name] = tcut
        if len(self.cutFlow) == 0 :
            self.cutFlow.append(tcut)
        else :
            previous_cut_index = len(self.cutFlow) - 1
            cutflow_before_this_cut = self.cutFlow[previous_cut_index]
            cutflow_with_this_cut = cutflow_before_this_cut + " && " + tcut
            self.cutFlow.append(cutflow_with_this_cut)

    def printCutflow(self) :
        '''
        Method to print the cutflow for debugging
        '''
        print self.cutDict
        print self.cutFlow


    def Print(self) :
        if self.is_cutflow :
            print 'Region "%s" ("%s") -- CutFlow -- %s'%(self.displayname, self.name, self.cutFlow[-1])
        else :
            print 'Region "%s" ("%s"): %s'%(self.displayname, self.name, self.tcut)
