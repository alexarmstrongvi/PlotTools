################################################################################
# Classes for making all ROOT objects needed for painting plots 
# such as TH1, TLegend, TGraph, etc..
#
# These get passed to the Plot classes to be turned into pretty plots
################################################################################
import re
from array import array

# Root data analysis framework
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)
r.RooStats.NumberCountingUtils
import ROOT.RooStats.NumberCountingUtils as ncu
# Turn off root ownership when creating classes
r.TEventList.__init__._creates = False
r.TH1F.__init__._creates = False
r.TGraphErrors.__init__._creates = False
r.TGraphAsymmErrors.__init__._creates = False
r.TLegend.__init__._creates = False
r.THStack.__init__._creates = False

from PlotTools.plot import Plot1D, Plot2D, Types
from PlotTools.region import Region
from PlotTools.YieldTable import UncFloat
import PlotTools.plot_utils as pu

#TODOs
# Move style stuff to ATLAS style files or plot.py and import

class HistBase:
    '''
    All histogram classes should be created within 'with' statments so that
    the owned root objects get deleted at the end
    '''

    def __init__(self):
        pass
    def __enter__(self):
        return self
    def __exit__(self, exception_type, exception_value, traceback):
        self.clear(self.list_of_root_objects())
    def list_of_root_objects(self):
        " Return list of all class members that are or contain root objects"
        #TODO: Get python to get all class variables that inherit fromTObject
        #      (e.g. use self.__class__.__dict__)
        raise NotImplementedError()

    def clear(self, list_of_root_objects):
        # Delete all underlying root objects
        for ro in list_of_root_objects:
            # Variables
            if not ro or isinstance(ro, r.THStack):
                # THStack is just a container for the hists in the stack and
                # therefore does not need to be deleted
                continue
            elif isinstance(ro, list):
                self.clear(ro)
            elif issubclass(type(ro), r.TObject):
                ro.Delete()
            else:
                print "ERROR :: Unknown primitive type", type(ro)

#TODO: Remover unused classes

class StackHist1D(HistBase):
    def __init__(self, plot, reg, samples) : 

        self.leg                    = r.TLegend()
        self.axis                   = self.make_plot1D_axis(plot)
        self.stack_hists            = make_stack_hists(plot, reg, samples)
        self.stack_hist             = make_stack_hist(self.stack_hists)
        self.total_hist             = make_stack_total_hist(self.mc_stack)
        self.total_error_graph      = make_error_graph(self.total_hist)
        self.total_error_dummy_hist = make_dummy_hist_for_legend("mcError") 
       
        if self.successful_setup() and plot.doNorm:
            self.normalize_hists(plot)

    def list_of_root_objects(self):
        return [
            self.leg,
            self.axis,
            self.stack_hists,
            self.stack_hist,
            self.total_hist,
            self.total_error_graph,
            self.total_error_dummy_hist,
    ]

    def successful_setup(self):
        return bool(self.stack_hists)

    def make_plot1D_axis(self, plot, name='axis'):
        return r.TH1F(name, "", plot.nbins, plot.bin_edges)

    def normalize_hists(self):
        norm_factor = 1.0/self.total_hist.Integral()
        pu.scale_thstack(self.stack_hist, norm_factor)
        self.mc_total.Scale(norm_factor)
        pu.scale_tgraph(self.total_error_graph, norm_factor)

def make_stack_hists(plot, reg, samples, sort=True):
    stack_hists = []
    for sample in samples:
        # Initilize histogram
        h = make_hist1D(plot, reg, sample)
        
        # Basic formatting
        h.SetFillColor(sample.color)
        h.SetFillStyle(1001)
        h.leg_name = sample.displayname

        if h.Integral() > 0:
            stack_hists.append(h)

    if not len(stack_hists):
        print "ERROR :: All hists are empty. Unable to make stack hist"
        return

    if sort:
        stack_hists = sorted(stack_hists, key=lambda h: h.Integral())

    return stack_hists

def make_stack_hist(stack_hists, name="h_stack"):
    h_stack = r.THStack(name)
    for h in stack_hists:
        h_stack.Add(h)
    return h_stack

def make_stack_total_hist(stack_hist, name="h_stack_total"):
    try:
        return stack_hist.GetStack().Last().Clone(name)
    except ReferenceError:
        print "ERROR :: Unable to access stack to get total"
        return

def make_error_graph(hist, symmetrize_errors=True):
    if not hist:
        print "ERROR :: Unable to make error graph from empty hists"
        return

    error_graph = pu.th1_to_tgraph(hist)

    if symmetrize_errors:
        for i in xrange(error_graph.GetN()) :
            ehigh = error_graph.GetErrorYhigh(i)
            elow  = error_graph.GetErrorYlow(i)

            error_sym = r.Double(0.0)
            error_sym = (ehigh + elow) / 2.0

            if ehigh != error_sym:
                print "INFO :: initial error (+%.2f,-%.2f), symmetrized = (+%.2f,-%.2f)"%(ehigh,elow, error_sym, error_sym)
                error_graph.SetPointEYhigh(i,0.0)
                error_graph.SetPointEYhigh(i, error_sym)
                error_graph.SetPointEYlow(i,0.0)
                error_graph.SetPointEYlow(i,error_sym)

    return error_graph

def make_dummy_hist_for_legend(name):
    return r.TH1F(name,"",1,0,1)


class CutScan1D(HistBase):
    def __init__(self, plot, reg, signals, backgrounds, xcut_is_max=True):
        '''
        Params:
            - xcut_is_max (bool) - Indicates cuts are of the form cut_var < cut_value
        '''
        self.axis = None
        self.hist_leg = None #ToDo
        self.signal = None
        self.bkgd = None
        self.eff_leg = None #ToDo
        self.signal_eff = None
        self.bkgd_rej = None
        self.s_over_b = None
        self.zn_sig = None
        self.roc_graph = None

        self.axis = make_plot1D_axis(plot)
        
        self.signal = make_hist1D(plot, reg, signals)
        self.bkgd = make_hist1D(plot, reg, backgrounds)
        
        self.hist_leg = make_basic_legend(plot)
        if len(signals) > 3:
            signal_name = "Signal"
        else:
            signal_name = "+".join([s.displayname for s in signals])
        self.hist_leg.AddEntry(self.signal, signal_name,'l')
        
        if len(backgrounds) > 3:
            bkgd_name = "MC Backgrounds"
        else:
            bkgd_name = "+".join([b.displayname for b in backgrounds])
        self.hist_leg.AddEntry(self.bkgd, bkgd_name, 'l')

        if plot.auto_set_ylimits:
            reformat_axis(plot, self.axis, [self.signal, self.bkgd])
        
        self.signal_eff = self.make_eff_hist(self.signal, xcut_is_max)
        self.bkgd_rej = self.make_eff_hist(self.bkgd, xcut_is_max, bkgd_rej=True)
        self.eff_leg = make_basic_legend(plot)
        self.eff_leg.AddEntry(self.signal_eff, "Signal Eff", 'l')
        self.eff_leg.AddEntry(self.bkgd_rej, "Background Rej", 'l')
        
        self.s_over_b = self.make_s_over_b_hist(self.signal, self.bkgd, xcut_is_max)
        self.zn_sig = self.make_zn_sig_hist(self.signal, self.bkgd, xcut_is_max)
        self.roc_graph = self.make_roc_curve()
    
    def list_of_root_objects(self):
        return [
            self.axis,
            self.hist_leg,
            self.signal,
            self.bkgd,
            self.eff_leg,
            self.signal_eff,
            self.bkgd_rej,
            self.s_over_b,
            self.roc_graph,
        ]

    def make_eff_hist(self, hist, xcut_is_max=True, bkgd_rej=False):
        eff_hist = GetCumulative1D(hist, xcut_is_max)
        total_bin = hist.GetNbinsX() if xcut_is_max else 1
        norm_factor = 1.0 / eff_hist.GetBinContent(total_bin)
        eff_hist.Scale(norm_factor)
        eff_hist.SetMaximum(1.4)
        eff_hist.SetMinimum(0)
        eff_hist.GetYaxis().SetTitle("%")
        
        format_y_axis(eff_hist)

        if bkgd_rej:
            for ibin in range(eff_hist.GetNbinsX()+1):
                bin_val = eff_hist.GetBinContent(ibin)
                eff_hist.SetBinContent(ibin, 1-bin_val)
        else:
            # Print cut values for specific efficiences
            nbins = eff_hist.GetNbinsX()
            for eff in [0.90, 0.95, 0.99]:
                if eff < 0 or 1 < eff: 
                    print "ERROR :: %f efficiency makes no sense. Try again" % (eff)
                if xcut_is_max:
                    ibin_up = eff_hist.FindFirstBinAbove(eff)
                    ibin_dn = ibin_up - 1
                    stating_eff = eff_hist.GetBinContent(nbins - 1)
                    ending_eff = eff_hist.GetBinContent(1)
                else:
                    ibin_dn = eff_hist.FindLastBinAbove(eff)
                    ibin_up = ibin_dn + 1
                    ending_eff = eff_hist.GetBinContent(nbins - 1)
                    stating_eff = eff_hist.GetBinContent(1)

                if ibin_up <= 0:
                    print "INFO :: Signal efficiency begins at %f. Expand hist range to find %f efficiency cut" % (starting_eff, eff)
                    continue
                elif ibin_up == 0 or ibin_dn == nbins-1:
                    print "INFO :: Signal efficiency ends at %f. Decrease bin width to find %f efficiency cut" % (ending_eff, eff)
                    continue
                cut_dn = eff_hist.GetBinLowEdge(ibin_dn)
                cut_up = eff_hist.GetBinLowEdge(ibin_up) + eff_hist.GetBinWidth(ibin_up)
                
                y_up = eff_hist.GetBinContent(ibin_up)
                y_dn = eff_hist.GetBinContent(ibin_dn)
                x_up = eff_hist.GetBinCenter(ibin_up)
                x_dn = eff_hist.GetBinCenter(ibin_dn)
                m = (y_up - y_dn) / (x_up - x_dn)
                cut_guess = (eff - y_dn) / m + x_dn 
                print "INFO :: %f efficiency achieved with cut value between [%f,%f]. Most likely %f" % (eff, cut_dn, cut_up, cut_guess)
            
            # Print cut values for specific efficiences
            for cut in [3.0, 5.0]:
                eff = eff_hist.GetBinContent(eff_hist.FindBin(cut))
                print "INFO :: Cut of %f leads to %f efficiency" % (cut, eff)

        return eff_hist

    def make_s_over_b_hist(self, signal, bkgd, xcut_is_max=True):
        signal_cm = GetCumulative1D(signal, xcut_is_max) 
        bkgd_cm = GetCumulative1D(bkgd, xcut_is_max) 

        s_over_b = signal_cm.Clone(signal.GetName() + "_s_over_b")
        #s_over_b.Divide(bkgd_cm)
        #s_over_b.GetYaxis().SetTitle("S/B")
        bkgd_cm.Add(signal_cm)
        s_over_b.Divide(bkgd_cm) # Purity
        s_over_b.GetYaxis().SetTitle("Purity")
        maxy = min(s_over_b.GetMaximum()*1.2,1)
        s_over_b.SetMaximum(maxy)
        s_over_b.SetMinimum(0)
       
        format_y_axis(s_over_b)
        #format_x_axis(s_over_b)
        xax = s_over_b.GetXaxis()
        xax.SetTitleSize(0.07)
        xax.SetLabelSize(0.08 * 0.81)
        xax.SetLabelOffset(1.15*0.02)
        #xax.SetTitleOffset(0.85 * xax.GetTitleOffset())
        
        signal_cm.Delete()
        bkgd_cm.Delete()

        return s_over_b

    def make_zn_sig_hist(self, signal, bkgd, xcut_is_max):
        signal_cm = GetCumulative1D(signal, xcut_is_max) 
        bkgd_cm = GetCumulative1D(bkgd, xcut_is_max) 

        zn_sig = signal_cm.Clone(signal.GetName() + "_zn_sig")
        if signal.Integral() > bkgd.Integral():
            last_bin = zn_sig.GetNbinsX()
            error = r.Double(0.0)
            bkgd.IntegralAndError(0,-1, error)
            num = 0.1*(signal.Integral()/bkgd.Integral() - 1)
            den = error/bkgd.Integral()
            dummy_scale = num/den 
            zn_sig.GetYaxis().SetTitle("Zn Sig. (scaled)")
        else:
            dummy_scale = 1
            zn_sig.GetYaxis().SetTitle("Zn Sig.")
        for xbin in range(1, zn_sig.GetNbinsX()+1):
            sig_yld = signal_cm.GetBinContent(xbin)
            bkd_yld = bkgd_cm.GetBinContent(xbin)
            bkd_unc = bkgd_cm.GetBinError(xbin)/bkd_yld if bkd_yld else 1
            bkd_unc *= dummy_scale
            zn = ncu.BinomialExpZ(sig_yld, bkd_yld, bkd_unc)
            if zn == float('inf'):
                zn = 0
            zn_sig.SetBinContent(xbin, zn)
            #print "TESTING :: sig_yld = %f, bkd_yld = %f, bkd_unc = %f, Zn = %f" % (sig_yld, bkd_yld, bkd_unc, zn)
        maxy = zn_sig.GetMaximum()
        maxy = maxy + 0.2*abs(maxy)
        miny = zn_sig.GetMaximum()
        miny = miny - 0.2*abs(miny)
        zn_sig.SetMaximum(maxy)
        zn_sig.SetMinimum(0)
       
        format_y_axis(zn_sig)
        xax = zn_sig.GetXaxis()
        xax.SetTitleSize(0.07)
        xax.SetLabelSize(0.08 * 0.81)
        xax.SetLabelOffset(1.15*0.02)
        #xax.SetTitleOffset(0.85 * xax.GetTitleOffset())
        
        signal_cm.Delete()
        bkgd_cm.Delete()

        return zn_sig

    def make_roc_curve(self):
        roc_plot = r.TGraph()
        skipped_bins = 0
        for xbin in range(1, self.signal.GetNbinsX()+1):
            for ybin in range(1, self.signal.GetNbinsY()+1):
                x_val = self.signal_eff.GetBinContent(xbin,ybin) 
                b_eff = 1 - self.bkgd_rej.GetBinContent(xbin,ybin) 
                # Ignore large spike in background rejection that occurs
                # after the cut has already removed 95% of signal
                # Also ignore cases where bkgd efficiency is negative because of negative weights
                if b_eff > 0 and x_val > 0.05:
                    y_val = 1.0 / (b_eff)
                    roc_plot.SetPoint(xbin-1-skipped_bins, x_val, y_val)
                else:
                    skipped_bins += 1
        roc_plot.GetXaxis().SetTitle("Signal Efficiency") 
        roc_plot.GetYaxis().SetTitle("Background Rejection") 
        format_y_axis(roc_plot)
        roc_plot.GetYaxis().SetLabelSize(0.06 * 0.81)
        roc_plot.GetYaxis().SetTitleSize(0.06 * 0.81)
        roc_plot.GetYaxis().SetTitleOffset(0.80 * 1.2)
        format_x_axis(roc_plot)
        roc_plot.GetXaxis().SetLabelSize(0.06 * 0.81)
        roc_plot.GetXaxis().SetTitleSize(0.06 * 0.81)
        
        return roc_plot
   
class CutScan2D(HistBase):
    def __init__(self, plot, reg, signals, backgrounds, 
            xcut_is_max=True, ycut_is_max=True, and_cuts=True, bkgd_rej=False):
        '''
        Params:
            - xcut_is_max (bool) - Indicates cuts are of the form cut_var < cut_value
        '''
        self.axis = None
        self.signal = None
        self.bkgd = None
        self.signal_eff = None
        self.bkgd_rej = None
        self.s_over_b = None
        self.roc_graph = None

        self.axis = make_plot2D_axis(plot)
        
        self.signal = make_hist2D(plot, reg, signals)
        self.bkgd = make_hist2D(plot, reg, backgrounds)
        
        self.signal_eff = self.make_eff_hist(self.signal, xcut_is_max, ycut_is_max, and_cuts)
        self.bkgd_rej = self.make_eff_hist(self.bkgd, xcut_is_max, ycut_is_max, and_cuts, bkgd_rej=True)
        
        self.s_over_b = self.make_s_over_b_hist(self.signal, self.bkgd, xcut_is_max, ycut_is_max, and_cuts)
        self.roc_graph = self.make_roc_curve()
        
        if len(signals) > 3:
            signal_name = "Signal"
        else:
            signal_name = "+".join([s.displayname for s in signals])
        self.signal.SetTitle(signal_name)
        
        if len(backgrounds) > 3:
            bkgd_name = "MC Backgrounds"
        else:
            bkgd_name = "+".join([b.displayname for b in backgrounds])
        self.bkgd.SetTitle(bkgd_name)

    def list_of_root_objects(self):
        return [
            self.axis,
            self.signal,
            self.bkgd,
            self.signal_eff,
            self.bkgd_rej,
            self.s_over_b,
            self.roc_graph,
        ]

    def make_eff_hist(self, hist, xcut_is_max=True, ycut_is_max=True, and_cuts=True, bkgd_rej=False):
        '''
        Params:
            xcut_is_max - xvar cuts are treated as a min (i.e. cut_var < cut_val)
        '''
        eff_hist = GetCumulative2D(hist, xcut_is_max, ycut_is_max, and_cuts)
        total_binx = hist.GetNbinsX()+1 if xcut_is_max else 0
        total_biny = hist.GetNbinsY()+1 if ycut_is_max else 0
        norm_factor = 1.0 / eff_hist.GetBinContent(total_binx, total_biny)
        eff_hist.Scale(norm_factor)
        eff_hist.SetMaximum(1.0)
        eff_hist.SetMinimum(0)
        eff_hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
        eff_hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
        eff_hist.GetZaxis().SetTitle("%")

        
        format_y_axis(eff_hist)

        if bkgd_rej:
            for xbin in range(eff_hist.GetNbinsX()+2):
                for ybin in range(eff_hist.GetNbinsY()+2):
                    bin_val = eff_hist.GetBinContent(xbin, ybin)
                    eff_hist.SetBinContent(xbin, ybin, 1-bin_val)

        return eff_hist

    def make_s_over_b_hist(self, signal, bkgd, xcut_is_max, ycut_is_max, and_cuts):
        signal_cm = GetCumulative2D(signal, xcut_is_max, ycut_is_max, and_cuts) 
        bkgd_cm = GetCumulative2D(bkgd, xcut_is_max, ycut_is_max, and_cuts) 

        s_over_b = signal_cm.Clone(signal.GetName() + "_s_over_b")
        #s_over_b.Divide(bkgd_cm)
        #s_over_b.GetYaxis().SetTitle("S/B")
        bkgd_cm.Add(signal_cm)
        s_over_b.Divide(bkgd_cm) # Purity
        s_over_b.GetZaxis().SetTitle("Purity")
        maxy = min(s_over_b.GetMaximum(),1)
        s_over_b.SetMaximum(maxy)
        s_over_b.SetMinimum(0)
       
        format_y_axis(s_over_b)
        format_x_axis(s_over_b)
        
        signal_cm.Delete()
        bkgd_cm.Delete()

        return s_over_b

    def make_roc_curve(self):
        roc_plot = r.TGraph()
        bin_counter = 0
        for xbin in range(1, self.signal.GetNbinsX()+1):
            for ybin in range(1, self.signal.GetNbinsY()+1):
                x_val = self.signal_eff.GetBinContent(xbin, ybin) 
                b_eff = 1 - self.bkgd_rej.GetBinContent(xbin, ybin) 
                # Ignore large spike in background rejection that occurs
                # after the cut has already removed 95% of signal
                # Also ignore cases where bkgd efficiency is negative because of negative weights
                if b_eff > 0 and x_val > 0.05: #Ignore large spike in background rejection
                    y_val = 1.0 / (b_eff)
                    roc_plot.SetPoint(bin_counter, x_val, y_val)
                    bin_counter += 1
        roc_plot.GetXaxis().SetTitle("Signal Efficiency") 
        roc_plot.GetYaxis().SetTitle("Background Rejection") 
        format_y_axis(roc_plot)
        roc_plot.GetYaxis().SetLabelSize(0.06 * 0.81)
        roc_plot.GetYaxis().SetTitleSize(0.06 * 0.81)
        roc_plot.GetYaxis().SetTitleOffset(0.80 * 1.2)
        format_x_axis(roc_plot)
        roc_plot.GetXaxis().SetLabelSize(0.06 * 0.81)
        roc_plot.GetXaxis().SetTitleSize(0.06 * 0.81)
        
        return roc_plot

def GetCumulative1D(hist, cut_is_max):
    '''
    Get cumulative histogram. Differs from TH1::GetCumulative by including overflow
    Params:
        cut_is_max (bool) - cuts on x-axis are treated as minimum cuts (x_var < x_cut_val)
    '''
    cumulative_hist = hist.Clone(hist.GetName() + "_cumulative")
    cumulative_hist.Reset()

    of_xbin = hist.GetNbinsX()+1 # overflow x-bin
    error = r.Double(0.0)
    for xbin in range(of_xbin+1):
        x1, x2 = (0, xbin) if cut_is_max else (xbin, of_xbin)
        val = hist.IntegralAndError(x1, x2, error)
        cumulative_hist.SetBinContent(xbin, val)
        cumulative_hist.SetBinError(xbin, error)

    return cumulative_hist

def GetCumulative2D(hist, xcut_is_max, ycut_is_max, and_cuts):
    '''
    Params:
        xcut_is_max (bool) - cuts on x-axis are treated as minimum cuts (x_var < x_cut_val)
        ycut_is_max (bool) - cuts on y-axis are treated as minimum cuts (y_var < y_cut_val)
        and_cuts (bool) - x and y cuts are AND'ed vs OR'ed
    '''
    cumulative_hist = hist.Clone(hist.GetName() + "_cumulative")
    cumulative_hist.Reset()

    of_xbin = hist.GetNbinsX()+1 # overflow x-bin
    of_ybin = hist.GetNbinsY()+1 # overfloy y-bin
    total = hist.Integral(0, -1, 0, -1)
    for xbin in range(of_xbin+1):
        for ybin in range(of_ybin+1):
            x1, x2 = (0, xbin) if xcut_is_max else (xbin, of_xbin)
            y1, y2 = (0, ybin) if ycut_is_max else (ybin, of_ybin)
            val = hist.Integral(x1,x2,y1,y2)

            if not and_cuts:
                # Cuts are OR'ed so add regions that only fail one criterion
                x1_flip, x2_flip = (xbin+1, -1) if xcut_is_max else (0, xbin-1) 
                y1_flip, y2_flip = (ybin+1, -1) if ycut_is_max else (0, ybin-1) 
                val_x = hist.Integral(x1_flip,x2_flip,y1,y2)
                val_y = hist.Integral(x1,x2,y1_flip,y2_flip)

                # If those regions are outide of the hist overflow, then skip
                if x2-x1 >= of_xbin: val_x = 0
                if y2-y1 >= of_ybin: val_y = 0
                
                val += val_x + val_y

            print "xbin = %d, ybin = %d, cumulative = %.2f" % (xbin, ybin, val)
            cumulative_hist.SetBinContent(xbin, ybin, val)

    return cumulative_hist

def format_y_axis(hist):
    yax = hist.GetYaxis()
    yax.SetTitleSize(0.10 * 0.83)
    yax.SetLabelSize(0.08 * 0.81)
    yax.SetLabelOffset(0.98 * 0.013 * 1.08)
    yax.SetTitleOffset(0.45 * 1.2)

def format_x_axis(hist):
    xax = hist.GetXaxis()
    xax.SetTitleSize(1.1 * 0.14)
    xax.SetLabelSize(0.08 * 0.81)
    xax.SetLabelOffset(1.15*0.02)
    #xax.SetTitleOffset(0.85 * xax.GetTitleOffset())

class RatioHist1D(HistBase) :
    ratio_ymax = 2.0
    ratio_ymin = 0.5

    def __init__(self, plot, num, den, ymax = 0, ymin = 0):
        self.axis = None
        self.ratio = None

        ymax = ymax if ymax else self.ratio_ymax
        ymin = ymin if ymin else self.ratio_ymin
        
        self.axis = make_ratio1D_axis(plot, ymin, ymax) 
        self.ratio = make_ratio_graph(num, den)

    def list_of_root_objects(self):
        return [
            self.axis,
            self.ratio
        ]

def make_ratio_graph(num, den):
    ratio = num.Clone("ratio_graph")
    ratio.Divide(den)
    ratio_graph = pu.th1_to_tgraph(ratio)
    return ratio_graph

def make_ratio1D_axis(plot, ymin = 0, ymax = 0):
    hax = r.TH1F("ratio_axis", "", plot.nbins, plot.bin_edges)
    hax.SetMinimum(ymin)
    hax.SetMaximum(ymax)

    xax = hax.GetXaxis()
    xax.SetTitle(plot.xlabel)

    yax = hax.GetYaxis()
    yax.SetTitle(plot.ylabel)

    return hax

class RegionCompare1D(HistBase):
    def __init__(self, regions, samples, plot, event_list_dir = "./"):
        self.axis = None
        self.leg = None
        self.hists = []
        
        self.axis = make_plot1D_axis(plot)
        self.leg = make_basic_legend(plot)
        if not regions or not samples: 
            print "WARNING :: Either regions or samples is empty for plot:", plot.var
            return
        for ii, reg in enumerate(regions):
            comb_hist = None
            for sample in samples:
                print "Setting EventLists for %s"%sample.name
                list_name = "list_" + reg.name + "_" + sample.name
                sample.set_event_list(reg.tcut, list_name, event_list_dir)
                hist = make_basic_hist(plot, sample, reg)
                if comb_hist:
                    comb_hist.Add(hist)
                    hist.Delete()
                else:
                    comb_hist = hist
            comb_hist.SetLineStyle(ii%10 + 1)
            self.leg.AddEntry(comb_hist, reg.displayname, "l")
            self.hists.append(comb_hist)
        
        #if plot.doNorm:
        if True:
            for hist in self.hists:
                normalize_hist(hist)
        
        if plot.auto_set_ylimits:
            reformat_axis(plot, self.axis, self.hists)

    def list_of_root_objects(self):
        return [
            self.leg,
            self.axis,
            self.hists
        ]

def reformat_axis(plot, axis, hists):
    ''' Reformat axis to fit content and labels'''
    # Get max y-value in hists
    maxy = 0
    miny = float("inf")
    for hist in hists:
        maxy_tmp = hist.GetMaximum()
        miny_tmp = hist.GetMinimum()    
        if maxy_tmp > maxy: maxy = maxy_tmp
        if miny_tmp < miny: miny = miny_tmp
    
    if maxy <= 0:
        print "WARNING :: Max value of plot is <= 0"
    
    plot_maxy, plot_miny = hist_yrange_to_plot_yrange(maxy, miny, plot.doLogY)

    # reformat the axis
    axis.SetMaximum(plot_maxy)
    axis.SetMinimum(plot_miny)

def hist_yrange_to_plot_yrange(hist_max, hist_min, log_plot=False, big_legend=False):
    if log_plot:
        plot_max = 10**(pu.get_order_of_mag(hist_max))
        plot_min = 1e-2
        #if hist_min > 0:
        #    plot_min = 10**(pu.get_order_of_mag(hist_min)-1)
        #else:
        #    plot_min = 10**(pu.get_order_of_mag(hist_max) - 4)
    else:
        plot_max = hist_max
        plot_min = 0

    # Get y-axis max multiplier to fit labels
    if log_plot:
        max_mult = 1e6 if big_legend else 1e5
    else:
        max_mult = 2.0 if big_legend else 1.8
    
    return plot_max * max_mult, plot_min

def normalize_hist(hist):
    norm_factor = 1.0/hist.Integral() if hist.Integral() else 1
    if isinstance(hist, r.TH1):
        hist.Scale(norm_factor)
    elif isinstance(hist, r.TGraph):
        pu.scale_tgraph(hist, norm_factor)
    elif isinstance(hist, r.THStack):
        pu.scale_thstack(hist, norm_factor)
    else:
        print "WARNING :: Unexpectd class type:", type(hist)

def make_hist1D(plot, reg, samples, samples_name=""):
    # Make simple hist
    if not isinstance(samples, list):
        samples = [samples]
    if not samples_name:
        samples_name = "_".join([x.name for x in samples])

    var_name = pu.strip_for_root_name(plot.variable)
    h_name = "h_"+reg.name+'_'+samples_name+"_"+ var_name
    #TODO Is setting the label here redundant since it probably is set later?
    h_label = ";%s;%s" % (plot.xlabel, plot.ylabel)
    h = r.TH1D(h_name,h_label, plot.nbins, plot.bin_edges)
    h.SetLineColor(samples[0].color)
    h.Sumw2

    draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
    
    for sample in samples:
        # Draw final histogram (i.e. selections and weights applied)
        if samples[0].weight_str:
            weight_str = "%s * %s"%(sample.weight_str, str(sample.scale_factor))
        else:
            weight_str = "1"

        sample.tree.Draw(draw_cmd, weight_str, "goff")

    # Rebin
    if plot.rebin_bins:
        new_bins = array('d', plot.rebin_bins)
        h = h.Rebin(len(new_bins)-1, h_name, new_bins)

    # Add overflow
    if plot.add_overflow:
        pu.add_overflow_to_lastbin(h)
    if plot.add_underflow:
        pu.add_underflow_to_firstbin(h)

    return h

class DataMCRatioHist1D(RatioHist1D) :
    def __init__(self, plot, reg, stack_hist):
        #TODO: take samples as input as opposed to stack hist
        self.axis = None
        self.errors = None
        self.ratio = None
        if not stack_hist.mc_stack: return

        self.ratio_label = 'Data/MC'
        self.get_ratio_axis(plot, stack_hist.mc_stack)
        self.get_ratio_errors(stack_hist.mc_errors)
        self.get_ratio_graph(stack_hist.data_hist, stack_hist.mc_errors)

    def list_of_root_objects(self):
        return [
            self.axis,
            self.errors,
            self.ratio
        ]

    def get_ratio_axis(self, plot, stack):
        # yaxis
        #if not stack:
        #    return
        self.axis = stack.GetStack().Last().Clone("h_sm")
        h_sm = self.axis #TODO: remove relabel
        yax = h_sm.GetYaxis()
        yax.SetRangeUser(self.ratio_ymin, self.ratio_ymax)
        yax.SetTitle(self.ratio_label)
        yax.SetTitleSize(0.11)
        yax.SetLabelSize(0.11)
        #yax.SetLabelOffset(0.98 * 0.013 * 1.08)
        yax.SetTitleOffset(0.5)
        yax.SetNdivisions(5)

        # xaxis
        xax = h_sm.GetXaxis()
        xax.SetTitleSize(0.11)
        xax.SetLabelSize(0.11)
        xax.SetLabelOffset(0.025)
        xax.SetTitleOffset(1.3)

        #h_sm.SetTickLength(0.06)

        if plot.bin_labels:
            plot.set_bin_labels(h_sm)

        return h_sm

    def get_ratio_errors(self, mc_errors):
        self.errors = r.TGraphAsymmErrors(mc_errors)
        pu.buildRatioErrorBand(mc_errors, self.errors)

    def get_ratio_graph(self, data_hist, mc_errors):
        g_data = pu.convert_errors_to_poisson(data_hist)
        #g_sm = pu.th1_to_tgraph(h_sm)
        #g_ratio = pu.tgraphAsymmErrors_divide(g_data, g_sm)

        # For Data/MC only use the statistical error for data
        # since we explicity draw the MC error band
        nominalAsymErrorsNoSys = r.TGraphAsymmErrors(mc_errors)
        for i in xrange(nominalAsymErrorsNoSys.GetN()) :
            nominalAsymErrorsNoSys.SetPointError(i-1,0,0,0,0) # TODO change i-1 to i
        ratio_raw = pu.tgraphAsymmErrors_divide(g_data, nominalAsymErrorsNoSys)
        #ratio_raw = pu.tgraphAsymmErrors_divide(nominalAsymErrorsNoSys,g_data)

        # Make final ratio plot
        self.ratio = r.TGraphAsymmErrors()
        x1, y1 = r.Double(0.0), r.Double(0.0)
        index = 0
        for i in xrange(ratio_raw.GetN()) :
            ratio_raw.GetPoint(i, x1, y1)
            if y1 <= 0. : y1 = r.Double(-2.0)

            self.ratio.SetPoint(index, x1, y1)
            xlo, xhi = ratio_raw.GetErrorXlow(i), ratio_raw.GetErrorXhigh(i)
            ylo, yhi = ratio_raw.GetErrorYlow(i), ratio_raw.GetErrorYhigh(i)
            self.ratio.SetPointError(index, xlo, xhi, ylo, yhi)
            index+=1

        # Format
        self.ratio.SetLineWidth(1)
        self.ratio.SetMarkerStyle(20)
        self.ratio.SetMarkerSize(1.5)
        self.ratio.SetLineColor(1)
        self.ratio.SetMarkerSize(1.5)

        # Clean up
        ratio_raw.Delete()
        nominalAsymErrorsNoSys.Delete()
        g_data.Delete()

def make_plot1D_axis(plot):
    hax = r.TH1F("axes", "", plot.nbins, plot.bin_edges)
    #hax = r.TH1F("axes", "", int(plot.nbins), plot.xmin, plot.xmax)
    hax.SetMinimum(plot.ymin)
    hax.SetMaximum(plot.ymax)
    xax = hax.GetXaxis()
    xax.SetTitle(plot.xlabel)
    xax.SetTitleSize(0.05)
    xax.SetTitleOffset(1)
    if plot.ptype == Types.ratio:
        xax.SetTitleOffset(-999)
        xax.SetLabelOffset(-999)

    yax = hax.GetYaxis()
    yax.SetTitle(plot.ylabel)
    yax.SetTitleOffset(1.7)
    yax.SetTitleSize(0.035)

    if plot.bin_labels:
        plot.set_bin_labels(hax)
        
    if plot.rebin_bins:
        new_bins = array('d', plot.rebin_bins)
        hax_rebinned = hax.Rebin(len(new_bins)-1, 'axes_rebinned', new_bins)
        hax.Delete()
        hax = hax_rebinned
    return hax

class SampleCompare1D(HistBase):
    def __init__(self, plot, reg, samples):
        self.axis = None
        self.leg = None
        self.hists = []

        self.axis = make_plot1D_axis(plot)
        self.leg = make_basic_legend(plot)
        for sample in samples:
            hist = make_basic_hist(plot, sample, reg)
            self.leg.AddEntry(hist, sample.displayname, "l")
            
            # Add overflow/underflow
            if plot.add_overflow:
                pu.add_overflow_to_lastbin(hist)
            if plot.add_underflow:
                pu.add_underflow_to_firstbin(hist)
            
            # Normalize
            if plot.doNorm: 
                normalize_hist(hist)
        
            self.hists.append(hist)

        if plot.auto_set_ylimits:
            reformat_axis(plot, self.axis, self.hists)

    def list_of_root_objects(self):
        return [
                self.axis,
                self.leg,
                self.hists
                ]

class DataMCStackHist1D(HistBase):
    def __init__(self, plot, reg, data, bkgds, sigs=None):
        # Get plotting primitives
        # Not all primitives are for drawing. Some are for preserving pointers
        self.leg_sig = None
        self.leg = None
        self.axis = None
        self.mc_stack = None
        self.mc_total = None
        self.data = None
        self.signals = []
        self.data_hist = None
        self.error_leg = None
        self.mc_errors = None
        self.histos_for_leg = []

        self.make_stack_legend(plot)
        self.axis = make_plot1D_axis(plot)
        self.add_stack_backgrounds(plot, reg, bkgds)
        if sigs: self.add_stack_signals(plot, reg, sigs)
        self.add_stack_data(plot, reg, data)
        self.add_stack_mc_errors(plot, reg)
        if plot.doNorm:
            self.normalize_hists()
        if plot.auto_set_ylimits:
            self.reformat_axis(plot)

    def list_of_root_objects(self):
        return [
            self.leg_sig,
            self.leg,
            self.axis,
            self.mc_stack,
            self.mc_total,
            self.data,
            self.data_hist,
            self.error_leg,
            self.mc_errors
        ]

    def make_stack_legend(self, plot):
        if plot.leg_is_left :
            self.leg = pu.default_legend(xl=0.2,yl=0.7,xh=0.47, yh=0.87)
        elif plot.leg_is_bottom_right :
            self.leg = pu.default_legend(xl=0.7, yl=0.17,xh=0.97,yh=0.41)
        elif plot.leg_is_bottom_left :
            self.leg = pu.default_legend(xl=0.2,yl=0.2,xh=0.47,yh=0.37)
        else :
            self.leg = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)

        self.leg.SetNColumns(2)
        # TODO: Incorporate signal legend
        self.leg_sig = pu.default_legend(xl=0.55, yl=0.60, xh=0.91, yh=0.71)
        self.leg_sig.SetNColumns(1)

    def add_stack_backgrounds(self, plot, reg, bkgds):

        # Initilize lists and defaults
        histos = []

        # Make sample hists
        for sample in bkgds :
            # Initilize histogram
            h = make_hist1D(plot, reg, sample)

            h.GetXaxis().SetLabelOffset(-999)

            if sample.isMC and sample.isSignal:
                h.SetLineWidth(2)
                h.SetLineStyle(2)
                h.SetFillStyle(0)
            else:
                h.SetFillColor(sample.color)
                h.SetFillStyle(1001)

            h.leg_name = sample.displayname #dynamic class members...ooo yeah!

            # Record all and non-empty histograms
            if sample.isMC and sample.isSignal:
                leg_sig.AddEntry(h, sample.displayname, "l")
            elif h.Integral() > 0:
                histos.append(h) 
        if not len(histos):
            print "ERROR (make_stack_background) :: All SM hists are empty. Skipping"
            return

        # Order the hists by total events
        histos = sorted(histos, key=lambda h: h.Integral())

        self.mc_stack = r.THStack("stack_"+plot.name, "")
        for h in histos :
            self.mc_stack.Add(h)

        h_leg = sorted(histos, key=lambda h: h.Integral(), reverse=True)
        self.histos_for_leg = self.arrange_histos_for_legend(h_leg)

        # draw the total bkg line
        self.mc_total = self.mc_stack.GetStack().Last().Clone("hist_sm")
        self.mc_total.SetLineColor(r.kBlack)
        self.mc_total.SetLineWidth(3)
        self.mc_total.SetLineStyle(1)
        self.mc_total.SetFillStyle(0)
        self.mc_total.SetLineWidth(3)

    def add_stack_signals(self, plot, reg, signals):
        # Make MC sample hists
        for sig_sample in signals:
            # Initilize histogram
            h = make_hist1D(plot, reg, sig_sample)
            
            h.GetXaxis().SetLabelOffset(-999)
            h.SetLineWidth(2)
            h.SetLineStyle(2)
            h.SetFillStyle(0)
            
            h.leg_name = sig_sample.displayname #dynamic class members...ooo yeah!

            # Record all and non-empty histograms
            self.leg_sig.AddEntry(h, sig_sample.displayname, "l")
            self.signals.append(h)
    
    def add_stack_data(self, plot, reg, data):
        if not data: return
        #TODO: Look for a way to combine with backgrounds
        if data.blinded and reg.isSR:
            data.scale_factor = 0
        self.data_hist = make_hist1D(plot, reg, data)
        self.data_hist.GetXaxis().SetLabelOffset(-999)

        self.data = pu.convert_errors_to_poisson(self.data_hist)
        #self.data.SetLineWidth(2)
        
        #uglify
        self.data.SetLineWidth(1)
        self.data.SetMarkerStyle(20)
        self.data.SetMarkerSize(1.5)
        self.data.SetLineColor(1)
        self.leg.AddEntry(self.data, "Data", "p")

    def add_stack_mc_errors(self, plot, reg):
        if not self.mc_stack:
            return
        r.gStyle.SetHatchesSpacing(0.9)

        self.error_leg = r.TH1F("mcError", "mcError", 2,0,2)
        self.error_leg.SetLineWidth(3)
        self.error_leg.SetFillStyle(3345)
        self.error_leg.SetFillColor(r.kBlue)
        self.error_leg.SetLineColor(r.kBlack)
        self.leg.AddEntry(self.error_leg, "Standard Model", "fl")

        # now add backgrounds to legend
        for h in self.histos_for_leg :
            if h.Integral(0, -1) <= 0: continue
            self.leg.AddEntry(h, h.leg_name, "f")

        try:
            totalSM = self.mc_stack.GetStack().Last().Clone("totalSM")
        except ReferenceError:
            print "WARNING (add_stack_mc_errors) :: Unable to access stack"
            return
        self.mc_errors = pu.th1_to_tgraph(totalSM)
        totalSM.Delete()
        self.mc_errors.SetMarkerSize(0)
        self.mc_errors.SetLineWidth(0)
        self.mc_errors.SetFillStyle(3345)
        self.mc_errors.SetFillColor(r.kBlue)

            # symmetrize the errors
        for i in xrange(self.mc_errors.GetN()) :
            ehigh = self.mc_errors.GetErrorYhigh(i)
            elow  = self.mc_errors.GetErrorYlow(i)


            error_sym = r.Double(0.0)
            error_sym = (ehigh + elow) / 2.0

            if ehigh != error_sym:
                print "initial error (+%.2f,-%.2f), symmetrized = (+%.2f,-%.2f)"%(ehigh,elow, error_sym, error_sym)


            self.mc_errors.SetPointEYhigh(i,0.0)
            self.mc_errors.SetPointEYhigh(i, error_sym)
            self.mc_errors.SetPointEYlow(i,0.0)
            self.mc_errors.SetPointEYlow(i,error_sym)
        return

    def normalize_hists(self):
        if self.mc_total and self.mc_total.Integral():
            mc_norm_factor = 1.0/self.mc_total.Integral()
            pu.scale_thstack(self.mc_stack, mc_norm_factor)
            self.mc_total.Scale(mc_norm_factor)
            pu.scale_tgraph(self.mc_errors, mc_norm_factor)
        for s in self.signals:
            sig_norm_factor = 1.0/s.Integral() if s.Integral() else 1
            s.Scale(sig_norm_factor)
        if self.data_hist and self.data_hist.Integral():
            data_norm_factor = 1.0/self.data_hist.Integral()
            pu.scale_tgraph(self.data, data_norm_factor)

    def reformat_axis(self, plot):
        ''' Reformat axis to fit content and labels'''
        if not self.mc_stack:
            return
        # Get maximum histogram y-value
        data_maxy = max(pu.get_tgraph_max(self.data), self.mc_stack.GetMaximum())
        data_miny = min(pu.get_tgraph_min(self.data), self.mc_stack.GetMinimum())
        mc_maxy = self.mc_stack.GetMaximum()
        mc_miny = self.mc_stack.GetMinimum("nostack")

        maxy = max(data_maxy, mc_maxy)
        miny = min(data_miny, mc_miny)

        if maxy <= 0:
            print "WARNING :: Max value of plot is <= 0"
            return

        if plot.doNorm:
            plot_maxy = plot.logy_norm_max
            plot_miny = plot.logy_norm_min
        else:
            plot_maxy, plot_miny = hist_yrange_to_plot_yrange(maxy, miny, plot.doLogY, bool(self.signals))

        # reformat the axis
        self.mc_stack.SetMaximum(plot_maxy)
        self.mc_stack.SetMinimum(plot_miny)
        self.axis.SetMaximum(plot_maxy)
        self.axis.SetMinimum(plot_miny)

    def arrange_histos_for_legend(self, histos) :
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

class ComparisonHist1D :
    def __init__(self):
        pass

################################################################################
# 2D Plots
################################################################################
def make_plot2D_axis(plot):
    axis = r.TH2D("axes", "", plot.nxbins, plot.xbin_edges, plot.nybins, plot.ybin_edges)
    #axis.SetMinimum(plot.zmin)
    #axis.SetMaximum(plot.zmax)

    xax = axis.GetXaxis()
    xax.SetTitle(plot.xlabel)
    xax.SetLabelSize(0.035)
    xax.SetTitleSize(0.041)
    xax.SetLabelOffset(0.023)
    xax.SetTitleOffset(2.1)

    #if plot.bin_labels:
    #    plot.set_bin_labels(axis)

    yax = axis.GetYaxis()
    yax.SetTitle(plot.ylabel)
    yax.SetTitleOffset(1.4)
    yax.SetLabelOffset(0.013)
    yax.SetLabelSize(0.042)
    yax.SetTitleSize(0.47)

    zax = axis.GetZaxis()
    zax.SetTitle(plot.zlabel)
    #zax.SetTitleOffset(1.4)
    #zax.SetLabelOffset(0.013)
    #zax.SetLabelSize(1.2 * 0.035)
    #zax.SetTitleSize(0.055 * 0.85)

    #if plot.bin_labels:
    #    plot.set_ybin_labels(axis)

    #if plot.rebin_xbins:
    #    new_bins = array('d', plot.rebin_xbins)
    #    axis = axis.RebinX(len(new_bins)-1, 'axes', new_bins)
    #if plot.rebin_ybins:
    #    new_bins = array('d', plot.rebin_ybins)
    #    axis = axis.RebinY(len(new_bins)-1, 'axes', new_bins)

    return axis

def make_hist2D(plot, reg, samples):
    #TODO: require region_name input instead of reg object
    if not isinstance(samples, list):
        samples = [samples]
    hist = r.TH2D(plot.name, "", plot.nxbins, plot.xmin, plot.xmax, plot.nybins, plot.ymin, plot.ymax)
    x_var = pu.strip_for_root_name(plot.xvariable)
    y_var = pu.strip_for_root_name(plot.yvariable)
    for sample in samples:
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

        hist.Add(h_tmp)

    zax = hist.GetZaxis()
    zax.SetTitle(plot.zlabel)
    zax.SetTitleOffset(1.5)
    zax.SetLabelOffset(0.013)
    zax.SetLabelSize(1.2 * 0.035)

    return hist

class Hist2D(HistBase) :
    def __init__(self, plot, reg, samples):
        #TODO: make simplest class take a single sample
        # and add derived class that takes lists of samples to add them
        # or write function that
        self.axis = None
        self.hist = None
        if not samples: return
        self.axis = make_plot2D_axis(plot)
        self.make_hist(plot, reg, samples)

    def list_of_root_objects(self):
        return [self.axis, self.hist]


    def make_hist(self, plot, reg, samples):
        self.hist = r.TH2D(plot.name, "", plot.nxbins, plot.xbin_edges, plot.nybins, plot.ybin_edges)
        for sample in samples:

            x_var = pu.strip_for_root_name(plot.xvariable)
            y_var = pu.strip_for_root_name(plot.yvariable)
            h_name_tmp = "h_"+reg.name+'_'+sample.name+"_"+x_var+"_"+y_var
            h_tmp = r.TH2D(plot.name, "", plot.nxbins, plot.xbin_edges, plot.nybins, plot.ybin_edges)
            # Draw final histogram (i.e. selections and weights applied)


            if not sample.isMC:
                weight_str = '1'
            else:
                weight_str = "%s * %s"%(sample.weight_str, str(sample.scale_factor))
            draw_cmd = "%s>>%s"%(plot.yvariable+":"+plot.xvariable, h_tmp.GetName())
            sample.tree.Draw(draw_cmd, weight_str, "goff")

            # Yield +/- stat error
            stat_err = r.Double(0.0)
            integral = h_tmp.IntegralAndError(0,-1,0,-1,stat_err)
            self.hist.Add(h_tmp)

        zax = self.hist.GetZaxis()
        zax.SetTitle(plot.zlabel)
        zax.SetTitleOffset(1.5)
        zax.SetLabelOffset(0.013)
        zax.SetLabelSize(1.2 * 0.035)

def make_basic_hist(plot, sample, reg, apply_cuts=False):
    var_name = pu.strip_for_root_name(plot.variable)
    h_name = "h_"+reg.name+'_'+sample.name+"_"+ var_name
    h = pu.th1d(h_name, "", int(plot.nbins),
                plot.xmin, plot.xmax,
                plot.xlabel, plot.ylabel)

    h.SetLineColor(sample.color)
    if sample.isMC and sample.isSignal:
        h.SetLineStyle(2)
    h.SetFillColor(0)
    h.Sumw2

    # Draw final histogram (i.e. selections and weights applied)
    if hasattr(sample, 'weight_str') and plot.variable != sample.weight_str:
        weight_str = "%s * %s"%(sample.weight_str, str(sample.scale_factor))
    else:
        weight_str = "1"

    if apply_cuts:
        weight = "(%s) * %s"%(reg.tcut, weight_str)
        draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
    else:
        weight = weight_str
        draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
    sample.tree.Draw(draw_cmd, weight, "goff")
    
    if plot.doNorm: normalize_hist(h)
            
    return h

def make_basic_legend(plot):
    if plot.leg_is_left :
        leg = pu.default_legend(xl=0.2,yl=0.7,xh=0.47, yh=0.87)
    elif plot.leg_is_bottom_right :
        leg = pu.default_legend(xl=0.7, yl=0.17,xh=0.97,yh=0.41)
    elif plot.leg_is_bottom_left :
        leg = pu.default_legend(xl=0.2,yl=0.2,xh=0.47,yh=0.37)
    else :
        leg = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)
    return leg
