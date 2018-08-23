import re
from array import array

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
        self.clear()
    def list_of_root_objects(self):
        " Return list of all class members that are or contain root objects"
        #TODO: Get python to get all class variables that inherit fromTObject
        #      (e.g. use self.__class__.__dict__)
        raise NotImplementedError()
    
    def clear(self):
        # Delete all underlying root objects
        # -
        for ro in self.list_of_root_objects():
            # Variables
            if not ro or isinstance(ro, r.THStack):
                # THStack is just a container for the hists in the stack and
                # therefore does not need to be deleted
                continue
            elif isinstance(ro, list):
                root_delete(ro)
            elif issubclass(type(ro), r.TObject):
                ro.Delete()
            else:
                print "ERROR :: Unknown primitive type", type(ro)

#TODO: Remover unused classes
class StackHist1D(HistBase):
    def __init__(self):
        pass

class RatioHist1D(HistBase) :
    ratio_ymax = 2.0

    def __init__(self, plot, reg, num, den):
        self.ratio_label = 'Num / Den'
        #self.get_ratio_axis(plot)
        #self.get_ratio_errors(plot)
        #self.get_ratio_graph(plot)


    def list_of_root_objects(self):
        return [
            self.axis,
            self.errors,
            self.ratio
        ]

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
        yax.SetRangeUser(0, self.ratio_ymax)
        yax.SetTitle(self.ratio_label)
        yax.SetTitleSize(0.14 * 0.83)
        yax.SetLabelSize(0.13 * 0.81)
        yax.SetLabelOffset(0.98 * 0.013 * 1.08)
        yax.SetTitleOffset(0.45 * 1.2)
        yax.SetLabelFont(42)
        yax.SetTitleFont(42)
        yax.SetNdivisions(5)

        # xaxis
        xax = h_sm.GetXaxis()
        xax.SetTitleSize(1.1 * 0.14)
        xax.SetLabelSize(yax.GetLabelSize())
        xax.SetLabelOffset(1.15*0.02)
        xax.SetTitleOffset(0.85 * xax.GetTitleOffset())
        xax.SetLabelFont(42)
        xax.SetTitleFont(42)

        h_sm.SetTickLength(0.06)

        if plot.bin_labels and plot.ptype == Types.ratio:
            plot.set_bin_labels(h_sm)

        #if plot.rebin_bins:
        #    print "WARNING :: rebinning not yet implemented"
        #    #TODO: Implement rebinning

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
            nominalAsymErrorsNoSys.SetPointError(i-1,0,0,0,0)
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
        #self.ratio.SetLineWidth(2)
        #uglify
        self.ratio.SetLineWidth(1)
        self.ratio.SetMarkerStyle(20)
        self.ratio.SetMarkerSize(1.5)
        self.ratio.SetLineColor(1)
        self.ratio.SetMarkerSize(1.5)

        # Clean up
        ratio_raw.Delete()
        nominalAsymErrorsNoSys.Delete()
        g_data.Delete()

def make_stack_axis(plot):
    hax = r.TH1F("axes", "", int(plot.nbins), plot.xmin, plot.xmax)
    hax.SetMinimum(plot.ymin)
    hax.SetMaximum(plot.ymax)
    xax = hax.GetXaxis()
    xax.SetTitle(plot.xlabel)
    xax.SetTitleFont(42)
    xax.SetLabelFont(42)
    xax.SetLabelSize(0.035)
    xax.SetTitleSize(0.048 * 0.85)
    if plot.ptype == Types.ratio:
        hax.GetXaxis().SetTitleOffset(-999)
        hax.GetXaxis().SetLabelOffset(-999)
    else:
        xax.SetLabelOffset(1.15 * 0.02)
        xax.SetTitleOffset(1.5 * xax.GetTitleOffset())

    yax = hax.GetYaxis()
    yax.SetTitle(plot.ylabel)
    yax.SetTitleFont(42)
    yax.SetLabelFont(42)
    yax.SetTitleOffset(1.4)
    yax.SetLabelOffset(0.013)
    yax.SetLabelSize(1.2 * 0.035)
    yax.SetTitleSize(0.055 * 0.85)

    if plot.bin_labels and plot.ptype == Types.stack:
        plot.set_bin_labels(hax)
    if plot.rebin_bins:
        new_bins = array('d', plot.rebin_bins)
        hax_rebinned = hax.Rebin(len(new_bins)-1, 'axes_rebinned', new_bins)
        hax.Delete()
        hax = hax_rebinned
    return hax

class DataMCStackHist1D(HistBase):
    def __init__(self, plot, reg, YIELD_TBL, data, bkgds, sig=None):
        # Get plotting primitives
        # Not all primitives are for drawing. Some are for preserving pointers
        #TOOD: Remove YIELD_TBL after validating
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

        self.make_stack_legend(plot, reg)
        self.axis = make_stack_axis(plot)
        self.add_stack_backgrounds(plot, reg, YIELD_TBL, bkgds)
        #self.add_stack_signals(plot, reg, sigs)
        self.add_stack_data(plot, reg, YIELD_TBL, data)
        self.add_stack_mc_errors(plot, reg)
        if plot.doNorm:
            self.normalize_stack(plot, reg)
        if plot.auto_set_ylimits:
            self.reformat_axis(plot, reg)

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

    def make_stack_legend(self, plot, reg):
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
        self.leg_sig = pu.default_legend(xl=0.55, yl=0.6, xh=0.91, yh=0.71)
        self.leg_sig.SetNColumns(1)

    def add_stack_backgrounds(self, plot, reg, YIELD_TBL, bkgds):

        # Initilize lists and defaults
        histos = []
        all_histos = []
        avoid_bkg = []
        hists_to_clear = []

        # Make MC sample hists
        for mc_sample in bkgds :
            # Initilize histogram
            var_name = pu.strip_for_root_name(plot.variable)
            h_name = "h_"+reg.name+'_'+mc_sample.name+"_"+ var_name
            hists_to_clear.append(h_name)
            h = pu.th1d(h_name, "", int(plot.nbins),
                        plot.xmin, plot.xmax,
                        plot.xlabel, plot.ylabel)

            h.SetLineColor(mc_sample.color)
            h.GetXaxis().SetLabelOffset(-999)
            if mc_sample.isSignal:
                h.SetLineWidth(2)
                h.SetLineStyle(2)
                h.SetFillStyle(0)
            else:
                h.SetFillColor(mc_sample.color)
                h.SetFillStyle(1001)
            h.Sumw2

            # Draw final histogram (i.e. selections and weights applied)
            if plot.variable != mc_sample.weight_str:
                weight_str = "%s * %s"%(mc_sample.weight_str, str(mc_sample.scale_factor))
            else:
                weight_str = 1
            cut = "(%s) * %s"%(reg.tcut, weight_str)
            cut = r.TCut(cut)
            sel = r.TCut("1")
            draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
            mc_sample.tree.Draw(draw_cmd, cut * sel, "goff")


            # Yield +/- stat error
            stat_err = r.Double(0.0)
            integral = h.IntegralAndError(0,-1,stat_err)

            # Rebin
            if plot.rebin_bins:
                new_bins = array('d', plot.rebin_bins)
                h = h.Rebin(len(new_bins)-1, h_name, new_bins)

            h.leg_name = mc_sample.displayname #dynamic class members...ooo yeah!

            # Add overflow
            if plot.add_overflow:
                pu.add_overflow_to_lastbin(h)
            if plot.add_underflow:
                pu.add_underflow_to_firstbin(h)

            # Record all and non-empty histograms
            if mc_sample.isSignal:
                leg_sig.AddEntry(h, mc_sample.displayname, "l")
                YIELD_TBL.signals[mc_sample.name] = UncFloat(integral, stat_err)
            else:
                all_histos.append(h)
                histos.append(h) if integral > 0 else avoid_bkg.append(mc_sample.name)
                YIELD_TBL.mc[mc_sample.name] = UncFloat(integral, stat_err)

        if not len(histos):
            print "ERROR (make_stack_background) :: All SM hists are empty. Skipping"
            return

        # Order the hists by total events
        histos = sorted(histos, key=lambda h: h.Integral())

        self.mc_stack = r.THStack("stack_"+plot.name, "")
        for h in histos :
            self.mc_stack.Add(h)

        h_leg = sorted(all_histos, key=lambda h: h.Integral(), reverse=True)
        self.histos_for_leg = self.arrange_histos_for_legend(h_leg)

        # draw the total bkg line
        self.mc_total = self.mc_stack.GetStack().Last().Clone("hist_sm")
        self.mc_total.SetLineColor(r.kBlack)
        self.mc_total.SetLineWidth(3)
        self.mc_total.SetLineStyle(1)
        self.mc_total.SetFillStyle(0)
        self.mc_total.SetLineWidth(3)

    def add_stack_data(self, plot, reg, YIELD_TBL, data):
        #TODO: Look for a way to combine with backgrounds
        var_name = pu.strip_for_root_name(plot.variable)
        hd_name = "h_"+reg.name+'_data_'+ var_name
        self.data_hist = pu.th1d(hd_name, "", int(plot.nbins),
                                  plot.xmin, plot.xmax,
                                  plot.xlabel, plot.ylabel)
        self.data_hist.Sumw2

        cut = "(" + reg.tcut + ")"
        cut = r.TCut(cut)
        draw_cmd = "%s>>%s"%(plot.variable, self.data_hist.GetName())

        data.tree.Draw(draw_cmd, cut, "goff")
        self.data_hist.GetXaxis().SetLabelOffset(-999)

        # print the yield +/- stat error
        stat_err = r.Double(0.0)
        integral = self.data_hist.IntegralAndError(0,-1,stat_err)
        YIELD_TBL.data[data.name] = UncFloat(integral, stat_err)

        # Rebin
        if plot.rebin_bins:
            new_bins = array('d', plot.rebin_bins)
            self.data_hist = self.data_hist.Rebin(len(new_bins)-1, hd_name, new_bins)

        # Add overflow
        if plot.add_overflow:
            pu.add_overflow_to_lastbin(self.data_hist)
        if plot.add_underflow:
            pu.add_underflow_to_firstbin(self.data_hist)

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

    def normalize_stack(self, plot, reg):
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

    def reformat_axis(self, plot, reg):
        ''' Reformat axis to fit content and labels'''
        if not self.mc_stack:
            return
        # Get maximum histogram y-value
        if self.data:
            maxy = max(pu.get_tgraph_max(self.data), self.mc_stack.GetMaximum())
            miny = min(pu.get_tgraph_min(self.data), self.mc_stack.GetMinimum())
        else:
            maxy = self.mc_stack.GetMaximum()
            miny = self.mc_stack.GetMinimum()
        if maxy <= 0:
            print "WARNING :: Max value of plot is <= 0"
            return

        # Get default y-axis max and min limits
        logy = plot.doLogY
        if logy:
            maxy = 10**(pu.get_order_of_mag(maxy))
            if miny > 0:
                miny = 10**(pu.get_order_of_mag(miny))
            else:
                miny = 10**(pu.get_order_of_mag(maxy) - 4)
        else:
            maxy = maxy
            miny = 0

        # Get y-axis max multiplier to fit labels
        if logy:
            max_mult = 1e6 if self.signals else 1e5
        else:
            max_mult = 2.0 if self.signals else 1.8

        # reformat the axis
        self.mc_stack.SetMaximum(max_mult*maxy)
        self.mc_stack.SetMinimum(miny)
        self.axis.SetMaximum(max_mult*maxy)
        self.axis.SetMinimum(miny)

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

class Hist2D(HistBase) :
    def __init__(self, plot, reg, YIELD_TBL, samples):
        #TODO: make simplest class take a single sample
        # and add derived class that takes lists of samples to add them
        # or write function that
        # TODO: Remove YIELD_TBL
        self.axis = None
        self.hist = None
        if not samples: return
        self.make_axis(plot)
        self.make_hist(plot, reg, YIELD_TBL, samples)

    def list_of_root_objects(self):
        return [self.axis, self.hist]

    def make_axis(self, plot):
        self.axis = r.TH2D("axes", "", plot.nxbins, plot.xmin, plot.xmax, plot.nybins, plot.ymin, plot.ymax)
        self.axis.SetMinimum(plot.zmin)
        self.axis.SetMaximum(plot.zmax)

        xax = self.axis.GetXaxis()
        xax.SetTitle(plot.xlabel)
        xax.SetTitleFont(42)
        xax.SetLabelFont(42)
        xax.SetLabelSize(0.035)
        xax.SetTitleSize(0.048 * 0.85)
        xax.SetLabelOffset(1.15 * 0.02)
        xax.SetTitleOffset(1.5 * xax.GetTitleOffset())

        #if plot.bin_labels:
        #    plot.set_bin_labels(self.axis)

        yax = self.axis.GetYaxis()
        yax.SetTitle(plot.ylabel)
        yax.SetTitleFont(42)
        yax.SetLabelFont(42)
        yax.SetTitleOffset(1.4)
        yax.SetLabelOffset(0.013)
        yax.SetLabelSize(1.2 * 0.035)
        yax.SetTitleSize(0.055 * 0.85)

        zax = self.axis.GetZaxis()
        zax.SetTitle(plot.zlabel)
        #zax.SetTitleFont(42)
        #zax.SetLabelFont(42)
        #zax.SetTitleOffset(1.4)
        #zax.SetLabelOffset(0.013)
        #zax.SetLabelSize(1.2 * 0.035)
        #zax.SetTitleSize(0.055 * 0.85)

        #if plot.bin_labels:
        #    plot.set_ybin_labels(self.axis)

        #if plot.rebin_xbins:
        #    new_bins = array('d', plot.rebin_xbins)
        #    self.axis = self.axis.RebinX(len(new_bins)-1, 'axes', new_bins)
        #if plot.rebin_ybins:
        #    new_bins = array('d', plot.rebin_ybins)
        #    self.axis = self.axis.RebinY(len(new_bins)-1, 'axes', new_bins)

    def make_hist(self, plot, reg, YIELD_TBL, samples):
        self.hist = r.TH2D(plot.name, "", plot.nxbins, plot.xmin, plot.xmax, plot.nybins, plot.ymin, plot.ymax)
        for sample in samples:

            x_var = pu.strip_for_root_name(plot.xvariable)
            y_var = pu.strip_for_root_name(plot.yvariable)
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

            # Yield +/- stat error
            stat_err = r.Double(0.0)
            integral = h_tmp.IntegralAndError(0,-1,0,-1,stat_err)
            if sample.isMC and sample.isSignal:
                YIELD_TBL.signals[sample.name] = UncFloat(integral, stat_err)
            elif sample.isMC and not sample.isSignal:
                YIELD_TBL.mc[sample.name] = UncFloat(integral, stat_err)
            elif not sample.isMC:
                YIELD_TBL.data[sample.name] = UncFloat(integral, stat_err)

            self.hist.Add(h_tmp)

        zax = self.hist.GetZaxis()
        zax.SetTitle(plot.zlabel)
        zax.SetTitleFont(42)
        zax.SetLabelFont(42)
        zax.SetTitleOffset(1.5)
        zax.SetLabelOffset(0.013)
        zax.SetLabelSize(1.2 * 0.035)
