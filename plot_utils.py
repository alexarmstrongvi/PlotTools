import ROOT
from math import sqrt, log10, isnan
import re
ROOT.gROOT.SetBatch(True)
from array import array
from collections import defaultdict, OrderedDict
from tabulate import tabulate

ROOT.gStyle.SetCanvasPreferGL(ROOT.kTRUE)

ROOT.TH1F.__init__._creates = False
ROOT.TH2F.__init__._creates = False
ROOT.TCanvas.__init__._creates = False
ROOT.TPad.__init__._creates = False
ROOT.TLine.__init__._creates = False
ROOT.TLegend.__init__._creates = False
ROOT.TGraphErrors.__init__._creates = False
ROOT.TGraphAsymmErrors.__init__._creates = False
ROOT.TLatex.__init__._creates = False
ROOT.THStack.__init__._creates = False


# TODO Check how many of these methods I need
# many can probably be removed

def get_order_of_mag(num):
    return int(log10(num))

def normalize_plot(hist):
    if hist and hist.Integral():
        hist.Scale(1.0/hist.Integral())

def strip_for_root_name(string):
    string = re.sub(r'[(){}[\]\ ]+','', string)
    string = re.sub(r'-','_minus_', string)
    string = re.sub(r'/','_divide_', string)
    string = re.sub(r'\.','_', string)
    return string

def determine_bin_edges(lo, hi, nbins):
    assert lo < hi, ("ERROR (plot_utils.determine_bin_edges) :: lo > hi", lo, hi)
    step = (hi-lo)/float(nbins)
    return [lo + x*step for x in range(0, nbins+1)]

# ----------------------------------------------
#  TH1D Methods
# ----------------------------------------------
def th1d(name, title, nbin, nlow, nhigh, xtitle, ytitle) :
    '''
    Book a TH1F and return it
    '''
    h = ROOT.TH1D(name, title, nbin, nlow, nhigh)
    h.SetFillColorAlpha(ROOT.kRed, 20)
    font = 42
    h.SetTitleFont(font)

    # x-axis
    xaxis = h.GetXaxis()
    xaxis.SetTitle(xtitle)
    xaxis.SetTitleOffset(1.2 * xaxis.GetTitleOffset())
    xaxis.SetTitleFont(font)
    xaxis.SetLabelFont(font)

    # y-axis
    yaxis = h.GetYaxis()
    yaxis.SetTitle(ytitle)
    yaxis.SetTitleOffset(1.2 * yaxis.GetTitleOffset())
    yaxis.SetTitleFont(font)
    yaxis.SetLabelFont(font)

    # 2 is better than 1
    h.Sumw2()
    return h

def get_tgraph_y_values(tgraph):
    '''
        Get the y-values stored in a TGraph

        One cannot loop over TGraph.GetY() as it causes a crash
    '''
    y_values = []
    for ibin in range(tgraph.GetN()):
        y_values.append(tgraph.GetY()[ibin])
    return y_values

def get_tgraph_max(tgraph):
    '''
    Getting maximum y-value of TGraph because doesn't make things easy

    TGraphs cannot calculate their own maximum for reasons discuessed here:
        https://root-forum.cern.ch/t/tgraph-getmaximum-getminimum/8867

    '''
    maxy = tgraph.GetMaximum()
    if maxy == -1111: #default
        maxy = max(get_tgraph_y_values(tgraph))
    return maxy

def get_tgraph_min(tgraph):
    '''
    Getting minimum y-value of TGraph because doesn't make things easy

    TGraphs cannot calculate their own minimum for reasons discuessed here:
        https://root-forum.cern.ch/t/tgraph-getmaximum-getminimum/8867
    '''
    miny = tgraph.GetMinimum()
    if miny == -1111: #default
        miny = min(get_tgraph_y_values(tgraph))
    return miny


def draw_text_on_top(text="", size=0.04, pushright=1.0, pushup=1.0) :
    s = size
    t = text
    top_margin = ROOT.gPad.GetTopMargin()
    left_margin = ROOT.gPad.GetLeftMargin()
    xpos = pushright * left_margin
    ypos = 1.0 - 0.85*top_margin
    ypos *= pushup
    #draw_text(x=xpos, y=1.0-0.85*top_margin, text=t, size=s)
    draw_text(x=xpos, y=ypos, text=t, size=s)

def th1_to_tgraph(hist) :
    '''
    The provided histogram is turned :q
    into a TGraphErrors object
    '''

    g = ROOT.TGraphAsymmErrors()

    # don't care about the underflow/overflow
    for ibin in xrange(1,hist.GetNbinsX()+1) :
        y = hist.GetBinContent(ibin)
        ey = hist.GetBinError(ibin)
        x = hist.GetBinCenter(ibin)
        ex = hist.GetBinWidth(ibin) / 2.0
        if y == 0: y = -2
        g.SetPoint(ibin-1,x,y)
        g.SetPointError(ibin-1,ex,ex,ey,ey)

    return g

def convert_errors_to_poisson(hist, scale=False) :
    '''
    Provided a histogram, convert the errors
    to Poisson errors
    '''
    # needed variables
    alpha = 0.158655
    beta = 0.158655

    g = ROOT.TGraphAsymmErrors()

    for ibin in xrange(1,hist.GetNbinsX()+1) :
        value = hist.GetBinContent(ibin)
        width = hist.GetBinWidth(ibin)
        value_for_error = value
        if scale :
            print "ASSUMING VARIABLE BIN WIDTH NORMALIZATION IN DATA ERRORS"
            print "value = %.2f  width = %.2f  value*width = %.3f"%(value, width, value*width)

            value_for_error = value*width
        if value != 0 :
            error_poisson_up = 0.5 * ROOT.TMath.ChisquareQuantile(1-beta,2*(value_for_error+1))-value_for_error
            error_poisson_down = value_for_error - 0.5*ROOT.TMath.ChisquareQuantile(alpha,2*value_for_error)
            ex = hist.GetBinWidth(ibin) / 2.0
            ex = 0.0 # removing x error bins

            if scale :
                error_poisson_up = error_poisson_up / width
                error_poisson_down = error_poisson_down / width

            g.SetPoint(ibin-1, hist.GetBinCenter(ibin), value)
            g.SetPointError(ibin-1, ex, ex, error_poisson_down, error_poisson_up)
        else :
            g.SetPoint(ibin-1, hist.GetBinCenter(ibin), -2)
            g.SetPointError(ibin-1, 0., 0., 0., 0.)

    return g

def tgraphErrors_divide(g1, g2) :
    '''
    Provided two TGraphErrors objects, divide them
    and return the resulting TGraphErrors object
    '''
    n1 = g1.GetN()
    n2 = g2.GetN()
    if n1 != n2 :
        print "traphErrors_divide ERROR    input TGraphErrors do not have same number of entries!"
    g3 = ROOT.TGraphErrors()

    iv = 0
    for i1 in xrange(n1) :
        for i2 in xrange(n2) :
            x1 = ROOT.Double(0.0)
            y1 = ROOT.Double(0.0)
            x2 = ROOT.Double(0.0)
            y2 = ROOT.Double(0.0)
            dx1 = ROOT.Double(0.0)
            dy1 = ROOT.Double(0.0)
            dy2 = ROOT.Double(0.0)

            g1.GetPoint(i1,x1,y1)
            g2.GetPoint(i2,x2,y2)

            if x1 != x2 : continue
                #print "test"
            else :
                dx1 = g1.GetErrorX(i1)
                if y1 > 0 : dy1 = g1.GetErrorY(i1)/y1
                if y2 > 0 : dy2 = g2.GetErrorY(i2)/y2

                if y1 <= 0. or y2 <= 0 :
                    g3.SetPoint(iv, x1, -10) # if the ratio is zero, don't draw point at zero (looks bad on ratio pad)
                else:
                    g3.SetPoint(iv, x1, y1/y2)

            e = ROOT.Double(0.0)

            if y1 > 0 and y2 > 0 :
                e = sqrt(dy1*dy1 + dy2*dy2)*(y1/y2)
            g3.SetPointError(iv,dx1,e)

            iv += 1

    return g3

def tgraphAsymmErrors_divide(g_num, g_den) :
    n_num = g_num.GetN()
    n_den = g_den.GetN()
    if n_num != n_den :
        print "tgraphAsymmErrors_divide ERROR    input TGraphs do not have same number of entries!"
    g3 = ROOT.TGraphAsymmErrors()

    iv = 0
    for inum in xrange(n_num) :
        for iden in xrange(n_den) :
            x_num = ROOT.Double(0.0)
            y_num = ROOT.Double(0.0)
            x_den = ROOT.Double(0.0)
            y_den = ROOT.Double(0.0)


            ex = ROOT.Double(0.0)
            ey_num_up = ROOT.Double(0.0)
            ey_num_dn = ROOT.Double(0.0)
            ey_den_up = ROOT.Double(0.0)
            ey_den_dn = ROOT.Double(0.0)

            g_num.GetPoint(inum, x_num, y_num)
            g_den.GetPoint(iden, x_den, y_den)

            if x_num != x_den : continue

            if y_num >0 :
                ey_num_up = g_num.GetErrorYhigh(inum)/y_num
                ey_num_dn = g_num.GetErrorYlow(inum)/y_num
            if y_den >0 :
                ey_den_up = g_den.GetErrorYhigh(iden)/y_den
                ey_den_dn = g_den.GetErrorYlow(iden)/y_den

            if y_num <= 0. or y_den <= 0.:
                g3.SetPoint(iv, x_num, -10)
            else:
                g3.SetPoint(iv, x_num, y_num/y_den)
            ex = g_num.GetErrorX(iv)

            e_up = ROOT.Double(0.0)
            e_dn = ROOT.Double(0.0)
            if y_num > 0 and y_den > 0 :
                e_up = sqrt(ey_num_up*ey_num_up + ey_den_up*ey_den_up)*(y_num/y_den)
                e_dn = sqrt(ey_num_dn*ey_num_dn + ey_den_dn*ey_den_dn)*(y_num/y_den)
            g3.SetPointError(iv, ex, ex, e_dn, e_up)

            iv += 1
    return g3


def buildRatioErrorBand(g_in, g_out) :
    g_out.SetMarkerSize(0)
    for bin in xrange(g_out.GetN()) :
        y_out = ROOT.Double(1.0)
        x_out = ROOT.Double(0.0)
        y_in = ROOT.Double(0.0)
        x_in = ROOT.Double(0.0)

        g_in.GetPoint(bin, x_in, y_in)
        g_out.SetPoint(bin, x_in, y_out)

        # set upper error
        if y_in > 0.0001 :
            g_out.SetPointEYhigh(bin, g_in.GetErrorYhigh(bin)/y_in)
           # g_out.GetErrorYhigh(bin) = g_in.GetErrorYhigh(bin) / y_in
        else :
            g_out.SetPointEYhigh(bin, 0.0)
           # g_out.GetErrorYhigh(bin) = 0.0

        # set lower error
        if y_in > 0.0001 :
            g_out.SetPointEYlow(bin, g_in.GetErrorYlow(bin)/y_in)
            #g_out.GetErrorYow(bin) = g_in.GetErrorYlow(bin) / y_in
        else :
            g_out.SetPointEYlow(bin, 0.0)
            #g_out.GetErrorYlow(bin) = 0.0

        if g_out.GetErrorYlow(bin) > 1. :
            g_out.SetPointEYlow(bin, 1.0)
            #g_out.GetErrorYlow(bin) = 1.
        if g_out.GetErrorYhigh(bin) > 1. :
            g_out.SetPointEYhigh(bin, 1.0)
            #g_out.GetErrorYhigh(bin) = 1.


def add_overflow_to_lastbin(hist) :
    '''
    Given an input histogram, add the overflow
    to the last visible bin
    '''
    # find the last bin
    ilast = hist.GetXaxis().GetNbins()

    # read in the values
    lastBinValue = hist.GetBinContent(ilast)
    lastBinError = hist.GetBinError(ilast)
    overFlowValue = hist.GetBinContent(ilast+1)
    overFlowError = hist.GetBinError(ilast+1)

    # set the values
    hist.SetBinContent(ilast+1, 0)
    hist.SetBinError(ilast+1, 0)
    hist.SetBinContent(ilast, lastBinValue+overFlowValue)
    hist.SetBinError(ilast, sqrt(lastBinError*lastBinError + overFlowError*overFlowError))

def add_underflow_to_firstbin(hist) :
    '''
    Given an input histogram, add the underflow
    to the last visible bin
    '''
    # find the last bin
    ifirst = 1

    # read in the values
    firstBinValue = hist.GetBinContent(ifirst)
    firstBinError = hist.GetBinError(ifirst)
    underFlowValue = hist.GetBinContent(ifirst-1)
    underFlowError = hist.GetBinError(ifirst-1)

    # set the values
    hist.SetBinContent(ifirst-1, 0)
    hist.SetBinError(ifirst-1, 0)
    hist.SetBinContent(ifirst, firstBinValue+underFlowValue)
    hist.SetBinError(ifirst, sqrt(firstBinError**2 + underFlowError**2))

def divide_histograms(hnum, hden, xtitle, ytitle) :
    '''
    Provide two histograms and divide hnum/hden.
    Converts the final result into a tgraph.
    '''
    nbins = hnum.GetNbinsX()
    xlow = hnum.GetBinCenter(1)
    xhigh = hnum.GetBinCenter(nbins+1)
    hratio = hnum.Clone("hratio")
    hratio.GetYaxis().SetTitle(ytitle)
    hratio.GetXaxis().SetTitle(xtitle)

    hratio.GetYaxis().SetTitleOffset(0.45 * hratio.GetYaxis().GetTitleOffset())
    hratio.GetYaxis().SetLabelSize(2 * hratio.GetYaxis().GetLabelSize())
    hratio.GetYaxis().SetTitleFont(42)
    hratio.GetYaxis().SetTitleSize(0.09)

    hratio.GetXaxis().SetTitleOffset(1.75 * hratio.GetYaxis().GetTitleOffset())
    hratio.GetXaxis().SetLabelSize(3 * hratio.GetXaxis().GetLabelSize())
    hratio.GetXaxis().SetTitleSize(0.15)
    hratio.GetXaxis().SetTitleFont(42)

    g = ROOT.TGraphAsymmErrors()
    gdata = th1_to_tgraph(hnum)
    for i in range(1,nbins+1) :
        c1 = float(hnum.GetBinContent(i))
        c2 = float(hden.GetBinContent(i))
        if c2 == 0 : continue
        c3 = c1 / c2 * 1.0
        if c3 == 0 : c3 = -99
        hratio.SetBinContent(i, c3)

        g.SetPoint(i-1, hnum.GetBinCenter(i), c3)
        error_up = ( ( c2 + gdata.GetErrorYhigh(i-1) ) / c2 ) - 1.0
        error_dn = 1.0 - ( ( c2 - gdata.GetErrorYlow(i-1)) / c2 )
        g.SetPointError(i-1, 0.5 * hratio.GetBinWidth(i), 0.5 * hratio.GetBinWidth(i), error_dn, error_up)
    return g

def add_to_band(g1, g2) : #, sys_name) :
    if g1.GetN()!=g2.GetN() :
        print "plot_utils::add_to_band WARNING    input graphs do not have the same number of points!"

  #  eyhigh = ROOT.Double(0.0)
  #  eylow  = ROOT.Double(0.0)

  #  x1 = ROOT.Double(0.0)
  #  y1 = ROOT.Double(0.0)
  #  y2 = ROOT.Double(0.0)
  #  y0 = ROOT.Double(0.0)

    for i in xrange(g1.GetN()) :
        eyhigh = ROOT.Double(0.0)
        eylow  = ROOT.Double(0.0)
        eyhigh = g2.GetErrorYhigh(i)
        eylow  = g2.GetErrorYlow(i)

        x1 = ROOT.Double(0.0)
        y1 = ROOT.Double(0.0)
        y2 = ROOT.Double(0.0)
        y0 = ROOT.Double(0.0)

        g1.GetPoint(i, x1, y1)
        g2.GetPoint(i, x1, y2)

        if y1 == 0 : y1 = 1
        if y2 == 0 : y2 = 1

        eyh = ROOT.Double(0.0)
        eyl = ROOT.Double(0.0)

        y0 = y1 - y2
        #print "    > y0 : ", y0
        if y0 != 0 :
            if y0 > 0 :
                eyh = eyhigh
                eyh = sqrt(eyh*eyh + y0*y0)
                #print "    > %s + "%sys_name, eyh
                g2.SetPointEYhigh(i,eyh)
            else :
                eyl = eylow
                eyl = sqrt(eyl*eyl + y0*y0)
                #print "    > %s - "%sys_name, eyl
                g2.SetPointEYlow(i,eyl)

def get_bin_edges(axis):
    return [axis.GetBinLowEdge(x) for x in range(1, axis.GetNbins()+2)]
        
         
# ----------------------------------------------
#  TH2F Methods
# ----------------------------------------------
def th2f(name, title, nxbin, xlow, xhigh, nybin, ylow, yhigh, xtitle, ytitle) :
    '''
    Book a TH2F and return it
    '''
    h = ROOT.TH2F(name, title, nxbin, xlow, xhigh, nybin, ylow, yhigh)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.Sumw2()
    return h

def make_rebinned_th2f(hist, xbins=None, ybins=None, name=''):
    #NOTE: overwrites hist if name is the same as hist.GetName() or name is not provided
    # else it returns a new histogram
    if not xbins:
        xbins = array('d',get_bin_edges(hist.GetXaxis()))
    if not ybins:
        ybins = array('d', get_bin_edges(hist.GetYaxis()))
    nx = len(xbins) - 1
    ny = len(ybins) - 1
    old_name = name if name else hist.GetName()
    name = name + 'tmp'
    hnew = ROOT.TH2F(name, hist.GetTitle(), nx,xbins,ny,ybins)
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    for ibin in range(0, xaxis.GetNbins()+1):
        for jbin in range(0, yaxis.GetNbins()+1):
            xval = xaxis.GetBinCenter(ibin)
            yval = yaxis.GetBinCenter(jbin)
            bin_val = hist.GetBinContent(ibin, jbin)
            hnew.Fill(xval ,yval, bin_val)

    if old_name == hist.GetName() or not name:
        hist.Delete()
    hnew.SetName(old_name)
    return hnew

def make_rebinned_th3(hist, xbins=None, ybins=None, zbins=None, name=''):
    #NOTE: overwrites hist if name is the same as hist.GetName() or name is not provided
    # else it returns a new histogram
    if not xbins:
        xbins = array('d',get_bin_edges(hist.GetXaxis()))
    if not ybins:
        ybins = array('d', get_bin_edges(hist.GetYaxis()))
    if not zbins:
        zbins = array('d', get_bin_edges(hist.GetZaxis()))
    nx = len(xbins) - 1
    ny = len(ybins) - 1
    nz = len(zbins) - 1
    old_name = name if name else hist.GetName()
    name = name + 'tmp'
    hnew = ROOT.TH3F(name, hist.GetTitle(), nx,xbins,ny,ybins,nz,zbins)
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    zaxis = hist.GetZaxis()
    for ibin in range(0, xaxis.GetNbins()+1):
        for jbin in range(0, yaxis.GetNbins()+1):
            for kbin in range(0, zaxis.GetNbins()+1):
                xval = xaxis.GetBinCenter(ibin)
                yval = yaxis.GetBinCenter(jbin)
                zval = zaxis.GetBinCenter(kbin)
                bin_val = hist.GetBinContent(ibin, jbin, kbin)
                hnew.Fill(xval, yval, zval, bin_val)

    if old_name == hist.GetName() or not name:
        hist.Delete()
    hnew.SetName(old_name)
    return hnew

def hist_to_dict(hist, add_overflow=False, bin_range_keys=False):
    # Returns a dict of the form dict[xbin][ybin] = bin value
    # Options to include overflow and underflow bins or use
    # bin range (e.g. "1.30 - 2.41"). Precision is currently hard coded 
    if isinstance(hist, ROOT.TH3):
        h_dict = OrderedDict()
        xaxis = hist.GetXaxis()
        yaxis = hist.GetYaxis()
        zaxis = hist.GetZaxis()
        nxbins = xaxis.GetNbins()+1
        nybins = yaxis.GetNbins()+1
        nzbins = zaxis.GetNbins()+1
        for xbin in range(0, nxbins + 1):
            for ybin in range(0, nybins + 1):
                for zbin in range(0, nzbins + 1):
                    if bin_range_keys:
                        xbin_low = xaxis.GetBinLowEdge(xbin)
                        xbin_high = xaxis.GetBinUpEdge(xbin)
                        ybin_low = yaxis.GetBinLowEdge(ybin)
                        ybin_high = yaxis.GetBinUpEdge(ybin)
                        zbin_low = zaxis.GetBinLowEdge(zbin)
                        zbin_high = zaxis.GetBinUpEdge(zbin)
                        if xbin == nxbins:
                            xkey = "> %.2f" % xbin_low 
                        elif xbin == 0:
                            xkey = "< %.2f" % xbin_high 
                        else:
                            xkey = "%.2f-%.2f" % (xbin_low, xbin_high) 
                        if ybin == nybins:
                            ykey = "> %.2f" % ybin_low 
                        elif ybin == 0:
                            ykey = "< %.2f" % ybin_high 
                        else:
                            ykey = "%.2f-%.2f" % (ybin_low, ybin_high) 
                        if zbin == nzbins:
                            zkey = "> %.2f" % zbin_low 
                        elif zbin == 0:
                            zkey = "< %.2f" % zbin_high 
                        else:
                            zkey = "%.2f-%.2f" % (zbin_low, zbin_high) 
                    else:
                        xkey = str(xbin)
                        ykey = str(ybin)
                        zkey = str(zbin)
                    bin_value = hist.GetBinContent(xbin, ybin, zbin)
                    underflow = (xbin == 0) or (ybin == 0) or (zbin == 0)
                    overflow = (xbin == nxbins) or (ybin == nybins) or (zbin == nzbins)
                    if (underflow or overflow) and not add_overflow: continue
                    if zkey not in h_dict: h_dict[zkey] = OrderedDict()
                    if ykey not in h_dict[zkey]: h_dict[zkey][ykey] = OrderedDict()
                    h_dict[zkey][ykey][xkey] = bin_value
    elif isinstance(hist, ROOT.TH2):
        h_dict = OrderedDict()
        xaxis = hist.GetXaxis()
        yaxis = hist.GetYaxis()
        nxbins = xaxis.GetNbins()+1
        nybins = yaxis.GetNbins()+1
        for xbin in range(0, nxbins + 1):
            for ybin in range(0, nybins + 1):
                if bin_range_keys:
                    xbin_low = xaxis.GetBinLowEdge(xbin)
                    xbin_high = xaxis.GetBinUpEdge(xbin)
                    ybin_low = yaxis.GetBinLowEdge(ybin)
                    ybin_high = yaxis.GetBinUpEdge(ybin)
                    if xbin == nxbins:
                        xkey = "> %.2f" % xbin_low 
                    elif xbin == 0:
                        xkey = "< %.2f" % xbin_high 
                    else:
                        xkey = "%.2f-%.2f" % (xbin_low, xbin_high) 
                    if ybin == nybins:
                        ykey = "> %.2f" % ybin_low 
                    elif ybin == 0:
                        ykey = "< %.2f" % ybin_high 
                    else:
                        ykey = "%.2f-%.2f" % (ybin_low, ybin_high) 
                else:
                    xkey = str(xbin)
                    ykey = str(ybin)
                bin_value = hist.GetBinContent(xbin, ybin)
                underflow = (xbin == 0) or (ybin == 0)
                overflow = (xbin == nxbins) or (ybin == nybins)
                if (underflow or overflow) and not add_overflow: continue
                if xkey not in h_dict: h_dict[xkey] = OrderedDict()
                h_dict[xkey][ykey] = bin_value
    elif isinstance(hist, ROOT.TH1):
        h_dict = OrderedDict()
        xaxis = hist.GetXaxis()
        nxbins = xaxis.GetNbins()+1
        for xbin in range(0, nxbins + 1):
            underflow = (xbin == 0)
            overflow = (xbin == nxbins)
            if bin_range_keys:
                xbin_low = xaxis.GetBinLowEdge(xbin)
                xbin_high = xaxis.GetBinUpEdge(xbin)
                if overflow:
                    key = "> %.2f" % xbin_low 
                elif underflow:
                    key = "< %.2f" % xbin_high 
                else:
                    key = "%.2f-%.2f" % (xbin_low, xbin_high)
            else:
                key = str(xbin)
            bin_value = hist.GetBinContent(xbin)
            if (underflow or overflow) and not add_overflow: continue
            h_dict[key] = bin_value
    else:
        print "Histogram type not recognized: ", type(hist)
    return h_dict

def print_hist(hist, tablefmt='psql'):
    hist_dict_bins = hist_to_dict(hist, add_overflow=True, bin_range_keys=False)
    hist_dict_range = hist_to_dict(hist, add_overflow=True, bin_range_keys=True)
    if not hist_dict_bins or not hist_dict_range: return ""

    if isinstance(hist, ROOT.TH3):
        tex_string = ""
        for z_bin, z_range in zip(hist_dict_bins, hist_dict_range):
            headers = [""]
            table = []
            if tablefmt == 'latex':
                tex_string += "(%s) %s\\\\\n" % (z_bin, z_range)
            else:
                tex_string += "(%s) %s\n" % (z_bin, z_range)
            # Fill headers and row lables
            if not table:
                for y_bin, y_range in zip(hist_dict_bins[z_bin], hist_dict_range[z_range]):
                    table.append(["(%s) %s" % (y_bin, y_range)])
                    for x_bin, x_range in zip(hist_dict_bins[z_bin][y_bin], hist_dict_range[z_range][y_range]):
                        header = "(%s) %s" % (x_bin, x_range)
                        if header not in headers: headers.append(header)

            # Fill table with values
            for row, (y_bin, y_range) in enumerate(zip(hist_dict_bins[z_bin], hist_dict_range[z_range])):
                for x_bin, x_range in zip(hist_dict_bins[z_bin][y_bin], hist_dict_range[z_range][y_range]):
                    table[row].append(hist_dict_bins[z_bin][y_bin][x_bin])
            tex_string += tabulate(table, headers=headers, tablefmt=tablefmt)        
            tex_string += '\n\n'
    elif isinstance(hist, ROOT.TH2):
        headers = [""]
        table = []
        for x_bin, x_range in zip(hist_dict_bins, hist_dict_range):
            headers.append("(%s) %s" % (x_bin, x_range))
            if not table:
                # Fill first column of table with y-labels
                for y_bin, y_range in zip(hist_dict_bins[x_bin], hist_dict_range[x_range]):
                    table.append(["(%s) %s" % (y_bin, y_range)])
            for row, (y_bin, y_range) in enumerate(zip(hist_dict_bins[x_bin], hist_dict_range[x_range])):
                table[row].append(hist_dict_bins[x_bin][y_bin])
        tex_string = tabulate(table, headers=headers, tablefmt=tablefmt)        
    elif isinstance(hist, ROOT.TH1):
        headers = []
        row = []
        for x_bin, x_range in zip(hist_dict_bins, hist_dict_range):
            headers.append("(%s) %s" % (x_bin, x_range))
            row.append(hist_dict_bins[x_bin])
        tex_string = tabulate([row], headers=headers, tablefmt=tablefmt)        
    return tex_string


def print_hist_old(hist):
    print_str = ""
    if isinstance(hist, ROOT.TH2):
        xaxis = hist.GetXaxis()
        yaxis = hist.GetYaxis()
        nxbins = xaxis.GetNbins()+1
        nybins = yaxis.GetNbins()+1
        print_str += "\t(xbin, ybin) = (xbin range, ybin range) = bin value\n"
        for xbin in range(0, nxbins + 1):
            for ybin in range(0, nybins + 1):
                #xbin_center = xaxis.GetBinCenter(xbin)
                #ybin_center = yaxis.GetBinCenter(ybin)
                xbin_low = xaxis.GetBinLowEdge(xbin)
                xbin_high = xaxis.GetBinUpEdge(xbin)
                ybin_low = yaxis.GetBinLowEdge(ybin)
                ybin_high = yaxis.GetBinUpEdge(ybin)
                bin_value = hist.GetBinContent(xbin, ybin)
                underflow = xbin == 0 or ybin == 0
                overflow = xbin == nxbins or ybin == nybins
                if (underflow or overflow) and not bin_value: continue
                print_str += "\t(%d, %d) = ([%.2f - %.2f], [%.2f - %.2f]) = %.3f\n" % (
                        xbin, ybin, xbin_low, xbin_high, ybin_low, ybin_high, bin_value)
    elif isinstance(hist, ROOT.TH1):
        xaxis = hist.GetXaxis()
        nxbins = xaxis.GetNbins()+1
        print_str += "\t(xbin) = [xbin range] = bin value\n"
        for xbin in range(0, nxbins + 1):
                xbin_low = xaxis.GetBinLowEdge(xbin)
                xbin_high = xaxis.GetBinUpEdge(xbin)
                bin_value = hist.GetBinContent(xbin)
                underflow = xbin == 0
                overflow = xbin == nxbins
                if (underflow or overflow) and not bin_value: continue
                print_str += "\t(%d) = [%.2f - %.2f] = %.3f\n" % (
                        xbin, xbin_low, xbin_high, bin_value)
    else:
        print "Histogram type not recognized: ", type(hist)
    return print_str

def scale_thstack(stack, scale_factor):
    '''
    Scale each hist in a thstack

    It is important to call Modified() to overwite the reset the stack object.
    Otherwise, the scaling will not be present when drawing the stack.
    '''
    for hist in stack.GetHists():
        hist.Scale(scale_factor)
    stack.Modified()

def scale_tgraph(tgraph, scale_factor):
    for ii in range(tgraph.GetN()):
        if isnan(tgraph.GetY()[ii]): continue
        tgraph.GetY()[ii] *= scale_factor
        if isinstance(tgraph, ROOT.TGraphErrors):
            tgraph.GetEY()[ii] *= scale_factor
        elif isinstance(tgraph, ROOT.TGraphAsymmErrors):
            tgraph.GetEYhigh()[ii] *= scale_factor
            tgraph.GetEYlow()[ii] *= scale_factor
        elif isinstance(tgraph, ROOT.TGraphBentErrors):
            print "WARNING :: No scaling implemented yet"

# ----------------------------------------------
#  TLegend Methods
# ----------------------------------------------
def default_legend(xl=0.7,yl=0.75,xh=0.9,yh=0.88) :
    leg = ROOT.TLegend(xl, yl, xh, yh)
   # leg.SetNDC()
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    return leg

# ----------------------------------------------
#  Text/Label Methods
# ----------------------------------------------
def draw_text(x=0.7, y=0.65, font=42, color=ROOT.kBlack, text="", size=0.04, angle=0.0) :
    '''
    Draw text on the current pad
    (coordinates are normalized to the current pad)
    '''
    l = ROOT.TLatex()
    l.SetTextSize(size)
    l.SetTextFont(font)
    l.SetNDC()
    l.SetTextColor(color)
    l.SetTextAngle(angle)
    l.DrawLatex(x, y, text)

def draw_atlas_label(status, analysis, region, move_x = 0, move_y = 0, scale = 1):
    left_edge = 0.18 + move_x
    status_indent = (0.13) * scale
    bottom_edge = (0.73 + move_y)
    vspacing = 0.05 * scale
    size = 0.04 * scale

    draw_text(text="ATLAS",
        x=left_edge,               y=bottom_edge+3*vspacing, size=size+0.01, font=72)
    draw_text(text=status,
        x=left_edge+status_indent, y=bottom_edge+3*vspacing, size=size+0.01, font=42)
    draw_text(text="#sqrt{s} = 13 TeV, 36.1 fb^{-1}",
        x=left_edge,               y=bottom_edge+2*vspacing, size=size)
    draw_text(text=analysis,
        x=left_edge,               y=bottom_edge+vspacing,   size=size)
    draw_text(text=region,
        x=left_edge,               y=bottom_edge,            size=size)

# ----------------------------------------------
#  TLine Methods
# ----------------------------------------------
def draw_line(xl=0.0,yl=0.0,xh=1.0,yh=1.0,color=ROOT.kBlack,width=2,style=1) :
    l = ROOT.TLine(xl,yl,xh,yh)
    l.SetLineColor(color)
    l.SetLineWidth(width)
    l.SetLineStyle(style)
    l.Draw()

# ----------------------------------------------
#  Style Methods
# ----------------------------------------------
def set_palette(name="", ncontours=999) :
    if name == "gray" or name == "grayscale" :
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
#    elif name == "redbluevector" :
#        #stops = [0.00, 0.34, 0.61, 0.84, 1.00]
#        stops = [0.00, 0.17, 0.61, 0.84, 1.00]
#        #stops = [0.00, 0.20, 0.5, 0.70, 1.00]
#        red   = [0.35, 0.29, 0.29, 0.89, 0.90]
#        green = [0.70, 0.57, 0.35, 0.22, 0.05]
#        blue  = [0.95, 0.88, 0.70, 0.45, 0.09]
    elif name == "redbluevector" :
        stops = [0.09,      0.18,       0.27,     0.36,     0.45,     0.54,     0.63,     0.72,    0.81,    0.90,     1.00]
        red   = [1/256.,    2/256.,     2/256.,   58/256.,  90/256.,  115/256., 138/256., 158/256.,180/256.,200/256., 219/256.]
        green = [158/256.,  150/256.,   140/256., 131/256., 120/256., 110/256., 97/256.,  82/256., 62/256., 34/256.,  2/256.]
        blue  = [237/256.,  224/256.,   213/256., 200/256., 190/256., 178/256., 167/256., 156/256.,146/256.,134/256., 123/256.]
    elif name == "bluetowhite" :
        stops = [0.0, 0.125, 0.25, 0.375, 0.5, 0.66, 0.83, 1.00, 1.00]
        red   = [23/256.,  46/256.,  69/256.,  92/256., 124/256., 157/256., 190/256., 222/256.]
        green = [32/256.,  63/256.,  95/256., 126/256., 152/256., 178/256., 203/256., 229/256.]
        blue  = [57/256., 115/256., 172/256., 229/256., 235/256., 240/256., 245/256., 250/256.]
    else :
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
    s = array('d', stops)
    R = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, R, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)
