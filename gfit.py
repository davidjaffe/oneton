#!/usr/bin/env python
'''
gaussian fitter
'''

import ROOT
import math

class gfit():
    def __init__(self):
        print 'gfit initialized'
        self.sr2pi = math.sqrt(2.*math.pi)
        self.sr2   = math.sqrt(2.)
        self.Nterm = 3
        return

    def fit(self,hname,nsigplus=1.5,nsigminus=2.0,debug=False,start_with_Chi2=False):
        # iterative gaussian fit to hist hname
        # final fit is to range (nsigminus*sigma+mean,nsigplus*sigma+mean) with
        # log-likelihood where sigma and mean are results of previous iteration
        noPopUp = True
        
        name = hname.GetName()
        ave = hname.GetMean()
        rms = hname.GetRMS()
        tot = hname.GetEntries()
        under = hname.GetBinContent(0)
        bin1 = hname.GetBinContent(1)
        over  = hname.GetBinContent(hname.GetXaxis().GetNbins()+1)
        ent = tot - under - over 
        xmi = hname.GetXaxis().GetXmin()
        xma = hname.GetXaxis().GetXmax()

        GoodFit = False

        # reject hists with no entries or with too few entries
        # in the meaningful bins (>1 and <Nbin)
        if tot<=0 or float(ent)/float(tot)<0.50 :
            return GoodFit, ave,rms, -1., -1.

        if debug :
            print 'gfit.fit: name,tot,under,bin1,over,ent',name,tot,under,bin1,over,ent

        fit_options = "RLQ" # Restricted range,Likelihood,don't plot fit,quiet
        if start_with_Chi2:
            fit_options = "RQ" # Restricted range,chi2,don't plot fit,quiet
        nsig = 3.
        mean = ave
        sgm  = rms
        for nsig in [3., 2., -1]:
            xlo = max(xmi,mean - nsig*sgm)
            xhi = min(xma,mean + nsig*sgm)
            if nsig<0:
                xlo = max(xmi,mean - nsigminus*sgm)
                xhi = min(xma,mean + nsigplus*sgm)
            g1 = ROOT.TF1("g1","gaus",xlo,xhi)
            if noPopUp : ROOT.gROOT.ProcessLine("gROOT->SetBatch()")
            hname.Fit("g1",fit_options)
            fit_options = "RLQ" # Restricted range,Likelihood,don't plot fit,quiet
            c = g1.GetParameter(0)
            mean = g1.GetParameter(1)
            emean= g1.GetParError(1)
            sgm  = g1.GetParameter(2)
            prob = g1.GetProb()
            if debug : print 'gfit.fit: name,nsig,mean,emean,sgm,prob',name,nsig,mean,emean,sgm,prob
            
        emean = g1.GetParError(1)
        GoodFit = True
        return GoodFit,mean,emean, sgm, prob
    def fitNGaus(self,hname,nsigplus=1.5,nsigminus=2.0,debug=False,start_with_Chi2=False):
        '''
        return GoodFit (true/false), mean, emean, sgm1, prob of fit, mupois to hist hname with
        function NGaus
        NOTE NUMBER OF RETURNED VARIABLES DIFFERS FROM gfit.fit
        may include: 
        iterative  fit to hist hname
        final fit is to range (nsigminus*sigma+mean,nsigplus*sigma+mean) with
        log-likelihood where sigma and mean are results of previous iteration
        '''
        noPopUp = True
        
        name = hname.GetName()
        ave = hname.GetMean()
        rms = hname.GetRMS()
        tot = hname.GetEntries()
        under = hname.GetBinContent(0)
        bin1 = hname.GetBinContent(1)
        over  = hname.GetBinContent(hname.GetXaxis().GetNbins()+1)
        ent = tot - under - over 
        xmi = hname.GetXaxis().GetXmin()
        xma = hname.GetXaxis().GetXmax()

        GoodFit = False

        # reject hists with no entries or with too few entries
        # in the meaningful bins (>1 and <Nbin)
        if tot<=0 or float(ent)/float(tot)<0.50 :
            return GoodFit, ave,rms, -1., -1.

        if debug :
            print 'gfit.fitNGaus: name,tot,under,bin1,over,ent',name,tot,under,bin1,over,ent

        default_fit_opt = "RLQB" # restricted range, likelihood, quiet, bounded parameters
        fit_options = default_fit_opt # Restricted range,Likelihood,don't plot fit,quiet
        if start_with_Chi2:
            fit_options = default_fit_opt.replace("L","") # Restricted range,chi2,don't plot fit,quiet
        nsig = 300.
        mean = ave
        sgm  = rms
        for nsig in [300.]: #####[3., 2., -1.]:
            xlo = max(xmi,mean - nsig*sgm)
            xhi = min(xma,mean + nsig*sgm)
            if nsig<0:
                xlo = max(xmi,mean - nsigminus*sgm)
                xhi = min(xma,mean + nsigplus*sgm)
            if debug : print 'gfit.fitNGaus: prefit. mean,sgm, xlo,xhi',mean,sgm,xlo,xhi
            g2 = ROOT.TF1("g2",self.NGaus,xlo,xhi,4)

            g2.SetParameters( ent, 0.1, mean, sgm)
            g2.SetParNames( 'Const', 'mu', 'mean', 'sigma')
            g2.SetParLimits(0, 0., 10.*ent)
            g2.SetParLimits(1, 0., 5.)
            g2.SetParLimits(2, 0., 10.*ave)
            g2.SetParLimits(3, 0., 10.*rms)
            
            if noPopUp : ROOT.gROOT.ProcessLine("gROOT->SetBatch()")

            hname.Fit("g2",fit_options)
            fit_options = default_fit_opt
            c = g2.GetParameter(0)
            mupois = g2.GetParameter(1)
            mean = g2.GetParameter(2)
            emean= g2.GetParError(2)
            sg1  = g2.GetParameter(3)
            prob = g2.GetProb()
            if debug : print 'gfit.fitNGaus: name,nsig,mean,emean,sg1,prob,mupois',name,nsig,mean,emean,sg1,prob,mupois
            
        emean = g2.GetParError(2)
        GoodFit = True
        return GoodFit,mean,emean, sg1, prob,mupois
    def NGaus(self,v,p):
        '''
        4 parameter fit to sum of N gaussians
        f(x) = C * ( pois(1,mu) * G(mu1,sg1) + pois(2,mu) * G(mu2,sg2)  + ...)
        where mu2 = 2*mu1, sg2^2 = 2 * sg1^2
        '''
        x = v[0]
        
        C = p[0]
        mu = p[1]
        mu1 = p[2]
        sg1 = p[3]
        poisSum = 0.
        pois = math.exp(-mu)
        func = 0.
        for i in range(self.Nterm):
            n = float(i+1)
            pois *= mu/n
            poisSum += pois
            s = math.sqrt(n)*sg1
            a  = ( x - n*mu1)/s
            func += pois * math.exp( - a*a/2.) / s /self.sr2pi
        if poisSum==0.:
            func = 1.e20
            print 'gfit.NGaus: avoid div-by-zero x,C,mu,mu1,sg1,poisSum',x,C,mu,mu1,sg1,poisSum
        else:
            func = C * func / poisSum
        #x = v[0]
        #a1 = (x - mu1)/2./sg1
        #a2 = (x - mu2)/2./sg2
        #func = C * ( f * math.exp(-a1*a1)/sg1/self.sr2pi + (1.-f) * math.exp(-a2*a2)/sg2/self.sr2pi )
        return func
        
