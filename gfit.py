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
    def fitNGaus(self,hname,nsigplus=1.5,nsigminus=2.0,debug=False,start_with_Chi2=False, inputPar=[None,0.1,None,None], inputLimits=[ [None,None], [None,None], [None,None], [None,None] ]):
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
        mupois = 0.1
        if inputPar[0] is not None: ent    = inputPar[0]
        if inputPar[1] is not None: mupois = inputPar[1]
        if inputPar[2] is not None: mean   = inputPar[2]
        if inputPar[3] is not None: sgm    = inputPar[3]
        for nsig in [300.]: #####[3., 2., -1.]:
            xlo = max(xmi,mean - nsig*sgm)
            xhi = min(xma,mean + nsig*sgm)
            if nsig<0:
                xlo = max(xmi,mean - nsigminus*sgm)
                xhi = min(xma,mean + nsigplus*sgm)
            if debug : print 'gfit.fitNGaus: prefit. mean,sgm, xlo,xhi',mean,sgm,xlo,xhi
            g2 = ROOT.TF1("g2",self.NGaus,xlo,xhi,4)

            g2.SetParameters( ent, mupois, mean, sgm)
            g2.SetParNames( 'Const', 'mu', 'mean', 'sigma')
            lo0,hi0  = inputLimits[0]
            if lo0 is None: lo0 = 0.
            if hi0 is None: hi0 = 10.*ent
            g2.SetParLimits(0, lo0, hi0)
            lo,hi = inputLimits[1]
            if lo is None: lo = 1.e-4
            if hi is None: hi = 1000.*mupois
            g2.SetParLimits(1, lo, hi)
            lo,hi = inputLimits[2]
            if lo is None: lo = 0.
            if hi is None: hi = 10.*ave
            g2.SetParLimits(2, lo, hi)
            lo3,hi3 =  inputLimits[3]
            if lo3 is None: lo3 = 1.e-4
            if hi3 is None: hi3 = 10.*rms
            g2.SetParLimits(3, lo3, hi3)
            
            if noPopUp : ROOT.gROOT.ProcessLine("gROOT->SetBatch()")

            g2.FixParameter(3, sgm) # fix width of gaussian for first pass
            g2.FixParameter(0, ent) # fix constant for first pass
            hname.Fit("g2",fit_options) # first pass
            
            g2.SetParLimits(0, lo0, hi0) # constraint on constant
            g2.SetParLimits(3, lo3, hi3) # constraint on width of gaussian
            hname.Fit("g2",fit_options)  # second pass
            
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
    def fakeNGaus(self,p,xlo=0.,xhi=500.):
        '''
        return histogram with a distribution given by NGaus parameters
        with limits xlo,xhi
        '''
        NG = ROOT.TF1("NG",self.NGaus,xlo,xhi,4)
        NG.SetParameters(p[0],p[1],p[2],p[3])
        title = 'Nevt {0} poisMu {1:.4f} gausMu {2:.2f} gausSg {3:.2f}'.format(int(p[0]),p[1],p[2],p[3])
        name = title.replace('.','_').replace(' ','_')
        H = ROOT.TH1D(name,title,500,xlo,xhi)
        H.FillRandom("NG",int(p[0]))
        H.Draw()
        return H
    def testFit(self,inputPoisMu=0.1,nevt=1000.):
        '''
        test fitting with fake NGaus distributions
        '''


        IC, Imupois, Imean, Isg1 = p = [nevt, inputPoisMu, 50., 25.] # C, poisMu, gausMu, gausSG
        IC, Imupois, Imean, Isg1 = p = [nevt, inputPoisMu, 50., 5.] # C, poisMu, gausMu, gausSG
        h = self.fakeNGaus(p)
        debug = False
        SwChi2= False
        GoodFit,mean,emean, sg1, prob,mupois = self.fitNGaus(h,debug=debug,start_with_Chi2=SwChi2)
        print 'INPUT:  mupois',Imupois,'mean',Imean,'sigma',Isg1
        print 'FREEFIT:mupois',mupois, 'mean', mean,'sigma',sg1
        if 0:
            GoodFit,mean,emean, sg1, prob,mupois = self.fitNGaus(h,inputPar = [None, Imupois, Imean, Isg1],debug=debug,start_with_Chi2=SwChi2)
            print 'Con1FIT:mupois',mupois, 'mean', mean,'sigma',sg1
            GoodFit,mean,emean, sg1, prob,mupois = self.fitNGaus(h,inputPar = [None, Imupois, Imean, None],debug=debug,start_with_Chi2=SwChi2)
            print 'Con2FIT:mupois',mupois, 'mean', mean,'sigma',sg1
            GoodFit,mean,emean, sg1, prob,mupois = self.fitNGaus(h,inputPar = [None, Imupois, None, None],debug=debug,start_with_Chi2=SwChi2)
            print 'Con3FIT:mupois',mupois, 'mean', mean,'sigma',sg1
            GoodFit,mean,emean, sg1, prob,mupois = self.fitNGaus(h,inputPar = [None, Imupois, None, Isg1],debug=debug,start_with_Chi2=SwChi2)
            print 'Con4FIT:mupois',mupois, 'mean', mean,'sigma',sg1
            
        inputLimits=[ [None,None], [1.e-4, 10.], [None,None], [.8*Isg1, 1.2*Isg1] ]
        GoodFit,mean,emean, sg1, prob,mupois = self.fitNGaus(h,inputPar = [None, Imupois, None, Isg1], inputLimits=inputLimits,debug=debug,start_with_Chi2=SwChi2)
        print 'Con5FIT:mupois',mupois, 'mean', mean,'sigma',sg1

        
        return h
if __name__ == '__main__' :
    G = gfit()
    hists = []
    ROOT.gROOT.ProcessLine("gROOT->SetBatch()")
    for nevt in [1000.]:
        for m in [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5., 6.]:
            G.Nterm = max(3,int(5.*m))
            h = G.testFit(inputPoisMu=m,nevt=nevt)
            print h.GetName()
            if h in hists: print 'appending duplicate hist',h.GetName()
            hists.append(h)
    rfn = 'testGfit.root'
    rf = ROOT.TFile(rfn,"RECREATE")
    for h in hists: rf.WriteTObject(h)
    rf.Close()
    print 'Wrote',len(hists),'hists to',rfn
    

