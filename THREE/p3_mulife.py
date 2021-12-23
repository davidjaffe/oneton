#!/usr/bin/env python
'''
fit 1d hists produced by Rong to measure muon lifetime 
20211221
'''
import numpy
#import graphUtils
import ROOT
from ROOT import TFile,TH1D,TH2D,gDirectory
#import gfit
#import pipath
import os
import datetime
import sys
import math

import Logger

class p3_mulife():
    def __init__(self,inputRFN=['DECAY_DATA/mc_decay_photontime.root'],useToffset=False):
        if inputRFN is None:
            sys.exit( 'p3_mulife.__init__ NO INPUT ROOT FILE LIST SPECIFIED')
        self.rfn = inputRFN

        self.useToffset = useToffset


        now = datetime.datetime.now()
        self.now = now.strftime('%Y%m%dT%H%M%S')

        parentDir = 'JOBS/'+self.now
        dirs = [parentDir]
        self.Figures = parentDir  + '/FIGURES/'
        self.logDir = parentDir
        dirs.append( self.Figures)
        dirs.append( self.logDir)

        for d in dirs:
            if not os.path.exists(d):
                os.makedirs(d)
                print('p3_mulife.__init__ create directory',d)

        lf = self.logDir + '/logfile.log'
        sys.stdout = Logger.Logger(fn=lf)
        print('p3_mulife.__init__ Output directed to stdout and',lf)

        
        self.Hists = {}

        
        print('p3_mulife.__init__ ROOT file: input',', '.join(self.rfn))
        print('p3_mulife.__init__ useToffset',useToffset)

        return
    def makeHists(self,func,nx=47,xlo=600.,xhi=2400.,Nhist=10,Nevt=1000,tau=2000.):
        '''
        make Nhists 1d histograms, range (xlo,xhi), with Nevt events each, 
        randomly drawn from function func with lifetime parameter tau
        '''
        self.Hists = {}
        self.initLife(func,A=1.,tau=tau)
        for ih in range(Nhist):
            hname = 'rh'+f'{ih:05}'
            title = hname +'_tau_'+f'{tau:.2f}'
            self.Hists[hname] = ROOT.TH1D(hname,title,nx,xlo,xhi)
            for i in range(Nevt):
                t = func.GetRandom()
                self.Hists[hname].Fill(t)
        print('p3_mulife.makeHists Created',Nhist,'hists with',Nevt,'events each, in range',xlo,'to',xhi,'ns, with tau(ns)=',tau)
        return
    def life(self,xlo=600.,xhi=2400.,name='life'):
        '''
        return TF1 object for lifetime fit
        '''
        func_form = '[0]*exp(-x/[1])'
        if self.useToffset:
            func_form = '[0]*exp(-(x-'+str(xlo)+')/[1])'
        print('p3_mulife.life functional form for fit is',func_form)
        func = ROOT.TF1(name,func_form,xlo,xhi,2)
        func.SetParName(0,'Const')
        func.SetParName(1,'Tau')
        self.initLife(func)
        return func
    def initLife(self,func,A=10.,tau=2000.):
        '''
        set initial parameter values for lifetime fit function
        '''
        func.SetParameter(0,A)
        func.SetParameter(1,tau)
        return
    def fitHist(self,func,histo,options='',makePDF=False,words=''):
        '''
        try to fit contents of histo for muon lifetime

        Try up to 3 times to get a good fit based on the fitResult.Status()
        '''
        xlo = histo.GetXaxis().GetXmin()
        xhi = histo.GetXaxis().GetXmax()

        if makePDF :
            c1 = self.makeCanvas()
            c1.Draw()
            
        fitResult = histo.Fit(func,options+'S')
        if fitResult.Status()!=0:
            stat1 = fitResult.Status()
            fitResult = histo.Fit(func,options+'S')
            stat2 = fitResult.Status()
            if stat2!=0 :
                fitResult = histo.Fit(func,options+'S')
                stat3 = fitResult.Status()
                if stat3!=0 : print('p3_mulife.fitHist fitResult status1',stat1,'status2',stat2,'status3',stat3,'fitResult.IsValid()',fitResult.IsValid(),'name',histo.GetName(),'words',words,fitResult)
        showFitResult = False
        if showFitResult : 
            print('p3_mulife.fitResult Correlation matrix and fit results follow')
            cormatrix = numpy.matrix([ [0.,0.], [0.,0.] ])
            covmatrix = numpy.matrix([ [0.,0.], [0.,0.] ])
            for i in range(2):
                for j in range(2):
                    cormatrix[i,j] = fitResult.Correlation(i,j)
                    covmatrix[i,j] = fitResult.CovMatrix(i,j)

            print('Correlation matrix',cormatrix,'\nCovariance matrix',covmatrix)
            print(fitResult)
        
        if makePDF :
            fname = self.Figures + histo.GetName() + '_' + words + '.pdf'
            c1.Print(fname)
            c1.IsA().Destructor(c1) # avoids seg fault?


        Const = func.GetParameter(0)
        Tau   = func.GetParameter(1)
        eConst= func.GetParError(0)
        eTau  = func.GetParError(1)
        fitProb=func.GetProb()
        return (Const,eConst),(Tau,eTau),fitProb
    def getHists(self,rfn,debug=False):
        '''
        fill dict with all hists from input root file
        return True if all objects inherit from TH1
        '''
        self.rf = rf = TFile.Open(rfn,"READ")
        keylist = rf.GetListOfKeys()
        AllTH1 = True
        
        for key in keylist:
            obj = key.ReadObj()

            if obj.IsA().InheritsFrom("TH1"):
                self.Hists[ obj.GetName() ] = obj
            else:
                AllTH1 = False
                print('this object',obj,'does not inherit from TH1')

        if debug : print('p3_mulife.getHists name, object for all TH1 found')
        for name in self.Hists:
            obj = self.Hists[name]
            if debug: print(name,obj)
        return  AllTH1
    def makeCanvas(self,name='c1',StatSize=None):
        '''
        return standard canvas
        StatSize controls size of text box
        '''
        c1 = ROOT.TCanvas(name)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(1111)
        ROOT.gStyle.SetTitleX(0.8)
        if StatSize is not None:
            ROOT.gStyle.SetStatW(StatSize) # size of stats box (and text?)
            ROOT.gStyle.SetStatH(StatSize)
        c1.SetGrid(1)
        c1.SetTicks(1)
        return c1
    def getHist(self,hn):
        hist = self.rf.Get(hn)
        return hist
    def fillDiagHists(self,fitResults,options,tauGen):
        '''
        fill diagnostic hists
        using fitResults[hname] = [ (A,eA), (tau,etau), prob] 
        tauGen = generated lifetime
        options = fit options
        '''
        words = options + '_tG' + str(int(tauGen)) # used to define output pdf
        
        A = [fitResults[hname][1][0]-tauGen for hname in fitResults]
        xma = max(abs(min(A)),max(A))
        xma = xma + 5.
        xmi = -xma
        diagHists = {}
        name  = 'fmg'
        title = 'Fitted - Generated lifetime (ns) '+words
        nx,xmi,xma = 100,xmi,xma
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)

        A =  [fitResults[hname][1][1] for hname in fitResults]
        xma = max(A) + 50.
        name  = 'efmg'
        title = 'uncertainty(Fitted - Generated lifetime) (ns) '+words
        nx,xmi,xma = 100,0.,xma
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
        
        name = 'fmgpull'
        title = 'Pull(fitted - generated lifetime) '+words
        nx,xmi,xma = 100,-10.,10.
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
        name = 'prob'
        title = 'Fit probability '+words
        nx,xmi,xma = 100,0.,1.
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)

        A = [fitResults[hname][2] for hname in fitResults]
        xmi = min(A)
        if xmi > 0. :
            xmi = math.log10(xmi) - 1.
        else:
            xmi = -30.
        xma = 0.
        name = 'log10prob'
        title = 'log10(Fit probability) '+words
        nx,xmi,xma = 100,xmi,xma
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
        

        for hname in fitResults:
            tau,etau = fitResults[hname][1]
            prob     = fitResults[hname][2]
            logprob = -30.+1.e-10
            if prob>0: logprob = math.log10(prob)
            diagHists['fmg'].Fill(tau-tauGen)
            diagHists['efmg'].Fill(etau)
            if etau>0: diagHists['fmgpull'].Fill((tau-tauGen)/etau)
            diagHists['prob'].Fill(prob)
            diagHists['log10prob'].Fill(logprob)

        histStats = {}
        for name in diagHists:
            histStats[name] = self.getStats(diagHists[name])            
            
            c1 = self.makeCanvas('cfD')
            ROOT.gStyle.SetOptStat(111111)
            ROOT.gStyle.SetOptFit(0)
            ROOT.gStyle.SetTitleX(0.5)
            c1.Draw()
            diagHists[name].Draw()
            fname = self.Figures + name + words +'.pdf'
            c1.Print(fname)
            c1.IsA().Destructor(c1) # avoids seg fault?
        return histStats
    def getStats(self,hist):
        '''
        return histogram stats
        '''
        entries = hist.GetEntries()
        mean    = hist.GetMean()
        emean   = hist.GetMeanError()
        return [entries, mean, emean]
    def plotFitStats(self,fitStats):
        '''
        plot fit statistics vs generated tau as extracted from histograms

        fitStats['OPT_tgXXXX'] = {'fmg' : [entries, mean, emean], 'efmg': [en,m,em] ...}
        fmg = fitted - generated
        efmg = uncertainty(fitted - generated)
        fmgpull = pull(fitted - generated)
        '''
        method= {'QL':'LogLike', 'Q':'Chisq'}
        qkeys = ['fmg', 'efmg', 'fmgpull']
        qnames= ['fitted - generated tau (ns)','uncertainty in fitted tau (ns)','pull(fitted - generated tau)']
        QvGen = {}
        for m in method:
            for q in qkeys:
                key = q+'_'+method[m]
                QvGen[key] = []
        for key in sorted(fitStats):
            ctauGen = key[-4:]
            tauGen = float(ctauGen)
            fitOpt = key.split('_')[0] # either 'Q' or 'QL'
            fitMethod = method[fitOpt]
            for qkey in qkeys:
                entries,mean,emean = fitStats[key][qkey]
                name = qkey+'_'+fitMethod
                QvGen[name].append( (tauGen,mean,emean) )

        graphs = {}
        for qkey in qkeys: graphs[qkey] = []

        colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kGreen, ROOT.kOrange+7, ROOT.kCyan, ROOT.kViolet-1]
        markers= [29,          20,         21,            22,          23,             33,         34]
        for iu,m in enumerate(method):
            for qkey,qtitle in zip(qkeys,qnames):
                name = qkey+'_'+method[m]
                title = qtitle + ' ' + method[m]
                A = sorted(QvGen[name])
                u = [x[0]+5.*float(iu) for x in A] # displace points slightly
                v = [x[1] for x in A]
                ev= [x[2] for x in A]
                graphs[qkey].append( self.makeTGraph(u,v,title,name,ev=ev) )
        for qkey,qtitle in zip(qkeys,qnames):
            tmg = self.makeTMultiGraph(qkey,tit='')
            lg  = ROOT.TLegend(.25,.85,.75,.99) # x1,y1,x2,y2 of corners
            for i,g in enumerate(graphs[qkey]):
                g.SetMarkerColor(colors[i])
                #g.SetLineColor(colors[i])
                g.SetMarkerStyle(markers[i])
                tmg.Add(g)
                lg.AddEntry(g,g.GetTitle(),"LP")
            canvas = ROOT.TCanvas(qkey)
            canvas.SetGrid(1)
            canvas.SetTicks(1)
            #print('p3_mulife.plotFitStats tmg.GetListOfGraphs()',tmg.GetListOfGraphs())
            tmg.Draw('ALP')
            tmg.GetXaxis().SetTitle('generated tau (ns)') ## this has to be after tmgDraw????
            tmg.GetYaxis().SetTitle(qtitle)
            gmi,gma = tmg.GetYaxis().GetXmin(),tmg.GetYaxis().GetXmax()
            dg = abs(gma-gmi)/5.
            gma += dg
            tmg.SetMaximum(gma) # increase maximum so that legend doesn't overlap points
            lg.Draw()
            canvas.Draw()
            canvas.cd()
            canvas.Modified() #needed?
            canvas.Update()   #needed?
            fname = self.Figures + qkey + '.pdf'
            canvas.Print(fname)
            canvas.IsA().Destructor(canvas) # needed?
        return
    def makeTGraph(self,u,v,title,name,eu=None,ev=None):
        '''
        make TGraph with axes u,v. eu,ev are uncertainties in u,v if not None
        '''

        dv = ev
        if ev is None: dv = [0. for x in range(len(u))]
        du = eu
        if eu is None: du = [0. for x in range(len(u))]
        g = ROOT.TGraphErrors(len(u),numpy.array(u),numpy.array(v),numpy.array(du),numpy.array(dv))
        g.SetTitle(title)
        g.SetName(name)
        if False : print('p3_mulife.makeTGraph Create obj',g,'name',name,'title',title,'len(u)',len(u),'u',u,'v',v,'du',du,'dv',dv)
        return g
    def makeTMultiGraph(self,name,tit=None,debug=False):
        title = tit
        if tit is None:title = name.replace('_',' ')
        tmg = ROOT.TMultiGraph()
        tmg.SetName(name)
        tmg.SetTitle(title)
        if debug:
            print('p3_mulife.makeTMultiGraph:name',name,'title',title,'object',tmg)
        return tmg

    def main(self,Nexpt=10,Nevt=1000):
        '''
        main module
        '''
        ROOT.gROOT.ProcessLine("gROOT->SetBatch()") # vague hope that this suppress TCanvas::Print messages
        
        func,nx,xlo,xhi = None,None,None,None

        GENERATE = False

        fitStats = {}
        for rfn in self.rfn:
            if 'generate' in rfn:
                GENERATE = True
                nx,xlo,xhi = 47,600.,2400.
                Nevt = Nevt
                Nhist = Nexpt
                thres = 1./float(Nhist)
                tauGen  = 2000.
                if 'tauGen' in rfn: tauGen = float(rfn.split()[2])
                func = self.life(xlo=xlo,xhi=xhi)
                self.makeHists(func,nx=nx,xlo=xlo,xhi=xhi,Nhist=Nhist,Nevt=Nevt,tau=tauGen)
            else:
                GENERATE = False
                words = ''
                OK = self.getHists(rfn,debug=True)
                if not OK : sys.exit('p3_mulife.main ERROR getHists return False for file '+self.rfn)

            for options in ['Q','QL']:
                if GENERATE : words = options + '_tG' + str(int(tauGen)) # used to define output pdf
                fitResults = {}
                for hname in self.Hists:
                    histo = self.Hists[hname]
                    if func is None :
                        xlo=histo.GetXaxis().GetXmin()
                        xhi=histo.GetXaxis().GetXmax()
                        nx =histo.GetNbins()
                        func = self.life(xlo=xlo,xhi=xhi)
                    y1 = histo.GetBinContent(1) ##### Add initialization of parameter A based on 1st bin in histogram
                    self.initLife(func,A=y1)
                    makePDF = numpy.random.random()<thres
                    fitResults[hname] = self.fitHist(func,histo,options=options,makePDF=makePDF,words=words)
                    if fitResults[hname][2]<1.e-10 :  # try refitting if prob(chi2)<.0001
                        prob_before = fitResults[hname][2]
                        fitResults[hname] = self.fitHist(func,histo,options=options,makePDF=True,words=words+'x2')
                        prob_after = fitResults[hname][2]
                        print('p3_mulife.main',hname,'prob(fit1)',prob_before,'prob(fit2)',prob_after)
                        
                        
                #print(fitResults)
                tau = numpy.array([fitResults[key][1][0] for key in fitResults])
                etau= numpy.array([fitResults[key][1][1] for key in fitResults])
                #print('tau',tau,'etau',etau)
                av = numpy.average(tau,weights=1/etau/etau)
                eav= numpy.sqrt(1./numpy.sum(1./etau/etau))
                print('fit options= {2}, weighted average {0:.2f}({1:.2f})'.format(av,eav,options))
                
                fitStats[words] = self.fillDiagHists(fitResults,options,tauGen)

            
            self.Hists = {}
            if not GENERATE : self.rf.Close()

        for key in sorted(fitStats):
            if False : print('p3_mulife.main key',key,'fitStats[key]',fitStats[key])
        self.plotFitStats(fitStats)
        return
if __name__ == '__main__' :

    inputRootFileNames = ['DECAY_DATA/mc_decay_photontime.root', 'DECAY_DATA/data_decay_photontime.root']

    taus = [1800., 1900., 2000., 2100., 2200.]
    inputRootFileNames = []
    for tauGen in taus:
        inputRootFileNames.append('generate tauGen '+str(tauGen))

    
    Nexpt = 10
    Nevt  = 1000
    useToffset = False
    if len(sys.argv)>1 :
        if 'help' in sys.argv[1].lower():
            print('USAGE: python p3_mulife.py Nexpt Nevt useToffset','\nDEFAULT: python p3_mulife.py',Nexpt,Nevt,useToffset)
            sys.exit('HELP WAS PROVIDED')
        Nexpt = int(sys.argv[1])
    if len(sys.argv)>2 : Nevt  = int(sys.argv[2])
    if len(sys.argv)>3 : useToffset = ('true' in sys.argv[3].lower()) 
        
    SC = p3_mulife(inputRFN=inputRootFileNames,useToffset=useToffset)

    SC.main(Nexpt=Nexpt,Nevt=Nevt)
    
