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
    def __init__(self,inputRFN=['DECAY_DATA/mc_decay_photontime.root']):
        if inputRFN is None:
            sys.exit( 'p3_mulife.__init__ NO INPUT ROOT FILE LIST SPECIFIED')
        self.rfn = inputRFN


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
        print('p3_mulife.makeHists Created',Nhist,'hists with',Nevt,'events each. tau(ns)=',tau)
        return
    def life(self,xlo=600.,xhi=2400.,name='life'):
        '''
        return TF1 object for lifetime fit
        '''
        func = ROOT.TF1(name,'[0]*exp(-x/[1])',xlo,xhi,2)
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
        '''
        xlo = histo.GetXaxis().GetXmin()
        xhi = histo.GetXaxis().GetXmax()

        if makePDF :
            c1 = self.makeCanvas()
            c1.Draw()
            
        histo.Fit(func,options)
        
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
        #print('p3_mulife.fillDiagHists',words,'tauGen',tauGen,'min(tauFit),max(tauFit)',min(A),max(A),'hist xmi,xma',xmi,xma)
        diagHists = {}
        name  = 'fmg'
        title = 'Fitted - Generated lifetime (ns) '+words
        nx,xmi,xma = 100,xmi,xma
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
        name = 'fmgpull'
        title = 'Pull(fitted - generated lifetime) '+words
        nx,xmi,xma = 100,-10.,10.
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)
        name = 'prob'
        title = 'Fit probability '+words
        nx,xmi,xma = 100,0.,1.
        diagHists[name] = ROOT.TH1D(name,title,nx,xmi,xma)

        for hname in fitResults:
            tau,etau = fitResults[hname][1]
            prob     = fitResults[hname][2]
            diagHists['fmg'].Fill(tau-tauGen)
            if etau>0: diagHists['fmgpull'].Fill((tau-tauGen)/etau)
            diagHists['prob'].Fill(prob)
        
        for name in diagHists:
            c1 = self.makeCanvas('cfD')
            ROOT.gStyle.SetOptStat(111111)
            ROOT.gStyle.SetOptFit(0)
            ROOT.gStyle.SetTitleX(0.5)
            c1.Draw()
            diagHists[name].Draw()
            fname = self.Figures + name + words +'.pdf'
            c1.Print(fname)
            c1.IsA().Destructor(c1) # avoids seg fault?
        return 
    def main(self,Nhist=10):
        '''
        main module
        '''
        #ROOT.gROOT.ProcessLine("gROOT->SetBatch()") # vague hope that this suppress TCanvas::Print messages
        
        func,nx,xlo,xhi = None,None,None,None

        GENERATE = False
        
        for rfn in self.rfn:
            if 'generate' in rfn:
                GENERATE = True
                nx,xlo,xhi = 47,600.,2400.
                Nevt = 1000
                Nhist = Nhist
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
            bias = {}    
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
                    self.initLife(func)
                    makePDF = numpy.random.random()<thres
                    fitResults[hname] = self.fitHist(func,histo,options=options,makePDF=makePDF,words=words)

                        
                #print(fitResults)
                tau = numpy.array([fitResults[key][1][0] for key in fitResults])
                etau= numpy.array([fitResults[key][1][1] for key in fitResults])
                #print('tau',tau,'etau',etau)
                av = numpy.average(tau,weights=1/etau/etau)
                eav= numpy.sqrt(1./numpy.sum(1./etau/etau))
                print('fit options= {2}, weighted average {0:.2f}({1:.2f})'.format(av,eav,options))
                self.fillDiagHists(fitResults,options,tauGen)


            
            self.Hists = {}
            if not GENERATE : self.rf.Close()
        return
if __name__ == '__main__' :
    #h2n = 'LED_WFD_At_vs_run_S0'



    inputRootFileNames = ['DECAY_DATA/mc_decay_photontime.root', 'DECAY_DATA/data_decay_photontime.root']

    inputRootFileNames = ['generate tauGen 1900.','generate tauGen 2000.','generate tauGen 2100.']
    Nhist = 10
    if len(sys.argv)>1 : Nhist = int(sys.argv[1])
    
    SC = p3_mulife(inputRFN=inputRootFileNames)

    SC.main(Nhist=Nhist)
    
