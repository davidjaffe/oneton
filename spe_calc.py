#!/usr/bin/env python
'''
fit 1d projects to get spe vs run
20160208
'''
import numpy
import graphUtils
import ROOT
from ROOT import TFile,TH1D,TH2D,gDirectory
import gfit
import pipath, os
import datetime
import sys
import math

class spe_calc():
    def __init__(self,inputRFN=None):
        if inputRFN is None:
            sys.exit( 'spe_calc.__init__ NO INPUT ROOT FILE SPECIFIED')
        inputRootFileName = inputRFN
        self.GFIT = gfit.gfit()
        self.pip  = pipath.pipath()
        self.gU   = graphUtils.graphUtils()

        fmt = '%Y%m%d_%H%M%S_%f'
        cnow = datetime.datetime.now().strftime(fmt)
        d = outdir = os.path.join( 'FitResults', cnow )

        if os.path.isdir(d):
            pass
        else:
            try:
                os.mkdir(d)
            except IOError,e:
                print 'spe_calc__init__',e
            else:
                print 'spe_calc__init__ created',d
        self.figdir = outdir
        
        self.rfn = self.pip.fix( inputRootFileName )
        self.outrfn = os.path.join( outdir, inputRootFileName.replace(os.sep,'_').replace('Second','SPE_FIT') )
        self.rf = None
        self.Hists = {}
        self.Graphs = {}
        self.tableOfResults = {}

        self.evtsPerRunHistName = 'LED_events_per_run'
        self.evtsPerRunHist = None
        
        print 'spe_calc.__init__ ROOT files: input',self.rfn,'output',self.outrfn

        return
    def histLoop(self,drawEachFit=False,singleRunMode=False,selectRun=None):
        '''
        loop over histograms of spe distribution vs run, fit and plot mean, sigma, etc. vs run
        if drawEachFit, then separate output file showing each fit is draw
        if singleRunMode, do not create 'vs run' plots
        '''
        DebugRunLoop = False
        
        parNames = ['mean','emean','sigma','prob','mupois']
        pre = 'LED_WFD_At_vs_run_'
        postNames = ['S'+str(x) for x in range(6)]
        singleRun = -1
        
        #postNames = ['S2','S3'] ############### TEMPORARY
        
        for icol,pN in enumerate(postNames):
            h2n = pre + pN
            print 'spe_calc.histLoop: Fit',h2n
            self.tableOfResults[pN] = results = self.runLoop(h2n,FitFun='NGaus',draw=drawEachFit,projName=pN+'_py',debug=DebugRunLoop,selectRun=selectRun)
            X = sRuns = sorted( results.keys() )
            eX= [0. for x in X]
            print 'spe_calc.histLoop: Plot',h2n,'for',len(X),'runs'
            if singleRunMode: singleRun = int(X[0])

            if not singleRunMode:

                for i,v in enumerate( parNames ):
                    j = i + 1
                    Y = []
                    eY= []
                    for run in sRuns:
                        Y.append( results[run][j] )
                        if v=='mean':
                            eY.append( results[run][j+1] )
                        else:
                            eY.append( 0. )

                    name = 'SPE_'+v+'_vs_run_'+pN
                    title = name.replace('_',' ')
                    g = self.Graphs[name] = self.gU.makeTGraph(X,Y,title,name,ex=eX,ey=eY)
                    self.gU.color(g,icol,icol,setMarkerColor=True)

        if not singleRunMode:
            for parName in parNames:
                name = parName + '_vs_run'
                tmg = self.gU.makeTMultiGraph(name)
                ss = '_' + parName + '_'

                for gname in sorted( self.Graphs.keys() ) :
                    if ss in gname: tmg.Add( self.Graphs[gname] )
                self.Graphs[name] = tmg
                self.gU.drawMultiGraph(tmg, abscissaIsTime=False, drawLines=False, figdir=self.figdir )
                self.gU.drawMultiGraph(tmg, abscissaIsTime=False, drawLines=False, figdir=self.figdir ,SetLogy=True)
        else:
            # PRINT RESULTS IN SINGLE RUN MODE
            print '\nSPE fit results for run',singleRun
            Q = ['Counter']
            for x in parNames: Q.append(x)
            print ' '.join('%7s' % x for x in Q),' %7s' % 'GoodFit'
            for pN in sorted( self.tableOfResults.keys() ):
                results = self.tableOfResults[pN][singleRun]
                GoodFit = results[0]
                Q = results[1:]

                print '{0:7s}'.format(pN),
                print ' '.join('%7.3f' % x for x in Q),' %7s' % str(GoodFit)
            
                
        return
        
    def runLoop(self,h2n,FitFun='NGaus',draw=False,debug=False,projName=None,selectRun=None):
        '''
        loop over xbins in 2d hist, making projections and fitting them
        '''

        if selectRun is not None: print 'spe_calc.runLoop selectRun',selectRun
        
        h2 = self.Hists[h2n]
        hepr = self.evtsPerRunHist
        
        print 'spe_calc.runLoop h2',h2
        if h2 is None:
            h2 = self.getHist(h2n)
            print 'spe_calc.runLoop try again h2',h2
        pyName = "py"
        if projName is not None: pyName = projName
            
        nx = h2.GetXaxis().GetNbins()
        nxe= hepr.GetXaxis().GetNbins()
        if nx!=nxe:
            sys.exit('spe_calc.runLoop ERROR bin # mismatch. ' + h2.GetName() + ' ' + str(nx) + ' bins and ' + hepr.GetName() + ' ' + str(nxe) + ' bins')
            
        results = {}
        for ix in range(nx):
            jx = ix+1 # 1st bin = bin 1, last bin = bin nx. project from only a single bin
            if jx<=nx:
                run = int(h2.GetXaxis().GetBinCenter(jx))
                #print 'found run',run
                if selectRun is None or run==selectRun:
                    print 'spe_calc.runLoop processing run',run,'ix',ix,'jx',jx
                    py = h2.ProjectionY(pyName,jx,jx)
                    LEDevts = int(hepr.GetBinContent(jx))
                    Projevts= py.GetEntries()
                    guessmupois = .1
                    print 'LEDevts',LEDevts,'Projevts',Projevts,
                    if LEDevts>0: guessmupois = -math.log(1. - float(Projevts)/float(LEDevts) )
                    print 'guessmupois',guessmupois
                    inputPar = [ None, guessmupois, 35., 15.]  # C, poisMu, gausMu, gausSG
                    inputLim = [ [None,None], [1.e-4, min(10.,guessmupois*5)], [20., 50.], [1.e-3, 25.]]
                    if Projevts>5: 
                        if FitFun=='NGaus':
                            GoodFit,mean,emean, sgm, prob, mupois = self.GFIT.fitNGaus(py,debug=debug, inputPar=inputPar, inputLimits=inputLim)
                        else:
                            mupois = 1.
                            GoodFit,mean,emean, sgm, prob = self.GFIT.fit(py)
                        if debug: print 'ix,run#,GoodFit,mean,emean,sgm,prob {0} {6} {1} {2:.2f} {3:.2f} {4:.2f} {5:.3f}'.format(ix,GoodFit,mean,emean, sgm, prob, run)
                        if draw: self.gU.drawFit(py,figdir=self.figdir,extraName=str(run))
                        results[run] = [GoodFit,mean,emean,sgm,prob,mupois]
                    else:
                        results[run] = [False, -1., -1., -1., -1., -1.]
                        
        return results
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
                print 'this object',obj,'does not inherit from TH1'

        print 'spe_calc.getHists name, object for all TH1 found'
        for name in self.Hists:
            obj = self.Hists[name]
            if debug: print name,obj
        return  AllTH1
    def getHist(self,hn):
        hist = self.rf.Get(hn)
        return hist
    def main(self,drawEachFit=False,singleRunMode=False,selectRun=None):
        '''
        main module
        '''
        rfn = self.rfn
        OK = self.getHists(rfn)

        self.evtsPerRunHist = self.getHist( self.evtsPerRunHistName) 

        print 'spe_calc.main drawEachFit = ',drawEachFit
        self.histLoop(drawEachFit=drawEachFit,singleRunMode=singleRunMode,selectRun=selectRun)

        outrf = TFile(self.outrfn,'RECREATE')
        for g in self.Graphs: outrf.WriteTObject( self.Graphs[g] )
        outrf.Close()
        print 'Wrote',len(self.Graphs),'objects to',self.outrfn
        return
if __name__ == '__main__' :
    #h2n = 'LED_WFD_At_vs_run_S0'
        
    inputRootFileName = 'Second/20160224_170053_982090/second.root' 
    inputRootFileName = 'Second/20160304_125953_037208/second.root'
    inputRootFileName = 'Second/20160304_214633_049558/second.root' # full processing, but >1 hit pb
    inputRootFileName = 'Second/20160308_093800_546984/second.root' # >1 hit pb fixed, runs1010-1019

    inputRootFileName = 'Second/20160308_130649_412369/second.root' # run1010, after merge 20160308
    
    SC = spe_calc(inputRFN=inputRootFileName)

    selectRun = 1010
    singleRunMode = selectRun is not None
    drawEachFit = singleRunMode
    SC.main(drawEachFit=drawEachFit,selectRun=selectRun,singleRunMode=singleRunMode)
    
