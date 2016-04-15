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
        self.projHists = {}

        self.evtsPerRunHistName = 'LED_events_per_run'
        self.evtsPerRunHist = None
        
        print 'spe_calc.__init__ ROOT files: input',self.rfn,'output',self.outrfn

        return
    def histLoop(self,drawEachFit=False,singleRunMode=False,selectRun=None,saveProj=False):
        '''
        loop over histograms of spe distribution vs run, fit and plot mean, sigma, etc. vs run
        if drawEachFit, then separate output file showing each fit is draw
        if singleRunMode, do not create 'vs run' plots
        '''
        DebugRunLoop = False
        '''
        order of variables returned by runLoop:
          [GoodFit, mean,emean, sgm,esgm, mupois,emupois, prob]
        '''
        parNames = ['mean','emean','sigma','esigma','mupois','emupois','prob']
        XYpairs  = [ ['prob','mean'], ['prob','sigma'], ['prob','mupois'], ['mupois','mean'] ]
        pre = 'LED_WFD_At_vs_run_'
        postNames = ['S'+str(x) for x in range(6)]
        singleRun = -1
        
        
        for icol,pN in enumerate(postNames):
            h2n = pre + pN
            print 'spe_calc.histLoop: Fit',h2n
            self.tableOfResults[pN] = results = self.runLoop(h2n,FitFun='NGaus',draw=drawEachFit,projName=pN+'_py',debug=DebugRunLoop,selectRun=selectRun,saveProj=saveProj)
            
            sRuns = sorted( results.keys() )
            print 'spe_calc.histLoop: Plot',h2n,'for',len(sRuns),'runs for GoodFits'
            if singleRunMode: singleRun = int(sRuns[0])

            if not singleRunMode:

                for i,v in enumerate( parNames ):  # parameters vs run

                    j = i + 1
                    Y,eY, X,eX = [],[],[],[]
                    for run in sRuns:
                        GoodFit = results[run][0]
                        if GoodFit:
                            X.append( float(run) )
                            eX.append( 0. )
                            Y.append( results[run][j] )
                            if v in ['mean','sigma','mupois']:
                                eY.append( results[run][j+1] )
                            else:
                                eY.append( 0. )

                    name = 'SPE_'+v+'_vs_run_'+pN
                    title = name.replace('_',' ')
                    g = self.Graphs[name] = self.gU.makeTGraph(X,Y,title,name,ex=eX,ey=eY)
                    self.gU.color(g,icol,icol,setMarkerColor=True)
                    
                for pair in XYpairs:    # parameter vs parameter
                    xpair,ypair = pair
                    ix = parNames.index(xpair)+1
                    iex,iey = -1,-1
                    if 'e'+xpair in parNames: iex = parNames.index('e'+xpair)+1
                    iy = parNames.index(ypair)+1
                    if 'e'+ypair in parNames: iey = parNames.index('e'+ypair)+1
                    X,eX,Y,eY = [],[],[],[]
                    for run in sRuns:
                        GoodFit = results[run][0]
                        if GoodFit:
                            X.append( results[run][ix] )
                            Y.append( results[run][iy] )
                            if iex>0: eX.append( results[run][iex] )
                            else:     eX.append( 0. )
                            if iey>0: eY.append( results[run][iey] )
                            else:     eY.append( 0. )
                    name = 'SPE_'+ypair+'_vs_'+xpair+'_'+pN
                    title = name.replace('_',' ')
                    g = self.Graphs[name] = self.gU.makeTGraph(X,Y,title,name,ex=eX,ey=eY)
                    self.gU.color(g,icol,icol,setMarkerColor=True)

        if not singleRunMode:
            for parName in parNames:
                name = parName + '_vs_run'
                tmg = self.gU.makeTMultiGraph(name)
                ss = '_vs_' + parName + '_'
                ss = '_' +  parName + '_vs_run_'

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
        
    def runLoop(self,h2n,FitFun='NGaus',draw=False,debug=False,projName=None,selectRun=None,saveProj=False):
        '''
        loop over xbins in 2d hist, making projections and fitting them
        if saveProj, then give projections unique names and save them in output root file
        '''

        if selectRun is not None: print 'spe_calc.runLoop selectRun',selectRun
        
        h2 = self.Hists[h2n]
        hepr = self.evtsPerRunHist
        
        if debug: print 'spe_calc.runLoop h2',h2
        if h2 is None:
            h2 = self.getHist(h2n)
            if debug: print 'spe_calc.runLoop try again h2',h2
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

                if selectRun is None or run==selectRun:
                    if debug: print 'spe_calc.runLoop processing run',run,'ix',ix,'jx',jx
                    if saveProj: pyName = projName + '_'+str(run)
                    py = h2.ProjectionY(pyName,jx,jx)
                    if saveProj: self.projHists[pyName] = py
                    LEDevts = int(hepr.GetBinContent(jx))
                    Projevts= py.GetEntries()
                    if debug: print 'spe_calc.runLoop LEDevts',LEDevts,'Projevts',Projevts,
                    GoodFit = False
                    if LEDevts>100 and Projevts>5 and Projevts<LEDevts: 
                        guessmupois = .1
                        if LEDevts>0: guessmupois = -math.log(1. - float(Projevts)/float(LEDevts) )
                        if debug: print 'guessmupois',guessmupois

                        inputPar = [ None, guessmupois, None, 15.]  # C, poisMu, gausMu, gausSG
                        if FitFun=='NGaus':
                            GoodFit,mean,emean, sgm,esgm, mupois,emupois, prob = self.GFIT.fitNGaus(py,debug=debug, inputPar=inputPar) 
                        else:
                            mupois = 1.
                            esgm = emupois = 0.
                            GoodFit,mean,emean, sgm, prob = self.GFIT.fit(py)
                            
                        if debug: print 'ix,run#,GoodFit,mean,emean,sgm,esgm,mupois,emupois,prob {0} {1} {2:.3f} {4} {5:.2f} {6:.2f} {7:.2f} {8:.2f} {9:.2f}'.format(ix,run,prob,GoodFit,mean,emean, sgm,esgm, mupois,emupois)
                            
                        if draw: self.gU.drawFit(py,figdir=self.figdir,extraName=str(run))
                        results[run] = [GoodFit, mean,emean, sgm,esgm, mupois,emupois, prob]
                    else:
                        results[run] = [GoodFit,   -1., -1.,   -1., -1., -1.,-1.,        -1.]
                    print 'spe_calc.runLoop',h2.GetName(),'run',run,'evts',Projevts,'trigs',LEDevts,'GoodFit',GoodFit
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
    def main(self,drawEachFit=False,singleRunMode=False,selectRun=None,saveProj=False):
        '''
        main module
        '''
        rfn = self.rfn
        OK = self.getHists(rfn)

        self.evtsPerRunHist = self.getHist( self.evtsPerRunHistName) 

        print 'spe_calc.main drawEachFit = ',drawEachFit
        self.histLoop(drawEachFit=drawEachFit,singleRunMode=singleRunMode,selectRun=selectRun,saveProj=saveProj)

        outrf = TFile(self.outrfn,'RECREATE')
        for g in self.Graphs: outrf.WriteTObject( self.Graphs[g] )
        nh = 0
        if saveProj:
            nh = len(self.projHists)
            for h in self.projHists: outrf.WriteTObject( self.projHists[h] )
        outrf.Close()
        print 'Wrote',len(self.Graphs)+nh,'objects to',self.outrfn
        return
if __name__ == '__main__' :
    #h2n = 'LED_WFD_At_vs_run_S0'
        
    inputRootFileName = 'Second/20160224_170053_982090/second.root' 
    inputRootFileName = 'Second/20160304_125953_037208/second.root'
    inputRootFileName = 'Second/20160304_214633_049558/second.root' # full processing, but >1 hit pb
    inputRootFileName = 'Second/20160308_093800_546984/second.root' # >1 hit pb fixed, runs1010-1019

    #inputRootFileName = 'Second/20160308_130649_412369/second.root' # run1010, after merge 20160308
    inputRootFileName = 'Second/20160309_110857_763170/second.root' # full processing runs 585-1346, 1hit pb fixed

    inputRootFileName = 'Second/20160310_165505_691918/second.root' # first run# 1385 last run# 1408
    inputRootFileName = 'Second/20160311_184040_239400/second.root' # first runs 1409-1428
    
    SC = spe_calc(inputRFN=inputRootFileName)

    selectRun = 1161
    selectRun = None
    singleRunMode = selectRun is not None
    drawEachFit = singleRunMode
    saveProj = True
    SC.main(drawEachFit=drawEachFit,selectRun=selectRun,singleRunMode=singleRunMode,saveProj=saveProj)
    
