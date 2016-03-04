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


class spe_calc():
    def __init__(self):
        self.GFIT = gfit.gfit()
        self.gU   = graphUtils.graphUtils()
        self.rfn = ' Second/20160224_170053_982090/second.root'
        self.rfn = 'Second/20160304_101658_998351/second.root'
        self.outrfn = self.rfn.replace('/','_').replace('Second','SPE_FIT')
        self.rf = None
        self.Hists = {}
        self.Graphs = {}

        return
    def histLoop(self,drawEachFit=False):
        '''
        loop over histograms of spe distribution vs run, fit and plot mean, sigma, etc. vs run
        '''

        parNames = ['mean','emean','sigma','prob','mupois']
        pre = 'LED_WFD_At_vs_run_'
        postNames = ['S'+str(x) for x in range(6)]

        
        #postNames = ['S2','S3'] ############### TEMPORARY
        
        for icol,pN in enumerate(postNames):
            h2n = pre + pN
            print 'spe_calc.histLoop: Fit',h2n
            results = self.runLoop(h2n,FitFun='NGaus',draw=drawEachFit,projName=pN+'_py')
            print 'spe_calc.histLoop: Plot',h2n
            X = sRuns = sorted( results.keys() )
            eX= [0. for x in X]

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

        for parName in parNames:
            name = parName + '_vs_run'
            tmg = self.gU.makeTMultiGraph(name)
            ss = '_' + parName + '_'

            for gname in sorted( self.Graphs.keys() ) :
                if ss in gname: tmg.Add( self.Graphs[gname] )
            self.Graphs[name] = tmg
            self.gU.drawMultiGraph(tmg, abscissaIsTime=False, drawLines=False, figdir='FitResults/')
            self.gU.drawMultiGraph(tmg, abscissaIsTime=False, drawLines=False, figdir='FitResults/',SetLogy=True)
            
                
        return
        
    def runLoop(self,h2n,FitFun='NGaus',draw=False,debug=False,projName=None):
        '''
        loop over xbins in 2d hist, making projections and fitting them
        '''
        h2 = self.Hists[h2n]
        print 'spe_calc.runLoop h2',h2
        if h2 is None:
            h2 = self.getHist(h2n)
            print 'spe_calc.runLoop try again h2',h2
        pyName = "py"
        if projName is not None: pyName = projName
        nx = h2.GetXaxis().GetNbins()
        results = {}
        for ix in range(nx):
            jx = ix+1
            if jx<=nx:
                run = int(h2.GetXaxis().GetBinCenter(ix))
                py = h2.ProjectionY(pyName,ix,jx)
                if py.GetEntries()>5: 
                    if FitFun=='NGaus':
                        GoodFit,mean,emean, sgm, prob, mupois = self.GFIT.fitNGaus(py,debug=debug)
                    else:
                        mupois = 1.
                        GoodFit,mean,emean, sgm, prob = self.GFIT.fit(py)
                    if debug: print 'ix,run#,GoodFit,mean,emean,sgm,prob {0} {6} {1} {2:.2f} {3:.2f} {4:.2f} {5:.3f}'.format(ix,GoodFit,mean,emean, sgm, prob, run)
                    if draw: self.gU.drawFit(py,figdir='FitResults/',extraName=str(run))
                    results[run] = [GoodFit,mean,emean,sgm,prob,mupois]
                        
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
if __name__ == '__main__' :
    #h2n = 'LED_WFD_At_vs_run_S0'
    SC = spe_calc()
    rfn = SC.rfn
    OK = SC.getHists(rfn)

    #SC.runLoop(h2n,draw=True)
    SC.histLoop(drawEachFit=False)

    outrf = TFile(SC.outrfn,'RECREATE')
    for g in SC.Graphs: outrf.WriteTObject( SC.Graphs[g] )
    outrf.Close()
    print 'Wrote',len(SC.Graphs),'objects to',SC.outrfn
