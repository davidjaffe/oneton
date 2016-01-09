#!/usr/bin/env python
'''
process 1ton data
20151231
'''
import h5py
import numpy
import graphUtils
from ROOT import TFile,TH1D,TH2D,gDirectory
import sys
import datetime
import reader
import math

class process():
    def __init__(self):
        self.R = reader.reader()
        self.gU= graphUtils.graphUtils()
        self.Hists = {}
        self.TDChistnames = {}
        self.refTDCs = ['S0','S2']
        
        
        return
    def startRun(self,file=None):
        self.R.start(fn=file)
        self.R.reportRunDetail()
        self.R.testGet()
        self.bookHists('TDC')
        return
    def finish(self,rfn='oneton.root'):
        rf = TFile(rfn,"RECREATE")
        nwrite = 0
        for h in self.Hists: rf.WriteTObject(self.Hists[h])
        nwrite += len(self.Hists)

        rf.Close()

        # draw selected hists
        runnum = self.R.getRunNum()
        srun = str(runnum)
        for ref in self.refTDCs:
            s = []
            for h in self.Hists:
                if 'dTDC' in h and h[-2:]==ref:
                    s.append(self.Hists[h])

            if len(s)>0:
                s.sort()
                self.gU.drawMultiHists(s,fname=srun+'_dTwrt'+ref,figdir='Figures/TDC/',setLogy=True)

        
        # next line deletes all hists/trees from memory according to Rene Brun 3Apr2002
        gDirectory.GetList().Delete()
        print 'process.finish Wrote',nwrite,'objects to',rfn,'and deleted them from memory'
        return
    def bookHists(self,kind):
        '''
        book hists based on kind of data
        '''
        Hists,TDChistnames = self.Hists,self.TDChistnames
        OK = False
        if kind.lower()=='tdc':
            OK = True
            md = self.R.getModuleData('TDC','ChanDesc')
            # 1d : raw and in ns
            for w in md:
                title = 'raw TDC ' + w.replace('+','or')
                name = title.replace(' ','_')
                nx = int(math.pow(2,12))
                xmi = 0.
                xma = xmi + float(nx)
                Hists[name] = TH1D(name,title,nx,xmi,xma)
                TDChistnames['raw'+w] = name

                title = 'ns TDC ' + w.replace('+','or')
                name = title.replace(' ','_')
                dx = self.R.getTDCConfig('width')
                nx = int(math.pow(2,12))
                xmi = 0.
                xma = xmi + float(nx)*dx
                Hists[name] = TH1D(name,title,nx,xmi,xma)
                TDChistnames['ns'+w] = name

                title = 'LED TDC ' + w.replace('+','or')
                name = title.replace(' ','_')
                Hists[name] = TH1D(name,title,nx,xmi,xma)
                TDChistnames['LED'+w] = name

                # time difference hists
                for ref in self.refTDCs:
                    if w[0]=='S' and w!=ref:
                        title = 'LED dTDC ' + w + '-' + ref
                        name  = title.replace(' ','_').replace('-','m')
                        nx = 500
                        xmi,xma = -20.,20.
                        Hists[name] = TH1D(name,title,nx,xmi,xma)
                        TDChistnames['dTDC'+w+'-'+ref] = name
                
            # 2d : hits vs hits
            nx = len(md)
            xmi = -0.5
            xma = xmi + float(nx)
            title = name = 'TDC_vs_TDC'
            Hists[name] = h = TH2D(name,title,nx,xmi,xma,nx,xmi,xma)
            for i,w in enumerate(md):
                h.GetXaxis().SetBinLabel(i+1,w)
                h.GetYaxis().SetBinLabel(i+1,w)
        if OK: print 'process.bookHists booked',kind,'hists'
        return
                
    def eventLoop(self,maxEvt=99999999):
        '''
        loop over events in order of increasing event number
        '''
        dumpOn = True
        
        CalData = self.R.getCalData()
        EvtDir  = self.R.getEvtDir()

        
        sEvtKeys = sorted(EvtDir.keys(),key=lambda b: int(b))
        for evtkey in sEvtKeys: # sorted(EvtDir.keys()):
            evtnum = int(evtkey)
            printEvt = int(evtnum)%1000==0
            dumpEvt = dumpOn and int(evtnum)%10000==1
            if evtnum>maxEvt: break
            if printEvt: 
                print 'Event#',evtkey,'_________'
            for x in EvtDir[evtkey]:
                X = EvtDir[evtkey][x]

                if dumpEvt:
                    print x,
                    if self.R.printDataset(X) :
                        for y in X:print y,
                        print ''
                    if 'TDC' in x: print 'triggers',self.R.unpackTrigger(X)
                if 'TDC' in x: self.analTDC(X)
                if 'Digitizer' in x: self.analWFD(X)
        return
    def analWFD(self,raw):
        return
    def analTDC(self,raw):
        '''
        analyze TDC data for a single event
        '''
        TDCmd   = self.R.getModuleData('TDC','ChanDesc')
        Hists,TDChistnames = self.Hists,self.TDChistnames
        TDC = self.R.unpackTDC(raw)
        TDCns = self.R.unpackTDC(raw,mode='inns')
        trigs = self.R.unpackTrigger(raw)
        LEDonly = trigs==['LED']
        TDChits = []
        TDCvalues=[]
        TDCon   = []
        for i,md in enumerate(TDCmd):
            for pair,pairns in zip(TDC,TDCns):
                if int(pair[0])==i:  # this TDC had a hit
                    TDChits.append(pair[0])    # channel #
                    TDCvalues.append(pairns[1])  # value in ns
                    TDCon.append(md)           # name
                    Hists[TDChistnames['raw'+md]].Fill(pair[1])
                    Hists[TDChistnames['ns'+md]].Fill(pairns[1])
        if LEDonly:
            for md,t in zip(TDCon,TDCvalues):
                Hists[TDChistnames['LED'+md]].Fill(t)
            for ref in self.refTDCs:
                if ref in TDCon:
                    tref = TDCvalues[TDCon.index(ref)]
                    #print ref,tref
                    for md,t in zip(TDCon,TDCvalues):
                        if md[0]=='S' and md!=ref:
                            #print md,t,t-tref
                            Hists[TDChistnames['dTDC'+md+'-'+ref]].Fill(t-tref)

        # correlations between TDC channels
        if len(TDChits)>0:
            for x in TDChits:
                for y in TDChits:
                    if y>x: Hists['TDC_vs_TDC'].Fill(x,y)

        return

if __name__ == '__main__' :
    nevt = big = 9999999999
    if len(sys.argv)>1: nevt = int(sys.argv[1])
    P = process()
    fnlist = ['/Users/djaffe/work/1TonData/Filled_151217/run462.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run463.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run464.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run493.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run494.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run495.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run496.h5']
    if nevt<big: fnlist = fnlist[0] # debug with one file
    for fn in fnlist:
        rfn = fn.split('/')[-1]
        rfn = 'Output/'+rfn.split('.')[0]+'.root'
        P.startRun(file=fn)
        P.eventLoop(nevt)
        P.finish(rfn)
