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
        
        
        return
    def startRun(self,file=None):
        self.R.start(fn=file)
        self.R.reportRunDetail()
        #self.R.testGet()
        self.bookHists('TDC')
        return
    def finish(self,rfn='oneton.root'):
        rf = TFile(rfn,"RECREATE")
        nwrite = 0
        for h in self.Hists: rf.WriteTObject(self.Hists[h])
        nwrite += len(self.Hists)

        rf.Close()
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
                dx = self.R.getTDCConfig('width')
                nx = int(math.pow(2,12))
                xmi = 0.
                xma = xmi + float(nx)*dx
                Hists[name] = TH1D(name,title,nx,xmi,xma)
                TDChistnames['LED'+w] = name
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
        dumpOn = False
        
        CalData = self.R.getCalData()
        EvtDir  = self.R.getEvtDir()
        TDCmd   = self.R.getModuleData('TDC','ChanDesc')
        Hists,TDChistnames = self.Hists,self.TDChistnames

        
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
                #print ' ',x,EvtDir[evtkey],'iterable=',hasattr(X,'__iter__'),'shape',X.shape,'dtype',X.dtype,'size',X.size
                if dumpEvt:
                    print x,
                    if self.R.printDataset(X) :
                        for y in X:print y,
                        print ''
                if 'TDC' in x:
                    TDC = self.R.unpackTDC(X)
                    #print 'TDC',TDC
                    TDCns = self.R.unpackTDC(X,mode='inns')
                    #print 'TDCns',TDCns
                    TDChits = []
                    TDCon   = []
                    for i,md in enumerate(TDCmd):
                        for pair,pairns in zip(TDC,TDCns):
                            if int(pair[0])==i:  # this TDC had a hit
                                TDChits.append(pair[0])
                                TDCon.append(md)
                                Hists[TDChistnames['raw'+md]].Fill(pair[1])
                                Hists[TDChistnames['ns'+md]].Fill(pairns[1])
                    if 'M' not in TDCon and 'CT' not in TDCon and 'CT+M+LED' in TDCon:
                        # arcane method to determine only LED trigger fired
                        for pairns in TDCns:
                            md = TDCmd[int(pairns[0])]
                            Hists[TDChistnames['LED'+md]].Fill(pairns[1])

                    if len(TDChits)>0:
                        for x in TDChits:
                            for y in TDChits:
                                if y>x: Hists['TDC_vs_TDC'].Fill(x,y)
        return

if __name__ == '__main__' :
    nevt = 9999999999
    if len(sys.argv)>1: nevt = int(sys.argv[1])
    P = process()
    fnlist = ['/Users/djaffe/work/1TonData/run462.h5','/Users/djaffe/work/1TonData/run463.h5','/Users/djaffe/work/1TonData/run464.h5']
    for fn in fnlist:
        rfn = fn.split('/')[-1]
        rfn = 'Output/'+rfn.split('.')[0]+'.root'
        P.startRun(file=fn)
        P.eventLoop(nevt)
        P.finish(rfn)
