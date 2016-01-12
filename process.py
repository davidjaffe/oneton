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
        self.WFHists = []
        
        
        return
    def startRun(self,file=None):
        print 'process.startRun',file
        self.R.start(fn=file)
        self.R.reportRunDetail()
        self.R.testGet()
        self.bookHists('TDC')
        return
    def endRun(self):
        self.R.closeHDF5File()
        return
    def finish(self,rfn='oneton.root'):
        self.R.summary()
        rf = TFile(rfn,"RECREATE")
        nwrite = 0
        for h in self.Hists: rf.WriteTObject(self.Hists[h])
        nwrite += len(self.Hists)
        for h in self.WFHists: rf.WriteTObject(h)
        nwrite += len(self.WFHists)
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
            Event = EvtDir[evtkey]
            evtnum = int(evtkey)
            printEvt = int(evtnum)%1000==0
            dumpEvt = dumpOn and (int(evtnum)%10000==1 or numpy.random.choice([True,False]))
            if evtnum>maxEvt: break
            if printEvt: 
                print 'Event#',evtkey,'_________'
            missing = self.R.unpackEvent(Event)
            if len(missing)==0: # OK
                triggers = self.R.unpackTrigger(Event['TDC'])
                self.analTDC(Event['TDC'])
                Draw = triggers==['CT']
                X1 = Event['Digitizer_1']
                X2 = Event['Digitizer_2']
                self.analWFD(X1,'Digitizer_1',evtnum)
                self.analWFD(X2,'Digitizer_2',evtnum)
                
            if dumpEvt:
                self.dumpEvent(Event,evtnum)
                self.evtDisplay(Event,evtnum)
        return
    def dumpEvent(self,raw,evtnum):
        print 'begin dump event#',evtnum,'--------------'
        for x in raw:
            X = raw[x]
            print x,
            if self.R.printDataset(X) :
                for y in X:print y,
                print ''

        print '---unpacked data with channel description---'
        for Q in self.R.TDCNames:
            U = self.R.unpackTDC(raw[Q])
            CD = self.R.addChanDesc(U,Q)
            print Q+"hits",
            for x,y in zip(CD,U): print x,y,
            print ''
        for Q in self.R.QDCNames:
            U = self.R.unpackQDC(raw[Q])
            CD = self.R.addChanDesc(U,Q)
            print Q+":",
            for x,y in zip(CD,U): print x,y,
            print ''
        for Q in self.R.WFDNames:
            U = self.R.unpackDigitizer(raw[Q])
            CD = self.R.addChanDesc(U,Q)
            print Q+":",
            for x,y in zip(CD,U): print x,U[y],
            print ''
        print '---------------- end dump event#',evtnum
        return
    def evtDisplay(self,raw,evtnum):
        '''
        concoct event display
        If both high and low gain results are available for QDC, use high gain value.
        '''

        # associate channel name with channel number and readout
        U = {}
        QDCcd = {}
        for Q in self.R.QDCNames:
            U[Q] = self.R.unpackQDC(raw[Q])
            cd  = self.R.addChanDesc(U[Q],Q)
            for x,y in zip(cd,U[Q]):
                if x in QDCcd:
                    if QDCcd[x][2]==0: QDCcd[x] = y
                else:
                    QDCcd[x] = y
                
        TDCcd,U = {},{}
        for Q in self.R.TDCNames:
            U[Q] = self.R.unpackTDC(raw[Q])
            cd  = self.R.addChanDesc(U[Q],Q)
            for x,y in zip(cd,U[Q]): TDCcd[x] = y

        U,WFDcd = {},{}
        channels = []
        for Q in self.R.WFDNames:
            if self.R.CalData is None:
                U[Q] = self.R.unpackDigitizer(raw[Q])
            else:
                U[Q] = self.R.unpackDigitizer(raw[Q],self.R.CalData[Q]['Pedestal'])
            cd  = self.R.addChanDesc(U[Q],Q)
            channels.extend(cd)
            for x,y in zip(cd,U[Q]): WFDcd[x] = y

        # create hist titles. first hist has triggers
        triggers = self.R.unpackTrigger(raw['TDC'])
        c = 'Trigs'
        for t in triggers: c += ' ' + t
                
        ipt = []
        for cd in channels:
            w = ''
            if channels.index(cd)==0: w = c
            w += ' ' + cd + ' '
            if cd in QDCcd:
                qdc = QDCcd[cd][1]
                w += 'QDC '+str(qdc)
            else:
                w += ' noQDC '
            if cd in TDCcd:
                tdc = TDCcd[cd][1]
                w += 'TDC '+str(tdc)
            else:
                w += ' noTDC '
            ipt.append(w)

        # create histograms. output file name contains run, event number
        i0 = 0
        H = []
        for Q in self.R.WFDNames:
            #print 'process.evtDisplay i0,ipt[i0:]',i0,ipt[i0:]
            h = self.R.displayWFD(U[Q],ipt=ipt[i0:])
            H.extend(h)
            i0 += len(h)


        for h in H:
            for a in ["x","y"]:
                s = h.GetLabelSize(a)
                h.SetLabelSize(2.*s,a)
        s = 'run'+str(self.R.getRunDetail('run'))+'evt'+str(evtnum)
        self.gU.drawMultiHists(H,fname=s,figdir='Figures/WF/',statOpt=0,dopt='Hist')
        
        for h in H: h.Delete()
                    
        return
    def AinB(self,A,B):
        for a in B:
            if a==A: return True
        return False
    def analWFD(self,raw,module,evtnum,Draw=False):
        '''
        unpack and pedestal subtract waveforms
        FIXME
        '''
        runnum = self.R.getRunDetail('run')

        WFD   = self.R.unpackDigitizer(raw)
        if self.R.CalData is not None:
            WFDps = self.R.unpackDigitizer(raw,self.R.CalData[module]['Pedestal'])

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
    fnlist = [                  '/Users/djaffe/work/1TonData/Filled_151217/run493.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run494.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run495.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run496.h5']

    fnlist = ['/Users/djaffe/work/1TonData/Filled_151217/run537.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run538.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run539.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run540.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run541.h5',
                '/Users/djaffe/work/1TonData/Filled_151217/run542.h5']
    if nevt<big: fnlist = [fnlist[0]] # debug with one file
    
    for fn in fnlist:
        rfn = fn.split('/')[-1]
        rfn = 'Output/'+rfn.split('.')[0]+'.root'
        P.startRun(file=fn)
        P.eventLoop(nevt)
        P.finish(rfn)
        P.endRun()
