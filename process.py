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
        self.refTDCs = ['S2','H0']
        self.WFHists = []

        self.QDCrange = [1.e20,-1.e20]
        
        
        return
    def start(self):
        self.bookHists('TDC')
        self.bookHists('QDC')
        return
    def startRun(self,file=None,tG=False):
        print 'process.startRun',file
        self.R.start(fn=file)
        self.R.reportRunDetail()
        if tG: self.R.testGet()
        return
    def endRun(self):
        self.R.closeHDF5File()
        return
    def finish(self,rfn='oneton.root'):
        '''
        execute at end of job
        '''
        self.R.summary()
        rf = TFile(rfn,"RECREATE")
        nwrite = 0
        for h in self.Hists: rf.WriteTObject(self.Hists[h])
        nwrite += len(self.Hists)
        for h in self.WFHists: rf.WriteTObject(h)
        nwrite += len(self.WFHists)
        rf.Close()

        print 'process.finish QDCrange',self.QDCrange

        # draw selected hists
        srun = self.getFilePrefix(rfn)
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
        TL = ['All']
        TL.extend(self.R.uniqueTriggers)
        OK = False
        if kind.lower()=='tdc':
            OK = True
            dx = self.R.getTDCConfig('width')
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

                # ns hists for all and each unique trigger
                for tl in TL:
                    nx = int(math.pow(2,12))
                    xmi = 0.
                    xma = xmi + float(nx)*dx
                    title = tl+' TDC ' + w.replace('+','or')
                    name = title.replace(' ','_')
                    Hists[name] = TH1D(name,title,nx,xmi,xma)
                    TDChistnames[tl+w] = name

                    # time difference hists
                    # differences between hodos and signal separately
                    for ref in self.refTDCs:
                        if w[0]==ref[0] and w!=ref:
                            title = tl+' dTDC ' + w + '-' + ref
                            name  = title.replace(' ','_').replace('-','m')
                            nx = 500
                            xmi,xma = -20.,20.
                            Hists[name] = TH1D(name,title,nx,xmi,xma)
                            TDChistnames[tl+'dTDC'+w+'-'+ref] = name
                
            # 2d : hits vs hits
            nx = len(md)
            xmi = -0.5
            xma = xmi + float(nx)
            title = name = 'TDC_vs_TDC'
            Hists[name] = h = TH2D(name,title,nx,xmi,xma,nx,xmi,xma)
            for i,w in enumerate(md):
                h.GetXaxis().SetBinLabel(i+1,w)
                h.GetYaxis().SetBinLabel(i+1,w)
            # 2d : channel vs time (show all in single plot)
            nx = int(math.pow(2,12))
            xmi = 0.
            xma = xmi + float(nx)*dx
            ny = len(md)
            ymi = -0.5
            yma = ymi + float(ny)
            for tl in TL:
                title = name = 'chan_vs_TDC_' + tl
                Hists[name] = h = TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
                for i,w in enumerate(md):
                    h.GetYaxis().SetBinLabel(i+1,w)

            #print 'process.bookHists',kind
            #for key in TDChistnames: print 'TDChistnames['+key+']',TDChistnames[key]
        if kind=='QDC':
            OK = True
            md = self.R.getModuleData('QDC','ChanDesc')
            clean = False
            while not clean:
                md.remove('N/C')
                clean = 'N/C' not in md
            for tl in TL:
                for w in md:
                    name = tl + '_QDC_'+w
                    title = name.replace('_',' ')
                    nx = 100
                    xmi,xma = 0.,4100.
                    #print 'process.bookHists name,title,nx,xmi,xma',name,title,nx,xmi,xma
                    Hists[name] = TH1D(name,title,nx,xmi,xma)
                    for r in ['lo','hi']: # lo=1,hi=0
                        name = tl+'_QDC'+r+'_'+w
                        title = name.replace('_',' ')
                        nx = 100
                        if r=='lo': xmi,xma = 0.,600.
                        if r=='hi': xmi,xma = 500.,4100.
                        Hists[name] = TH1D(name,title,nx,xmi,xma)
        if OK: print 'process.bookHists booked',kind,'hists'
        return
                
    def eventLoop(self,maxEvt=99999999):
        '''
        loop over events in order of increasing event number
        '''
        dumpOn,dumpThres = True,0.999
        
        CalData = self.R.getCalData()
        EvtDir  = self.R.getEvtDir()

        sEvtKeys = sorted(EvtDir.keys(),key=lambda b: int(b))
        for evtkey in sEvtKeys: # sorted(EvtDir.keys()):
            Event = EvtDir[evtkey]
            evtnum = int(evtkey)
            if evtnum>maxEvt: break

            if int(evtnum)%1000==0: print 'Event#',evtkey,'_________'

            dumpEvt = dumpOn and (int(evtnum)%10000==1 or numpy.random.uniform()>dumpThres)
            
            missing = self.R.unpackEvent(Event)
            if len(missing)==0: # OK
                triggers = self.R.unpackTrigger(Event['TDC'])
                self.analTDC(Event['TDC'])
                self.analQDC(Event,evtnum)
                Draw = triggers==['CT']
                if 0:
                    X1 = Event['Digitizer_1']
                    X2 = Event['Digitizer_2']
                    self.analWFD(X1,'Digitizer_1',evtnum)
                    self.analWFD(X2,'Digitizer_2',evtnum)
                
            if dumpEvt:
                self.dumpEvent(Event,evtnum)
                self.evtDisplay(Event,evtnum)
        return
    def dumpEvent(self,raw,evtnum):
        triggers = self.R.unpackTrigger(raw['TDC'])
        T = 'Triggers='
        for t in triggers: T += ' ' + t
        print '\n begin dump event#',evtnum,T,'--------------'
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
        QDCcd = self.R.assocQDC(raw)
        TDCcd = self.R.assocTDC(raw,mode='ns')
        WFDcd = self.R.assocWFD(raw,pedsub=True)
        channels = sorted(WFDcd.keys())


        # create hist titles. first hist has triggers
        triggers = self.R.unpackTrigger(raw['TDC'])
        c = 'Trigs'
        for t in triggers: c += ' ' + t

        H = []
        for cd in channels:
            w = ''
            if channels.index(cd)==0: w = c
            if cd in QDCcd:
                w += ' QDC '
                for y in QDCcd[cd]:  w += '{0:.2f} '.format(y[1])
            else:
                w += ' noQDC '
            if cd in TDCcd:
                w += ' TDC {0:.2f} '.format( TDCcd[cd][1] )
            else:
                w += ' noTDC '
            h = self.R.displayWFDnew(WFDcd[cd],cd,pretitle=w)
            H.append(h)

        # increase size of x,y labels
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
    def analQDC(self,raw,evtnum):
        '''
        unpack and analyze QDCs
        Histogram of lo,hi range
        Compare lo,hi range when they overlap
        '''
        runnum = self.R.getRunDetail('run')
        moniker = 'Run ' + str(runnum) + ' Event ' + str(evtnum)
        TL = ['All']
        TL.extend(self.R.unpackTrigger(raw['TDC']))

        
        # create a dictionary with key = channel, value = list of QDC values (can have 2 entries, 1 for lo- and 1 for hi-range)
        QDC = self.R.assocQDC(raw)
        cd = sorted(QDC.keys())
        rng = ['hi','lo'] # low,high range: hilo =1,0
        for x in cd:
            for y in QDC[x]:
                #print 'process.analQDC x,QDC[x],y:',x,QDC[x],y
                ch,v,hilo,ovfl = y # chan#, value, hirange or lowrange, overflow
                self.QDCrange[0]=min(self.QDCrange[0],v)
                self.QDCrange[1]=max(self.QDCrange[1],v)
                for tl in TL:
                    name = tl+'_QDC'+rng[hilo]+'_'+x
                    self.Hists[name].Fill(v)
                    name = tl+'_QDC_'+x
                    self.Hists[name].Fill(v)
        
        for x in QDC:
            
            if len(QDC[x])>1 and x!='N/C':
                print 'process.QDC',moniker,'overlap',x,QDC[x]
                    
        return
        
    def analWFD(self,raw,module,evtnum,Draw=False):
        '''
        unpack and pedestal subtract waveforms
        FIXME
        '''
        display = False
        
        runnum = self.R.getRunDetail('run')
        pretitle = 'r'+str(runnum)+'e'+str(evtnum)

        WFD = self.R.unpackDigitizer(raw)

        if self.R.CalData is not None:
            WFDps = self.R.unpackDigitizer(raw,self.R.CalData[module]['Pedestal'])
            if display:
                cd  = self.R.addChanDesc(raw,module)
                ipt = []
                for x in cd: ipt.append(pretitle+x)
                h = self.R.displayWFD(WFDps,ipt=ipt)
                self.WFHists.extend(h)



        return
    def analTDC(self,raw):
        '''
        analyze TDC data for a single event
        '''
        TDCmd   = self.R.getModuleData('TDC','ChanDesc')
        Hists,TDChistnames = self.Hists,self.TDChistnames
        TL = ['All']
        TL.extend(self.R.unpackTrigger(raw))

        # unpack TDCs, get description of each channel with data
        # then form dict with key=channel_description,value=TDC value
        TDC = self.R.unpackTDC(raw)
        TDCns = self.R.unpackTDC(raw,mode='inns')
        cd = self.R.addChanDesc(raw,'TDC')
        rawCD, nsCD = {},{}
        for x,y in zip(cd,TDC)   : rawCD[x]=y
        for x,y in zip(cd,TDCns) : nsCD[x] =y

        TDChits = []


        for x in cd:
            i = TDCmd.index(x)
            TDChits.append(i)
            Hists[TDChistnames['raw'+x]].Fill(rawCD[x][1])
            for tl in TL:
                Hists[TDChistnames[tl+x]].Fill(nsCD[x][1])
        
        for ref in self.refTDCs:
            if ref in cd:
                tref = nsCD[ref][1]
                for x in cd:
                    if x!=ref and ref[0]==x[0]:
                        dt = nsCD[x][1]-tref
                        for tl in TL:
                            Hists[TDChistnames[tl+'dTDC'+x+'-'+ref]].Fill(dt)

        for tl in TL:
            name = 'chan_vs_TDC_' + tl
            for x in cd:
                y = float(TDCmd.index(x))
                Hists[name].Fill(nsCD[x][1],y)
                            
        # correlations between TDC channels
        if len(TDChits)>0:
            for x in TDChits:
                for y in TDChits:
                    if y>x: Hists['TDC_vs_TDC'].Fill(x,y)

        return
    def getFilePrefix(self,fn):
        '''
        get file name prefix from full path name
        '''
        f = fn.split('/')[-1]
        return f.split('.')[0]

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
    fnlist = ['/Users/djaffe/work/1TonData/Filled_151217/run585.h5',
            '/Users/djaffe/work/1TonData/Filled_151217/run586.h5']
    if nevt<big: fnlist = [fnlist[0]] # debug with one file

    # create name of output root file
    f = sorted(fnlist)
    f1,f2 = P.getFilePrefix(f[0]),P.getFilePrefix(f[-1])
    rfn = 'Output/'
    g = [f1]
    if f1!=f2: g.append(f2)
    for q in g:
        rfn += q
        if len(g)>1 and '-' not in rfn: rfn += '-'
    rfn += '.root'
    print 'Output root file name is',rfn

    first = True
    for fn in fnlist:
        P.startRun(file=fn)
        if first: P.start()
        P.eventLoop(nevt)
        P.endRun()
        first = False
    P.finish(rfn)
