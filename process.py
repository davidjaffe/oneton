#!/usr/bin/env python
'''
process 1ton data
20151231
'''
import h5py
import numpy
import graphUtils
import ROOT
from ROOT import TFile,TH1D,TH2D,gDirectory
import sys
import datetime
import reader
import math
import os
#import Logger
import wfanal

class process():
    def __init__(self):
        self.R = reader.reader()
        self.W = wfanal.wfanal()
        self.gU= graphUtils.graphUtils()
        self.Hists = {}
        self.TDChistnames = {}
        self.refTDCs = ['S2','H0']
        self.WFHists = []



        now = datetime.datetime.now()
        self.start_time = now
        fmt = '%Y%m%d_%H%M%S'
        cnow = now.strftime(fmt)
        parentDir = 'Results/'+cnow+'/'
        self.logdir = parentDir + 'Log/'
        self.figdir = parentDir + 'Figures/'
        self.WFfigdir = self.figdir + 'WF/'
        self.TDCfigdir= self.figdir + 'TDC/'
        self.outputdir= parentDir + 'Output/'
        dirs = [parentDir, self.logdir, self.figdir, self.WFfigdir, self.TDCfigdir, self.outputdir]
        for d in dirs:
            if os.path.isdir(d):
               pass
            else:
                try:
                    print("Skipping mkdir")
                    #os.mkdir(d)
                except IOError,e:
                    print 'process__init__',e
                else:
                    print 'process__init__ created',d
        lfn = self.logdir + cnow + '.log'
        #sys.stdout = Logger.Logger(fn=lfn)
        print 'process__init__ Output directed to terminal and',lfn
        print 'process__init__ Job start time',self.start_time.strftime('%Y/%m/%d %H:%M:%S')
        
        return
    def start(self):
        self.bookHists('TDC')
        self.bookHists('QDC')
        self.bookHists('WFD')
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
        Print summary output
        Open and fill ROOT file
        Draw some hists
        Delete all hists/trees in memory
        Print job timing stats
        '''
        self.R.summary()
        rf = TFile(rfn,"RECREATE")
        nwrite = 0
        for h in self.Hists: rf.WriteTObject(self.Hists[h])
        nwrite += len(self.Hists)
        for h in self.WFHists: rf.WriteTObject(h)
        nwrite += len(self.WFHists)
        rf.Close()

        # draw selected hists
        srun = self.getFilePrefix(rfn)
        for ref in self.refTDCs:
            s = []
            for h in self.Hists:
                if 'dTDC' in h and h[-2:]==ref:
                    s.append(self.Hists[h])

            if len(s)>0:
                s.sort()
                self.gU.drawMultiHists(s,fname=srun+'_dTwrt'+ref,figdir=self.TDCfigdir,setLogy=True)
        
        # next line deletes all hists/trees from memory according to Rene Brun 3Apr2002
        gDirectory.GetList().Delete()
        print 'process.finish Wrote',nwrite,'objects to',rfn,'and deleted them from memory'

        # how long did this job take?
        endtime = datetime.datetime.now()
        print 'process.finish Job end time',endtime.strftime('%Y/%m/%d %H:%M:%S')
        elapsed = endtime-self.start_time
        s = elapsed.total_seconds()
        h = int(s/60./60.)
        m = int((s-60.*60.*h)/60.)
        sc= s - 60.*60.*h - 60.*m
        dt = str(h)+'h'+str(m)+'m'+str(sc)
        print 'process.finish Job duration(s)',s,'or',dt
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
                # 1-D hists
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
                # 2-D summary hists
                ny = len(md)
                ymi = -0.5
                yma = ymi + float(ny)
                nx,xmi,xma = 120,0.,600.
                name = 'chan_vs_QDC_'+tl
                title = name.replace('_',' ')
                Hists[name] = h = TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
                #print 'process.bookHists name,Hists[name]',name,Hists[name]
                for i,w in enumerate(md):
                    h.GetYaxis().SetBinLabel(i+1,w)
        if kind=='WFD':
            OK = True
            md = self.R.getModuleData('Digitizer','ChanDesc')
            for tl in TL:
                ny = len(md)
                ymi = -0.5
                yma = ymi + float(ny)
                for pren in ['time','area','ped','pedsd','npulse','nsubp']:
                    name = pren + '_vs_WFD_' + tl
                    nx,xmi = 4096,-0.5
                    xma = xmi + float(nx)
                    if pren=='area':
                        nx,xma = 1000,100.
                        xmi = xma - float(nx)
                    if pren=='ped':
                        nx,xmi = 200,1950.
                        xma = xmi + 0.25*float(nx)
                    if pren=='pedsd':
                        nx,xmi = 50,-0.5
                        xma = xmi + float(nx)/5.
                    if pren=='npulse' or pren=='nsubp':
                        nx,xmi = 11,-0.5
                        xma = xmi + float(nx)
                    title = name.replace('_',' ')
                    Hists[name] = h = TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
                    for i,w in enumerate(md):
                        h.GetYaxis().SetBinLabel(i+1,w)
        if OK:
            print 'process.bookHists booked',kind,'hists'
        else:
            print 'process.bookHists Invalid input',kind
        return
    def eventLoop(self,maxEvt=99999999,dumpAll=False,dumpThres=0.9999):
        '''
        loop over events in order of increasing event number
        '''
        dumpOn = dumpThres>0 
        
        CalData = self.R.getCalData()
        EvtDir  = self.R.getEvtDir()

        self.dumpCalib()
        
        sEvtKeys = sorted(EvtDir.keys(),key=lambda b: int(b))
        for evtkey in sEvtKeys: # sorted(EvtDir.keys()):
            Event = EvtDir[evtkey]
            evtnum = int(evtkey)
            if evtnum>maxEvt: break

            if int(evtnum)%1000==0: print 'Event#',evtkey,'_________'

            dumpEvt = dumpAll or \
                (dumpOn and (int(evtnum)%10000==1 or numpy.random.uniform()>dumpThres))
            
            missing = self.R.unpackEvent(Event)
            if len(missing)==0: # OK
                triggers = self.R.unpackTrigger(Event['TDC'])
                self.analTDC(Event['TDC'])
                self.analQDC(Event,evtnum)
                self.analWFD(Event,evtnum)
                
            if dumpEvt:
                self.dumpEvent(Event,evtnum)
                if not dumpAll: self.evtDisplay(Event,evtnum)
        return
    def dumpCalib(self):
        '''
        dump calibration data
        '''
        CalData = self.R.getCalData()
        if CalData is None:
            print '=====> NO CALIBRATION DATA <==========='
            return
        print '+++++++++++++ dump Calibration data +++++++++'            
        for x in CalData:
            print x, # Digitizer_1,Digitizer_2
            for y in CalData[x]: # Address, Pedestal
                print y,
                g = CalData[x][y][()]
                if y=='Address':
                    print g
                elif y=='Pedestal':
                    for i in range(4):
                        print "Chan",i,g[i,:]
                else:
                    print 'Unexpected value ......'
                    
        print '++++++++++++++ end Calibration data dump ++++'
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

        NX = 480
        XMI= 19.5
        XMA= XMI+float(NX)
            
        H,HR = [],[]
        self.hWFD,self.hrWFD = {},{}
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
            self.hWFD[cd] = h = self.R.displayWFDnew(WFDcd[cd],cd,pretitle=w)
            H.append(h)
            self.hrWFD[cd] = hr= self.R.displayWFDnew(WFDcd[cd],cd,pretitle=w+' Z',NX=NX,XMI=XMI,XMA=XMA)
            HR.append(hr)

        # increase size of x,y labels
        for h,hr in zip(H,HR):
            for a in ["x","y"]:
                s = h.GetLabelSize(a)
                h.SetLabelSize(2.*s,a)
                s = hr.GetLabelSize(a)
                hr.SetLabelSize(2.*s,a)
        s = 'run'+str(self.R.getRunDetail('run'))+'evt'+str(evtnum)

        self.diagnoseWFD(fname=s,figdir=self.WFfigdir)

        self.gU.drawMultiHists(H,fname=s,figdir=self.WFfigdir,statOpt=0,dopt='Hist')
        for h in H: h.Delete()

        s += 'zoomed'
        self.gU.drawMultiHists(HR,fname=s,figdir=self.WFfigdir,statOpt=0,dopt='Hist')
        for hr in HR: hr.Delete()
        
                    
        return
    def diagnoseWFD(self,fname='diagnoseWFD',figdir=''):
        '''
        overlay histograms of waveforms with results of wfanal
        '''
        md = sorted(self.hWFD.keys())
        for i,DICT in enumerate([self.hWFD,self.hrWFD]):
            # arcane bullshit to draw lines...???
            lines = {}
            for icd,cd in enumerate(md):

                V = self.aWFD[cd] #### [ped,iPulse,subPperP,pArea,pTime]
                h = DICT[cd]
                ymi,yma = h.GetMinimum(),h.GetMaximum()
                xmi,xma = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
                h.SetStats(0)
                ped = V[0]
                L =  ROOT.TLine(xmi,ped,xma,ped)
                L.SetLineColor(ROOT.kGreen)
                lines[cd] = [ L ]
                pedsd = V[1]
                for x in [-1.,1.]:
                    y = ped+x*pedsd
                    L = ROOT.TLine(xmi,y,xma,y)
                    L.SetLineColor(ROOT.kGreen)
                    L.SetLineStyle(3)
                    lines[cd].append( L )
                iPulse = V[2]
                for pair in iPulse:
                    for ib in pair:
                        L = ROOT.TLine(float(ib),ymi,float(ib),yma)
                        L.SetLineColor(ROOT.kBlue)
                        L.SetLineStyle(3) # dotted
                        lines[cd].append( L )
                pTime = V[5]
                for t in pTime:
                    L = ROOT.TLine(t,ymi,t,yma)
                    L.SetLineColor(ROOT.kRed)
                    L.SetLineStyle(2) # dashed
                    lines[cd].append( L )

            pdf = fname + '_'+str(i) 
            if figdir!='': pdf = figdir + pdf
            ps = pdf + '.ps'
            pdf= pdf + '.pdf'
            ROOT.gROOT.ProcessLine("gROOT->SetBatch()") # no pop up
            xsize,ysize = 1100,850 # landscape style
            canvas = ROOT.TCanvas(pdf,fname,xsize,ysize)
            canvas.Divide(2,4)
            for icd,cd in enumerate(md):
                canvas.cd(icd+1)
                h = DICT[cd]
                h.SetStats(0)
                h.Draw("hist")
                for L in lines[cd]:
                    L.Draw()
                    canvas.Modified() # this one?
                    canvas.Update() # is this needed?
                    #print L.Print()
            if 0:
                canvas.cd()
                canvas.Update()
                canvas.Modified()    
                canvas.Print(ps,'Landscape')
                os.system('ps2pdf ' + ps + ' ' + pdf)
                if os.path.exists(pdf): os.system('rm ' + ps)
            self.gU.finishDraw(canvas,ps,pdf,ctitle=fname)
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

        md = self.R.getModuleData('QDC','ChanDesc')
        while 'N/C' in md: md.remove('N/C')
        
        # create a dictionary with key = channel, value = list of QDC values (can have 2 entries, 1 for lo- and 1 for hi-range)
        QDC = self.R.assocQDC(raw)
        cd = sorted(QDC.keys())
        rng = ['hi','lo'] # low,high range: hilo =1,0
        for x in cd:
            i = md.index(x)
            for y in QDC[x]:
                #print 'process.analQDC x,QDC[x],y:',x,QDC[x],y
                if len(y)==3: # early data
                    ch,v,hilo = y
                else:
                    ch,v,hilo,ovfl = y # chan#, value, hirange or lowrange, overflow
                for tl in TL:
                    name = tl+'_QDC'+rng[hilo]+'_'+x
                    self.Hists[name].Fill(v)
                    name = tl+'_QDC_'+x
                    self.Hists[name].Fill(v)
                    name = 'chan_vs_QDC_'+tl
                    self.Hists[name].Fill(v,float(i))
        
        for x in QDC:
            if len(QDC[x])>1 and x!='N/C':
                print 'process.QDC',moniker,'overlap',x,QDC[x]
                    
        return
        
    def analWFD(self,raw,evtnum):
        '''
        unpack and pedestal subtract waveforms
        FIXME
        '''
        display = False

        TL = ['All']
        TL.extend(self.R.unpackTrigger(raw['TDC']))
                
        runnum = self.R.getRunDetail('run')
        pretitle = 'r'+str(runnum)+'e'+str(evtnum)
        WFD = self.R.assocWFD(raw,pedsub=True)
        cd = sorted(WFD.keys())
        prenames = ['ped','pedsd','npulse','nsubp','area','time']
        debug = 0
        #if evtnum%100==3: debug = 1
        self.aWFD = {}
        for x in cd:
            ped,pedsd,iPulse,subPperP,pArea,pTime = self.W.pulseAnal(WFD[x],x,debug=debug)
            self.aWFD[x] = V = [ped,pedsd,iPulse,subPperP,pArea,pTime]
            y = float(cd.index(x))
            for tl in TL:
                for i,pren in enumerate(prenames):
                    name = pren+'_vs_WFD_'+tl
                    #print 'i,pren,V[i]',i,pren,V[i]
                    if pren=='ped' or pren=='pedsd':
                        self.Hists[name].Fill(V[i],y)
                    elif pren=='npulse':
                        self.Hists[name].Fill(float(len(V[i])),y)
                    else:
                        for v in V[i]:
                            self.Hists[name].Fill(float(v),y)
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
    def getRawDataList(self,rawDataDir):
        '''
        return list of full path name of files in rawDataDir
        cribbed from http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
        '''
        onlyfiles = [os.path.join(rawDataDir,f) for f in os.listdir(rawDataDir) if os.path.isfile(os.path.join(rawDataDir, f))]
        return onlyfiles
    def seeCalibData(self,rawDataDir):
        '''
        dump calib data for all runs in rawDataDir
        '''
        onlyfiles = self.getRawDataList(rawDataDir)
        for fn in onlyfiles:
            self.R.start(fn=fn)
            self.dumpCalib()
            self.R.closeHDF5File()
        return

if __name__ == '__main__' :
    nevt = big = 9999999999
    runstring = None
    dumpThres = -1.
    dumpAll = False
    if len(sys.argv)>1: nevt = int(sys.argv[1])
    if len(sys.argv)>2: runstring = str(sys.argv[2])
    if len(sys.argv)>3: dumpThres = float(sys.argv[3])
    if len(sys.argv)>4: dumpAll = True
    P = process()
    

    rawDataDir = "/Users/djaffe/work/1TonData/Filled_151217/"
    ######## enable to dump calib data for all files in rawDataDir
    if 0: 
        P.seeCalibData(rawDataDir)
        sys.exit('End of calibration data dump')

    
    fnlist = ["/Users/djaffe/work/1TonData/Filled_151217/run585.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run586.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run587.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run588.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run589.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run590.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run591.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run592.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run593.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run594.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run595.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run596.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run597.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run598.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run599.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run600.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run601.h5",
            "/Users/djaffe/work/1TonData/Filled_151217/run602.h5"]

    if runstring is not None:
        newlist = []
        for fn in fnlist:
            if runstring in fn: newlist.append(fn)
        if len(newlist)==0:
            # look for file in rawDataDir
            flist = P.getRawDataList(rawDataDir)
            for fn in flist:
                if runstring in fn: newlist.append(fn)
        if len(newlist)==0:
            sys.exit( 'ERROR Input runstring '+runstring+'yields zero matches')
            
        fnlist = newlist
        
    if nevt<big and runstring is None: fnlist = [fnlist[7]] # debug with one file

    # create name of output root file
    f = sorted(fnlist)
    f1,f2 = P.getFilePrefix(f[0]),P.getFilePrefix(f[-1])
    rfn = P.outputdir
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
        P.eventLoop(nevt,dumpAll=dumpAll)
        P.endRun()
        first = False
    P.finish(rfn)
