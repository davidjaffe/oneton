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
import Logger
import wfanal
import calibThermocouple
import writer
from optparse import OptionParser
import pipath
import SaveAllEvts

class process():
    def __init__(self,makeDirs=True):
        self.R = reader.reader()
        self.W = wfanal.wfanal()
        self.writer = writer.writer()
        self.cT = calibThermocouple.calibThermocouple()
        self.gU= graphUtils.graphUtils()
        self.pip= pipath.pipath()

        self.writeRecon = None
        self.overlaps = 0
        self.lastEventDumped = None
        
        self.Hists = {}
        self.TDChistnames = {}
        self.refTDCs = ['S2','H0']
        self.WFHists = []
        self.Times = []
        self.rawTemps = []
        self.calTemps = []
        self.AllTrigsfname = None
        self.rootpyEvts = None
        self.NoHists = False

        if makeDirs:
            now = datetime.datetime.now()
            self.start_time = now
            fmt = '%Y%m%d_%H%M%S_%f'
            cnow = now.strftime(fmt)
            parentDir = 'Results/'+cnow+'/'
            self.logdir = parentDir + 'Log/'
            self.figdir = parentDir + 'Figures/'
            self.WFfigdir = self.figdir + 'WF/'
            self.TDCfigdir= self.figdir + 'TDC/'
            self.outputdir= parentDir + 'Output/'
            # create list of platform-independent path
            dirs = self.pip.fix( [parentDir, self.logdir, self.figdir, self.WFfigdir, self.TDCfigdir, self.outputdir] )
            parentDir, self.logdir, self.figdir, self.WFfigdir, self.TDCfigdir, self.outputdir = dirs
            for d in dirs:
                if os.path.isdir(d):
                    pass
                else:
                    try:
                        os.makedirs(d)
                    except IOError,e:
                        print 'process__init__',e
                    else:
                        print 'process__init__ created',d
            lfn = self.pip.fix( self.logdir + '/' + cnow + '.log' )
            sys.stdout = Logger.Logger(fn=lfn)
            print 'process__init__ Output directed to terminal and',lfn
            print 'process__init__ Job start time',self.start_time.strftime('%Y/%m/%d %H:%M:%S')

        return
    def start(self, StoreAllTriggers=False):
        if not self.NoHists:
            self.bookHists('TDC')
            self.bookHists('QDC')
            self.bookHists('WFD')
        if StoreAllTriggers:
            self.InitTree()
            print('Initialized the tree.')
        return
    def InitTree(self):
        if self.AllTrigsfname is None:
            print('process.InitTree: the filename has not been defined!')
            return
        self.rootpyEvts = SaveAllEvts.rootpyEvts(self.AllTrigsfname)
        return
    def startRun(self,file=None,tG=False):
        '''
        start run processing
        report run details
        if writing to hdf5, set the run number and write parts of run record
        return True for successful start
        '''
        print 'process.startRun',file
        OK = self.R.start(fn=file)
        if OK: 
            self.R.reportRunDetail()
            if tG: self.R.testGet()
            if self.writeRecon:
                runnum = self.R.getRunDetail('run')
                self.writer.setRunNum(runnum)
                self.writer.writeRunData('RunNumber', self.R.getRunDetail('run') )
                self.writer.writeRunData('RunType',   self.R.getRunDetail('type') )
                self.writer.writeRunData('Material',  self.R.getRunDetail('Material') )
                self.writer.writeRunData('StartTime', self.R.getTime('Start_Time_str') )
                self.writer.writeRunData('Comments',  self.R.getRunDetail('Comments') )
        return OK
    def endRun(self):
        self.R.closeHDF5File()
        if self.overlaps>0:
            print 'process.endRun recorded',self.overlaps,'overlaps of lo,hi range of an ADC channel'
        self.overlaps = 0
        if self.rootpyEvts is not None:
            self.rootpyEvts.rfile.Flush()
        return
    def makeTvTgraphs(self):
        '''
        make graphs of raw and calib temperature vs time
        '''
        X = self.Times
        Y = self.rawTemps

        title = 'raw temp vs time'
        name = title.replace(' ','_')
        graw = self.gU.makeTGraph(X,Y,title,name)
        Y = self.calTemps
        title = 'calib temp vs time'
        name = title.replace(' ','_')
        gcal = self.gU.makeTGraph(X,Y,title,name)
        tmg = self.gU.makeTMultiGraph('temp_vs_time')
        tmg.Add(graw)
        tmg.Add(gcal)
        return [gcal,graw,tmg]
    def finish(self,rfn='oneton.root',Draw=True):
        '''
        execute at end of job
        Print summary output
        close hdf5 output file, if necessary
        Open and fill ROOT file
        Draw some hists, if requested
        Delete all hists/trees in memory
        Print job timing stats
        return name of hdf5 outfile, if it exists
        '''
        evtsSeen = self.R.summary()
        hdf5OutputFileName = None
        if self.writeRecon:
            hdf5OutputFileName = self.writer.closeFile()
        if self.rootpyEvts is not None:
            self.rootpyEvts.rfile.Write()
            self.rootpyEvts.rfile.Close()

        nwrite = 0
        if evtsSeen > 0:
            tvtgraphs  = self.makeTvTgraphs()
            tvtmg = tvtgraphs[-1]
            rf = TFile(rfn,"RECREATE")

            for h in self.Hists: rf.WriteTObject(self.Hists[h])
            nwrite += len(self.Hists)
            for h in self.WFHists: rf.WriteTObject(h)
            nwrite += len(self.WFHists)
            for g in tvtgraphs: rf.WriteTObject(g)
            nwrite += len(tvtgraphs)
            rf.Close()

            if Draw:
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

                # draw temp vs time
                if 0: #### DO NOT DRAW THESE GRAPHS. CAUSES PYTHON TO CRASH
                    for g in tvtgraphs[0:-1]:
                        self.gU.fixTimeDisplay(g,showDate=True)
                        self.gU.color(g,2,2)
                        self.gU.drawGraph(g,figDir=self.figdir)
                    self.gU.drawMultiGraph(tvtmg,figdir=self.figdir,xAxisLabel='Time',yAxisLabel='Temperature (C)')

            else:
                print 'process.finish Skip drawing hists and graphs'

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

        return hdf5OutputFileName
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
            if self.R.TDCchdesc is None:
                self.R.TDCchdesc = self.R.getModuleData('TDC','ChanDesc')
            md = self.R.TDCchdesc
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
            if self.R.QDCchdesc is None:
                self.R.QDCchdesc = self.R.getModuleData('QDC','ChanDesc')
            md = self.R.QDCchdesc
            clean = 'N/C' not in md
            while not clean:
                md.remove('N/C')
                clean = 'N/C' not in md
            for tl in TL:
                # 1-D hists
                for w in md:
                    name = tl + '_QDC_'+w
                    title = name.replace('_',' ')
                    nx = 4092
                    xmi = 0.
                    xma = xmi + float(nx)
                    #print 'process.bookHists name,title,nx,xmi,xma',name,title,nx,xmi,xma
                    Hists[name] = TH1D(name,title,nx,xmi,xma)
                    for r in ['lo','hi']: # lo=1,hi=0
                        name = tl+'_QDC'+r+'_'+w
                        title = name.replace('_',' ')
                        nx = 4092
                        if r=='lo': xmi,xma = 0.,float(nx)/8.
                        if r=='hi': xmi,xma = 0.,float(nx)
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
    def eventLoop(self,maxEvt=99999999,dumpAll=False,dumpThres=0.9999,
                  dumpNZcode=False,timeTempOnly=True,StoreAllTriggers=False):
        '''
        loop over events in order of increasing event number
        '''
        dumpOn = dumpThres>0 and not timeTempOnly
        dumpFirstEventOfEachRun = False
        
        CalData = self.R.getCalData()
        EvtDir  = self.R.getEvtDir()

        if dumpAll: self.dumpCalib()

        self.evtCode = {}

        if self.rootpyEvts is not None:
            self.rootpyEvts.tree.rnum = int(self.R.RunNumber)
            
        integerKeys = []
        for key in EvtDir.keys():
            if key!='evt_table': integerKeys.append(key)
        
        sEvtKeys = sorted(integerKeys,key=lambda b: int(b))
        for evtkey in sEvtKeys: # sorted(EvtDir.keys()):
            
            Event = EvtDir[evtkey]
            evtnum = int(evtkey)
            if evtnum>maxEvt: break

            if int(evtnum)%1000==0: print 'Event#',evtkey,'_________'

            dumpEvt = dumpAll or \
                (dumpOn and ( (int(evtnum)%10000==1 and dumpFirstEventOfEachRun) or numpy.random.uniform()>dumpThres))

            if dumpEvt:
                self.dumpEvent(Event,evtnum)


                            
            missing = self.R.unpackEvent(Event)
            if len(missing)==0: # OK
                self.analTvT(Event,evtnum)

                if not timeTempOnly:
                    triggers = self.R.unpackTrigger(Event['TDC'])
                    self.evtCode['TDC'] = tdcCode = self.analTDC(Event['TDC'])
                    self.evtCode['QDC'] = qdcCode = self.analQDC(Event,evtnum)
                    self.evtCode['WFD'] = wfdCode = self.analWFD(Event,evtnum)

                    if dumpAll or (dumpNZcode and wfdCode>0):
                        print 'process.eventLoop evtnum',evtnum,'wfdCode',wfdCode
                        self.dumpaWFD()
                        self.evtDisplay(Event,evtnum)  # must be called after analWFD
                        self.dumpEvent(Event,evtnum)
                
                  
            if self.writeRecon:
                dets = ['TDC','QDC']
                for det,DICT in zip(dets,[self.aTDC,self.aQDC]):
                    ks = sorted(DICT.keys())
                    #print 'process.eventLoop det,ks',det,ks
                    labels,data = [],[]
                    for key in ks:
                        labels.append(  det + '/' + key )
                        data.append(  DICT[key] )
                    EN = self.writer.writeEvent(labels,data,evtNo=evtnum)
                # WFD data structure more complicated
                # suppress storing of data when no pulses found
                det = 'WFD'
                ks = sorted(self.aWFD.keys())
                del ks[ ks.index('prenames') ]
                del ks[ ks.index('RunEventNumber') ]
                prenames = self.aWFD['prenames']
                use = {}
                ipn = prenames.index('time')
                for key in ks:
                    use[key] = len(self.aWFD[key][ipn])>0
                labels,data = [],[]
                for ipn,pn in enumerate(prenames):
                    for key in ks:
                        if use[key]:
                            labels.append( det + '/' + pn + '/' + key )
                            data.append(self.aWFD[key][ipn])
                EN = self.writer.writeEvent(labels,data,evtNo=evtnum)
                label,data = 'Time',self.aTime                
                EN = self.writer.writeEvent(label,data,evtNo=evtnum)
                label,data = 'Temp',self.aTemp
                EN = self.writer.writeEvent(label,data,evtNo=evtnum)
            
            if self.rootpyEvts is not None:
                thetree = self.rootpyEvts.tree
                thetree.evtnum = int(evtkey)
                thetree.temp = self.cT.getCalibTC(Event['Event_Temp'].value)
                thetree.time = Event['Event_Time'].value
                thetree.QDC1 = [ch[1] for ch in Event['QDC_1']]
                thetree.QDC2 = [ch[1] for ch in Event['QDC_2']]
                thetree.TDCnumch = len(Event['TDC'])
                thetree.scaler = [int(val) for val in Event['Scaler']]
                trigstring = ''                
                for trig in triggers: trigstring += trig + ' '
                thetree.trigtype = trigstring.ljust(9)
                
                for i in range(len(Event['TDC'])):
                    thetree.TDCch[i] = Event['TDC'][i][0]
                    thetree.TDCval[i] = Event['TDC'][i][1]
                
                if not timeTempOnly:
                    #note: channel names are assumed to be sampled from:
                    # ['S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7']
                    #self.aWFD[channel name] goes like prenames:
                    #[pedestal, pedestal stdev, npulse, nsubp, area, time]
                    #I'm not sure what npulse and nsubp are, so I'm skipping them.
                    thekeys = [k for k in self.aWFD.keys()
                                if k is not 'RunEventNumber'
                                and k is not 'prenames']
                    for k in thekeys:
                        areas = [0]*32
                        times = [0]*32
                        theidx = int(k[1])
                        thetree.WFDped[theidx] = self.aWFD[k][0]
                        thetree.WFDpedsd[theidx] = self.aWFD[k][1]
                        for i in range(len(self.aWFD[k][4])):
                            if i>31:
                                print("WARNING: tried to save event with > 32 pulses. Only saving first 32...")
                                therun, theevt = self.aWFD['RunEventNumber']                                
                                print("Run = {0}, Event = {1}, channel = {2}".format(therun, theevt, k))
                                break
                            areas[i] = self.aWFD[k][4][i]
                            times[i] = self.aWFD[k][5][i]
                        setattr(thetree, k+'_area', areas)
                        setattr(thetree, k+'_time', times)

                self.rootpyEvts.tree.Fill()
                self.reinitROOTf()                    
        return

    def reinitROOTf(self):
        '''
        This function reinitializes the root tree and should be called after
        filling the tree.
        '''
        thetree = self.rootpyEvts.tree
        thetree.TDCch = [0]*16
        thetree.TDCval = [0]*16
        thetree.QDC1 = [0]*16
        thetree.QDC2 = [0]*16
        thetree.WFDped = [0]*8
        thetree.WFDpedsd = [0]*8        
        for i in range(8):
            prefix = 'S{0}_'.format(i)
            setattr(thetree, prefix + 'area', [0]*32)
            setattr(thetree, prefix + 'time', [0]*32)
        return

    def analTvT(self,raw,evtnum):
        '''
        analysis time, temperature data
        '''
        ti = self.R.unpackTime(raw)
        self.aTime = ti
        rawtemp = self.R.unpackTemperature(raw)
        caltemp = self.cT.getCalibTC(rawtemp)
        self.aTemp = [rawtemp, caltemp]
        self.Times.append(ti)
        self.rawTemps.append(rawtemp)
        self.calTemps.append(caltemp)
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
        '''
        dump trigger info, unpacked data for an event
        avoid duplicate dumps
        '''
        if self.lastEventDumped is not None:
            if self.lastEventDumped==evtnum:
                print 'process.dumpEvent evt#',evtnum,'has already been dumped.'
                return
        self.lastEventDumped = evtnum

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

            
        H = []
        self.hWFD = {}
        mostPulses = -1
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
            mostPulses = max(mostPulses,len(self.aWFD[cd][2]))
                

        # increase size of x,y labels
        for h in H:
            for a in ["x","y"]:
                s = h.GetLabelSize(a)
                h.SetLabelSize(2.*s,a)
        s = 'run'+str(self.R.getRunDetail('run'))+'evt'+str(evtnum)

        self.diagnoseWFD(fname=s,figdir=self.WFfigdir,zoomPulse=0) # add 'zoom' on first pulse
        if 'WFD' in self.evtCode:
            if self.evtCode['WFD']>0 and mostPulses>-1:
                for zP in range(mostPulses):
                    self.diagnoseWFD(fname=s,figdir=self.WFfigdir,zoomPulse=zP)
                    
        return
    def diagnoseWFD(self,fname='diagnoseWFD',figdir='',zoomPulse=-1):
        '''
        overlay histograms of waveforms with results of wfanal
        if zoomPulse is valid index of pulse start, stop then zoom in around pulse
        '''
        md = sorted(self.hWFD.keys())
        DICT = self.hWFD

        # arcane bullshit to draw lines...???
        lines = {}
        for icd,cd in enumerate(md):

            V = self.aWFD[cd] #### [ped,pedsd,iPulse,subPperP,pArea,pTime]
            h = DICT[cd]

            ymi,yma = h.GetMinimum(),h.GetMaximum()
            xmi,xma = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
            h.SetStats(0)
            ped = V[0]
            L =  ROOT.TLine(xmi,ped,xma,ped)
            L.SetLineColor(ROOT.kGreen) # solid green = global baseline
            lines[cd] = [ L ]
            pedsd = V[1]
            for x in [-1.,1.]:
                y = ped+x*pedsd
                L = ROOT.TLine(xmi,y,xma,y)
                L.SetLineColor(ROOT.kGreen) # dashed green = +- rms on global baseline
                L.SetLineStyle(3)
                lines[cd].append( L )
            iPulse = V[2]
            for pair in iPulse:
                for ib in pair:
                    L = ROOT.TLine(float(ib),ymi,float(ib),yma)
                    L.SetLineColor(ROOT.kBlue) # dotted blue = start, end of pulse
                    L.SetLineStyle(3) # dotted
                    lines[cd].append( L )
            pTime = V[5]
            for t in pTime:
                L = ROOT.TLine(t,ymi,t,yma)
                L.SetLineColor(ROOT.kRed)   # dashed red = fitted time of pulse
                L.SetLineStyle(2) # dashed
                lines[cd].append( L )

        path = fname #+ '_'+str(i)
        if zoomPulse>-1: path += 'zoomPulse'+str(zoomPulse)
        if figdir!='':
            if figdir[-1] != os.path.sep:
                path = self.pip.fix(figdir + '/' + path)
            else:
                path = figdir + path
        ps = path + '.ps'
        pdf= path + '.pdf'
        ROOT.gROOT.ProcessLine("gROOT->SetBatch()") # no pop up
        xsize,ysize = 1100,850 # landscape style
        canvas = ROOT.TCanvas(pdf,fname,xsize,ysize)
        canvas.Divide(2,4)
        for icd,cd in enumerate(md):
            canvas.cd(icd+1)
            h = DICT[cd]

            # zoom?
            xlo=xhi=0
            if zoomPulse>-1:
                V = self.aWFD[cd] #### [ped,pedsd,iPulse,subPperP,pArea,pTime]
                iPulse = V[2]
                if len(iPulse)>zoomPulse:
                    xlo = max(0,iPulse[zoomPulse][0]-30)
                    xhi = min(iPulse[zoomPulse][1]+100,4092)
            h.GetXaxis().SetRange(xlo,xhi)

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
            if os.path.exists(pdf): os.remove(ps)
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
        return code = 0 = ok, otherwise not
        '''
        code = 0
        runnum = self.R.RunNumber
        moniker = 'Run ' + str(runnum) + ' Event ' + str(evtnum)
        TL = ['All']
        TL.extend(self.R.unpackTrigger(raw['TDC']))

        if self.R.QDCchdesc is None:
            self.R.QDCchdesc = self.R.getModuleData('QDC','ChanDesc')
        md = self.R.QDCchdesc
        while 'N/C' in md: md.remove('N/C')
        
        # create a dictionary with key = channel, value = list of QDC values (can have 2 entries, 1 for lo- and 1 for hi-range)
        QDC = self.R.assocQDC(raw)
        cd = sorted(QDC.keys())
        rng = ['hi','lo'] # low,high range: hilo =1,0
        self.aQDC = {}
        for x in cd:
            if len(QDC[x]) == 1:
                self.aQDC[x] = QDC[x]
            elif len(QDC[x]) == 2:
                #probably both ranges triggered, use the low one unless there is an overflow.
                val = [data for data in QDC[x] if data[2]==1 and not data[3]==1]
                if len(val) is not 0:
                    self.aQDC[x] = val[0]
                else:
                    self.aQDC[x] = [data for data in QDC[x] if data[2] != 1]
            else:
                print('WARNING: more than 2 entries in a single QDC channel! Writing out zeroes to channel -1')
                QDC[x][0] = (-1,0,0,0)
                self.aQDC[x] = QDC[x][0]
            i = md.index(x)
            for y in QDC[x]:
                #print 'process.analQDC x,QDC[x],y:',x,QDC[x],y
                if len(y)==3: # early data
                    ch,v,hilo = y
                else:
                    ch,v,hilo,ovfl = y # chan#, value, hirange or lowrange, overflow
                if not self.NoHists:
                    for tl in TL:
                        name = tl+'_QDC'+rng[hilo]+'_'+x
                        self.Hists[name].Fill(v)
                        name = tl+'_QDC_'+x
                        self.Hists[name].Fill(v)
                        name = 'chan_vs_QDC_'+tl
                        self.Hists[name].Fill(v,float(i))
        
        for x in QDC:
            if len(QDC[x])>1 and x!='N/C':
                self.overlaps += 1
                if self.overlaps<=5:
                    print 'process.QDC',moniker,'overlap',x,QDC[x]
                if self.overlaps==5:
                    print 'process.QDC Report of further overlaps will be suppressed'
                    
        return code
        
    def analWFD(self,raw,evtnum):
        '''
        unpack and pedestal subtract waveforms
        return code = 0 = all ok, != 0 => not ok
        '''
        display = False

        code = 0

        TL = ['All']
        TL.extend(self.R.unpackTrigger(raw['TDC']))
                
        runnum = self.R.getRunDetail('run')
        pretitle = 'r'+str(runnum)+'e'+str(evtnum)
        WFD = self.R.assocWFD(raw,pedsub=False)#pedsub=True)
        cd = sorted(WFD.keys())
        prenames = ['ped','pedsd','npulse','nsubp','area','time']
        debug = 0
        #if evtnum%100==3: debug = 1
        self.aWFD = {}
        self.aWFD['RunEventNumber'] = [runnum, evtnum]
        self.aWFD['prenames'] = prenames
        for x in cd:
            ped,pedsd,iPulse,subPperP,pArea,pTime = self.W.pulseAnal(WFD[x],x,debug=debug)
            self.aWFD[x] = V = [ped,pedsd,iPulse,subPperP,pArea,pTime]
            for K in iPulse:
                if K[0]<0 or K[1]<0:
                    code = 1
                    print 'process.analWFD ',x,'setting code',code,'iPulse',iPulse,'K',K
            for t in pTime:
                if t<0. :
                    code = 2
                    print 'process.analWFD ',x,'setting code',code,'pTime',pTime
            #if code>0: a1,a2,a3,a4,a5,a6 = self.W.pulseAnal(WFD[x],x,debug=1)   ###### SPECIAL #####
            y = float(cd.index(x))
            if not self.NoHists:
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
        return code
    def dumpaWFD(self):
        '''
        dump current contents of results of waveform analysis
        '''
        run,evt = None,None
        sREN = 'RunEventNumber'
        run,event = self.aWFD[sREN]
        print 'process.dumpaWFD run',run,'event',event,' +++++++++++++++++'
        spn = 'prenames'
        prenames = self.aWFD[spn]
        
        for x in self.aWFD:
            if x!=sREN and x!=spn:
                print x,
                for w,v in zip(prenames,self.aWFD[x]): print w,v,
                print ''
        print 'process.dumpaWFD ++++++++++++++++++++++++ end'
        return
    def analTDC(self,raw):
        '''
        analyze TDC data for a single event
        return code = 0 = ok, otherwise not
        '''
        code = 0
        if self.R.TDCchdesc is None:
            self.R.TDCchdesc = self.R.getModuleData('TDC','ChanDesc')
        TDCmd = self.R.TDCchdesc
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
        self.aTDC = {}
        for x,y in zip(cd,TDCns) :
            nsCD[x] =y
            self.aTDC[x] = y

        TDChits = []

        if not self.NoHists:
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

        return code
    def getFilePrefix(self,fn):
        '''
        get file name prefix from full path name
        20160304 make platform independent
        '''
        #f = fn.split('/')[-1]
        f = os.path.basename(fn) 
        return os.path.splitext(f)[0]#Won't fail on periods in the filename.
    def getRawDataList(self,rawDataDir, H5Files=False, excldirs=[]):
        '''
        return list of full path name of files in rawDataDir
        cribbed from http://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
        and http://stackoverflow.com/questions/13571134/python-how-to-recursively-go-through-all-subdirectories-and-read-files

        Arguments:
        - Base directory.
        - Flag to use HDF5 files or zip files.
        - A list of directories to exclude from the search.
        '''
        if H5Files:
            ext = '.h5'
        else:
            ext = '.zip'
        datafiles = []
        for folder, subfolders, files in os.walk(rawDataDir, topdown=True):
            subfolders[:] = [d for d in subfolders if d not in excldirs]
            for thefile in files:
                if os.path.splitext(thefile)[1]==ext and 'RunInfo.h5' not in thefile:
                    datafiles.append(os.path.join(folder, thefile))
        return datafiles
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
    def defaultFileList(self):
        '''
        predetermined list of input files
        '''
        fnlist = self.pip.fix ( ["/Users/djaffe/work/1TonData/Filled_151217/run585.h5",
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
                "/Users/djaffe/work/1TonData/Filled_151217/run602.h5"] )
        return fnlist

if __name__ == '__main__' :
    big = 9999999999
    parser = OptionParser()
    parser.add_option("-N","--Nevents",default=big,type=int,
                      help="Number of events per input file to process. [default %default]")
    parser.add_option("-C","--MakeInputFileList",action="store_true",
                      help="Make input file list from all appropriate files in raw data directory")
    parser.add_option("-s","--RunString",default=None,type=str,
                      help="String to match when selecting input files. If None, then use default file list [Default %default]")
    parser.add_option("-O","--OneRunOnly",action="store_true",
                      help="Process only the first run in the filelist")
    parser.add_option("-t","--DumpThreshold",default=0.9999,type=float,
                      help="Dump threshold. If random>threshold then dump event [default %default]")
    parser.add_option("-A","--DumpAll",action="store_true",
                      help="Dump every event")
    parser.add_option("-X","--DumpNonZeroWFDcode",action="store_true",
                      help="Dump events with non-zro WFDcode")
    parser.add_option("-W","--WriteRecon",action="store_true",
                      help="Write out reconstruction/calibrated events to file")
    parser.add_option("-T","--TemperatureVsTime",action="store_true",
                      help="Just do temperature vs time analysis")
    parser.add_option("-Z","--UseCompression",default=None,type=str,
                      help="Use compression in creation of output hdf5 file. options=lzf,gzip [default %default]")
    parser.add_option("--DoSPECalculation",action="store_true",
                      help="Do SPE calculation")
    parser.add_option("-F", "--H5Files", action="store_true", dest="H5Files", default=False,
                      help="Option to process .h5 files (rather than zipped HDF5 files)")
    parser.add_option("-R", "--StoreAllTriggers", action="store_true",
                      dest="StoreAllTriggers", default=False,
                      help="Optionally store all triggers as a tree in the \
                      ROOT file. This can make the output VERY big!")
    parser.add_option("-E", "--ExcludeDirs", action="store", dest="excldirs",
                      default = "", type="string",
                      help="Exclude sub-dirs from the data file search if using -R. \
                      Use paths to the base search directory, separated by a semicolon.")
    parser.add_option('-H', '--NoHists', action='store_true',
                      help='Don\'t calculate any histograms')

    (options, args) = parser.parse_args(args=sys.argv)
    #print 'options',options
    nevt = options.Nevents
    dumpAll = options.DumpAll
    dumpThres=options.DumpThreshold
    dumpNZcode=options.DumpNonZeroWFDcode
    compalg = options.UseCompression

    Draw = True
    if options.DoSPECalculation:
        print 'process.main Do SPE Calculation. Three stage job'
        Draw = False # drawing sometimes causes seg fault
    
    
    P = process()
    P.writeRecon = options.WriteRecon
    P.NoHists=options.NoHists
    # get list of potential input files
    if len(args) > 1:
        rawDataDir = P.pip.fix(args[1])
    else:
        rawDataDir = os.getcwd()#Get current directory
    if options.MakeInputFileList: 
        excldirs = options.excldirs.split(";")
        fnlist = P.getRawDataList(rawDataDir, options.H5Files, excldirs)
    else:
        fnlist = P.defaultFileList()

    # winnow list of input files. 
    if options.RunString is not None:
        newlist = []
        for fn in fnlist:
            if options.RunString in fn: newlist.append(fn)
        if len(newlist)==0:
            sys.exit("process.main ERROR RunString " + options.RunString + " produced zero length file list")
        fnlist = newlist
    if options.OneRunOnly: fnlist = [fnlist[0]]

    print 'fnlist',fnlist
        
    # create name of output root file
    f = sorted(fnlist)
    f1,f2 = P.getFilePrefix(f[0]),P.getFilePrefix(f[-1])
    rfn = P.outputdir
    if len(rfn)>0 and rfn[-1]!=os.path.sep : rfn += os.path.sep # ensure dir name ends in separator
    #Make filename for output    
    g = [f1]
    if f1!=f2: g.append(f2)
    for q in g:
        rfn += q
        if len(g)>1 and '-' not in rfn: rfn += '-'
    outputFilePrefix = rfn
    rfn += '.root'
    if options.StoreAllTriggers:
        P.AllTrigsfname = outputFilePrefix + '_alltrigs.root'
    print 'process.main Output root file name is',rfn
    reconfn = outputFilePrefix + '.h5'
    if options.WriteRecon:
        print 'process.main Output recon file name is',reconfn
        P.writer.openFile(reconfn,compalg=compalg)


    # begin processing
    first = True
    for fn in fnlist:
        if P.startRun(file=fn):
            if first: P.start(StoreAllTriggers=options.StoreAllTriggers)
            P.eventLoop(nevt,dumpAll=dumpAll,dumpThres=dumpThres,
                        timeTempOnly=options.TemperatureVsTime,
                        dumpNZcode=dumpNZcode,
                        StoreAllTriggers=options.StoreAllTriggers)
            P.endRun()
            first = False
    hdf5OutputFileName = P.finish(rfn=rfn, Draw=Draw)

    print 'back in __main__'

    # begin second stage of processing
    if options.DoSPECalculation:
        print '\n\n Begin second stage of processing for SPE calculation \n'
        import second
        S = second.second(useLogger=False)
        outrfn = S.main(maxevt=nevt,LEDonly=True,inputFile=hdf5OutputFileName)

        print '\n\n Begin third stage of processing for SPE calculation \n'
        import spe_calc
        SC = spe_calc.spe_calc(inputRFN=outrfn)
        SC.main(drawEachFit=True,singleRunMode=True)

        

    sys.exit('Normal exit')  ############EXIT EXIT EXIT ######################
