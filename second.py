#!/usr/bin/env python
'''
process recon/calib 1ton data
20160208
'''
import h5py
import numpy
import graphUtils
import ROOT
from ROOT import TFile,TH1D,TH2D,gDirectory
import sys
import datetime
#import reader
import math
import os
import Logger
#import wfanal
#import calibThermocouple
import writer
#from optparse import OptionParser
import process
import gzip,shutil

class second():
    def __init__(self):
        self.f = None
        self.writer = writer.writer()
        self.P = process.process(makeDirs=False)
        self.gU= graphUtils.graphUtils()

        now = datetime.datetime.now()
        self.start_time = now
        fmt = '%Y%m%d_%H%M%S_%f'
        cnow = now.strftime(fmt)
        self.parentDir = 'Second/'+cnow+'/'
        dirs = [self.parentDir]
        for d in dirs:
            if os.path.isdir(d):
                pass
            else:
                try:
                    os.mkdir(d)
                except IOError,e:
                    print 'second__init__',e
                else:
                    print 'second__init__ created',d

        # set time range for hists. taken from log files
        self.tStart = '12:05:49 01/11/16' # run537
        self.tEnd   = '23:59:59 02/22/16' # in future...
        self.timeFormat = '%H:%M:%S %m/%d/%y'

                                        

        self.CTtimes = {} # for WFD
        for c in ['S0','S1','S2','S3','S4','S5','S6','S7']:
            self.CTtimes[c] = [120.,200.]
            if c in ['S4','S5']: self.CTtimes[c] = [90.,160.]

        # used in getHodosInCoinc
        self.goodCoincNsigma = 4.
        self.goodCoinc = {}
        dTrange = {} # wrt H0
        dTrange['H1'] = [-6.54, 0.96 ] # mean, rms
        dTrange['H2'] = [-6.20, 1.12 ] # mean, rms
        dTrange['H3'] = [-5.07, 1.08 ] # mean, rms
        dTrange['H4'] = [-6.81, 3.46 ] # mean, rms
        dTrange['H5'] = [-6.77, 2.97 ] # mean, rms
        self.goodCoinc['CT'+'H0'] = dTrange
        dTrange = {} # wrt H2
        dTrange['H0'] = [ 6.20, 1.12] # mean, rms
        dTrange['H1'] = [-0.425,1.22] # mean, rms
        dTrange['H3'] = [ 1.14, 1.23] # mean, rms
        dTrange['H4'] = [-0.62, 3.80] # mean, rms
        dTrange['H5'] = [-0.43, 2.99] # mean, rms
        self.goodCoinc['CT'+'H2'] = dTrange

        # used in goodHodoTiming
        self.goodHodo = {}
        TDCrange = {}
        TDCrange['H0'] = [70.,85.] # upper, lower limit
        TDCrange['H1'] = [65.,75.]
        TDCrange['H2'] = [65.,80.]
        TDCrange['H3'] = [70.,75.]
        TDCrange['H4'] = [60.,85.]
        TDCrange['H5'] = [65.,85.]
        self.goodHodo['CT'] = TDCrange
        
        #### used in getTrajectories
        T = {}
        for hi in ['H0','H2','H0.H2']:
            for mid in ['H1','H3','H1.H3']:
                c = hi + '.' + mid
                T[c] = c.split('.')
                for lo in ['H4','H5','H4.H5']:
                    c = hi + '.' + mid + '.' + lo
                    T[c] = c.split('.')
        self.Trajectories = T

        sKeys = sorted(self.Trajectories.keys())
        #print 'second__init__ Trajectories sKeys sort'
        #for c in sKeys:
        #    print c,self.Trajectories[c]
        #print ''

        tKeys = sorted(sKeys,key = lambda b:len(b))
        if 0:
            print 'second__init__ Trajectories tKeys sort'
            for c in tKeys:
                print c,self.Trajectories[c]
            print ''
        self.tKeys = tKeys
        
        
        return
    def open(self,fn):
        '''
        open hdf5 file. handle gzipped files
        return True if file was opened successfully
        '''
        OK = True
        bn = os.path.basename(fn)
        bnsuf = bn.split('.')[-1]
        if bnsuf=='gz':
            h5f = fn.replace('.gz','')
            print 'second.open gunzip',fn,'to',h5f
            with gzip.open(fn,'rb') as fin, open(h5f,'wb') as fout:
                    shutil.copyfileobj(fin,fout)
            self.f = h5py.File(h5f,'r')
            self.gzf = fn
        elif bnsuf=='h5':
            self.gzf = None
            self.f = h5py.File(fn,'r')
        else:
            OK = False
            print 'second.open ERROR processing ' + fn +  ' UNKNOWN suffix ' + bnsuf
  
        if OK : print 'second.open Opened',self.f.filename
        return OK
    def close(self):
        '''
        close hdf5 file and delete it if gzipped version exists
        '''
        fn = self.f.filename
        self.f.close()
        words = 'second.close closed ' + fn
        if self.gzf is not None:
            if os.path.isfile(fn) and os.path.isfile(self.gzf):
                os.remove(fn)
                words += ' and deleted it since gzipped version exists'
        print words

        return
    def loop(self,maxevt=99999999):
        '''
        loop over all events in file
        up to event number maxevt
        '''
        for run in self.f['Run']:
            runnum = int(run)
            print 'second.loop run',runnum
            Events = self.f['Run/'+run+'/Event']
            for event in Events:
                evtnum = int(event)
                if evtnum%1000==0  : print 'second.loop event',evtnum
                if evtnum>maxevt: break
                data = Events[event]
                TDC = self.getTDC(data['TDC'])
                QDC = self.getQDC(data['QDC'])
                WFDtime = WFDarea = WFDnp = WFDped = WFDpsd = None
                if 'WFD' in data: WFDtime,WFDarea,WFDnp,WFDped,WFDpsd = self.getWFD(data['WFD'])
                Time = data['Time']
                Temp = data['Temp']
                self.analyze(TDC,QDC,WFDtime,WFDarea,WFDnp,WFDped,WFDpsd,Time,Temp)
        return
    def book(self):
        self.Hists = {}
        tS = self.gU.getTDatime(self.tStart,self.timeFormat)
        tE = self.gU.getTDatime(self.tEnd,self.timeFormat)
        for itrig,trig in enumerate(['CT','M','LED']):
            for x in ['S0','S1','S2','S3','S4','S5','S6','S7']:
                name,nx,xmi,xma = trig + '_WFD_dt_'+x,100,0.,400.
                title = name.replace('_',' ')
                self.Hists[name] = TH1D(name,title,nx,xmi,xma)

                name,nx,xmi,xma = trig + '_WFD_time_'+x,100,0.,400.
                title = name.replace('_',' ')
                self.Hists[name] = TH1D(name,title,nx,xmi,xma)

                name,nx,xmi,xma = trig + '_WFD_area_'+x,200,0.,2000.
                title = name.replace('_',' ')
                self.Hists[name] = TH1D(name,title,nx,xmi,xma)

                name,nx,xmi,xma = trig + '_WFD_Area_tcut_'+x,200,0.,2000.
                title = name.replace('_',' ')
                self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                if trig=='CT':
                    nx,xmi,xma = len(self.tKeys),-0.5,-0.5+float(len(self.tKeys))
                    ny,ymi,yma = 50,0.,1000.
                    name = trig + '_WFD_At_vs_Traj_'+x
                    title = name.replace('_',' ')
                    self.Hists[name] = TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
                    for i,t in enumerate(self.tKeys):
                        self.Hists[name].GetXaxis().SetBinLabel(i+1,t)
                    self.Hists[name].LabelsOption("v","X")

                    nx,xmi,xma = 100,tS,tE
                    ny,ymi,yma = 50,0.,1000.
                    name = trig + '_WFD_At_vs_time_'+x
                    title = name.replace('_',' ')
                    self.Hists[name] = TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
                    self.Hists[name].GetXaxis().SetCanExtend(False)
                    self.Hists[name].GetYaxis().SetCanExtend(False)
                    ##self.gU.reportHist(self.Hists[name])  #### FOR DEBUG

                if itrig==0: # trigger-independent figures
                    name,nx,xmi,xma = 'npulse_'+x,21,-0.5,20.5
                    title = name.replace('_',' ')
                    self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                    name,nx,xmi,xma = 'pedestal_'+x,50,-5.,5.
                    title = name.replace('_',' ')
                    self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                    name,nx,xmi,xma = 'ped_stdev_'+x,50,0.,5.
                    title = name.replace('_',' ')
                    self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                    
                    
            if 'LED'!=trig:
                for Href in ['H0','H2']:
                    for H in ['H0','H1','H2','H3','H4','H5']:
                        if H!=Href:
                            name,nx,xmi,xma = trig+'_dT_'+H+'wrt'+Href,100,-20.,20.
                            title = name.replace('_',' ')
                            self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                        if Href=='H0':
                            name,nx,xmi,xma = trig+'_TDC_'+H,100,0.,200.
                            title = name.replace('_',' ')
                            self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                            
        # distribution of hodo coincidences
        nx = len(self.tKeys)
        xmi = -0.5
        xma = xmi + float(xmi)
        name = title = 'CT_Trajectories'
        self.Hists[name] = TH1D(name,title,nx,xmi,xma)
        for i,t in enumerate(self.tKeys):
            self.Hists[name].GetXaxis().SetBinLabel(i+1,t)
        self.Hists[name].LabelsOption("v","X")

        
        
        
        print 'second.book booked hists'
        return
    def getTriggers(self,TDC):
        '''
        return list of strings giving triggers that fired
        '''
        triggers = []
        for t in ['CT','M']:
            if t in TDC: triggers.append(t)
        if len(triggers)==0 and 'CT+M+LED' in TDC: triggers.append('LED')
        return triggers
            
    def analyze(self,TDC,QDC,WFDtime,WFDarea,WFDnp,WFDped,WFDpsd,Time,Temp):
        '''
        analysis
        get triggers based on TDC info
        obtain trajectories defined by hodoscopes for cosmic triggers
        plot WFD time and area
        '''
        triggers = self.getTriggers(TDC)
        
        #print 'second.analysis Time',Time[()],
        dt = datetime.datetime.fromtimestamp(Time[()])
        #print 'dt',dt,
        xt = datetime.datetime.strftime(dt,self.timeFormat)
        #print 'xt',xt,
        evtTime = self.gU.getTDatime(xt,fmt=self.timeFormat)
        #print 'evtTime',evtTime
        
        Traj = []
        if 'CT' in triggers:
            if 0:
                print '----- TDCs',[x+' '+str(TDC[x][1]) for x in TDC]
                Hits = []
                for x in TDC:
                    if 'H' in x: Hits.append(x)
                print 'TDC hodos',Hits
                Hodos = self.goodHodoTiming(TDC,trig='CT')
                print 'good timing',Hodos
            Hodos = self.getHodosInCoinc(TDC,trig='CT')
            if 0: print 'in coinc',Hodos
            Traj  = self.getTrajectories(Hodos)
            if 0: print 'trajs',Traj
            name = 'CT_Trajectories'
            for t in Traj:
                i = self.tKeys.index(t)
                self.Hists[name].Fill(float(i))

        
        if WFDtime is not None:
            #print 'second.analyze WFDnp,WFDped,WFDpsd',WFDnp,WFDped,WFDpsd
            
            for x in WFDnp:
                name = 'npulse_'+x
                #print name,len(WFDnp[x]),
                self.Hists[name].Fill(float(len(WFDnp[x])))
            #print ''
            for x in WFDped:
                name = 'pedestal_'+x
                #print name,WFDped[x],
                self.Hists[name].Fill(WFDped[x][()])
            #print ''
            for x in WFDpsd:
                name = 'ped_stdev_'+x
                #print name,WFDpsd[x],
                self.Hists[name].Fill(WFDpsd[x][()])
            #print ''

            
            for trig in triggers:
                for x in WFDtime:
                    ts = WFDtime[x]
                    for t in ts:
                        name = trig+'_WFD_time_'+x
                        self.Hists[name].Fill(t)
                    if len(ts)>1:
                        name = trig+ '_WFD_dt_'+x
                        for i,t in enumerate(ts):
                            if i+1<len(ts):
                                dt = ts[i+1]-t
                                self.Hists[name].Fill(dt)
                for x in WFDarea:
                    for areas in WFDarea[x]:
                        name = trig+'_WFD_area_'+x
                        for a in self.makeList(areas):
                            self.Hists[name].Fill(abs(a))
                if trig=='CT':
                    for x in WFDtime:
                        for ts,areas in zip(WFDtime[x],WFDarea[x]):
                            for t,a in zip(self.makeList(ts),self.makeList(areas)):
                                if self.CTtimes[x][0]<=t and t<=self.CTtimes[x][1]:
                                    name = trig + '_WFD_Area_tcut_'+x
                                    self.Hists[name].Fill(abs(a))
                                    for Tj in Traj:
                                        name = trig + '_WFD_At_vs_Traj_'+x
                                        i = self.tKeys.index(Tj)
                                        self.Hists[name].Fill(float(i),abs(a))

                                    name = trig + '_WFD_At_vs_time_'+x
                                    self.Hists[name].Fill(evtTime,abs(a))

                                
        for trig in triggers:
            if 'LED'!=trig:
                for Href in ['H0','H2']:
                    tref = None
                    if Href in TDC: tref = TDC[Href][1]
                    for H in ['H0','H1','H2','H3','H4','H5']:
                        if H in TDC:
                            t = TDC[H][1]
                            #print 'H,t',H,t
                            if Href=='H0':
                                name = trig+'_TDC_'+H
                                self.Hists[name].Fill(t)
                            if H!=Href and tref is not None:
                                name = trig+'_dT_'+H+'wrt'+Href
                                self.Hists[name].Fill(t-tref)
        return
    def getTrajectories(self,Hodos):
        '''
        return trajectories based on coincident Hodos
        '''

        #print 'second.getTrajectories:Hodos',Hodos
        T = []
        for c in self.Trajectories:
            #print 'second.getTrajectories:c,T[c]',c,self.Trajectories[c]
            #print 'second.getTrajectories:set(self.Trajectories[c])<=set(Hodos)',set(self.Trajectories[c]) <= set(Hodos)
            if set(self.Trajectories[c]) <= set(Hodos) : T.append(c)
        return T
    def goodHodoTiming(self,TDC,trig='CT'):
        '''
        return list of hodoscopes with 'good' timing for trig
        '''

        Hodos = []
        if trig in self.goodHodo:
            TDCrange = self.goodHodo[trig]
            for x in TDC:
                t = TDC[x][1]
                if x in TDCrange:
                    tlo,thi = TDCrange[x]
                    if self.btwn(t,tlo,thi): Hodos.append( x )
        return Hodos
    def getHodosInCoinc(self,TDC,trig='CT'):
        '''
        return list of hodoscopes in coincidence
        require each hodo to be in 'good' timing range for the input trigger
        then apply requirements on the relative timing wrt H0 or H2
        '''
        
        Hodos = self.goodHodoTiming(TDC,trig=trig)
        nsig = self.goodCoincNsigma
        for href in ['H0','H2']:
            if href in Hodos:
                key = trig+href
                tref = TDC[href][1]
                dTrange = self.goodCoinc[key]
                #print 'second.getHodosInCoinc: href,tref',href,tref
                for ih,h in enumerate(Hodos):
                    if h!=href:
                        t = TDC[h][1]
                        mean,rms = dTrange[h]
                        #print 'second.getHodosInCoinc: h,t,t-tref,mean,rms,nsig,near',h,t,t-tref,mean,rms,nsig,self.near(t-tref,mean,rms,nsig)
                        if not self.near(t-tref,mean,rms,nsig): del Hodos[ih]
                return Hodos
        return []
    def near(self,x,mean,sg,nsig):
        return self.btwn(x,mean-nsig*sg,mean+nsig*sg)
    def btwn(self,x,xlo,xhi):
        return xlo<=x and x<=xhi
    def makeList(self,q):
        if isinstance(q, list) : return q
        return [q]
    def getTDC(self,data):
        TDC = {}
        for x in data:
            TDC[x] = data[x]
        return TDC
    def getQDC(self,data):
        QDC = {}
        for x in data:
            QDC[x] = data[x]
        return QDC
    def getWFD(self,data):
        WFDtime,WFDarea,WFDnp,WFDped,WFDpsd = {},{},{},{},{}
        #print 'second.getWFD data',data,data.name
        #self.writer.show('xx',group=data.name,h5pyf=self.f)
        for x in data['time']:
            WFDtime[x] = data['time'][x]
            WFDarea[x] = data['area'][x]
            WFDnp[x]   = data['npulse'][x]
            WFDped[x]  = data['ped'][x]
            WFDpsd[x]  = data['pedsd'][x]
        return WFDtime,WFDarea,WFDnp,WFDped,WFDpsd
    def endroot(self,rfn):
        '''
        make summary plots
        write out results to root file
        delete all hists/trees
        '''
        p1list = ['CT_dT','CT_dT','CT_TDC','WFD_At_vs_Traj','WFD_At_vs_time','npulse','pedestal','ped_stdev']
        p2list = ['wrtH0','wrtH2',''      ,''              ,''              ,''      ,''        ,'']
        for p1 in ['WFD_time','WFD_dt','WFD_area','WFD_Area_tcut']:
            for t in ['CT','M','LED']:
                p1list.append( t + '_' + p1)
                p2list.append('')

        for phrase1,phrase2 in zip( p1list,p2list ):
            hlist = []
            for h in self.Hists:
                if (phrase1) in h and (phrase2 in h): hlist.append(self.Hists[h])
            if len(hlist)>0:
                hlist.sort()
                fname = phrase1
                if len(phrase2)>0: fname += '_'+phrase2
                abscissaIsTime = 'vs_time' in phrase1.lower()
                dopt = ''
                Logy = True
                sopt = 1111111
                if '_vs_' in phrase1.lower():
                    dopt = 'box'
                    Logy = False
                    sopt = 0
                self.gU.drawMultiHists(hlist,fname=fname,figdir=self.parentDir,setLogy=Logy,abscissaIsTime=abscissaIsTime,dopt=dopt,statOpt=sopt)
                if 0 and 'vs_time' in phrase1.lower():
                    for x in hlist:  self.gU.reportHist(x)  #### FOR DEBUG

        
        rf = TFile(rfn,"RECREATE")
        for h in self.Hists:rf.WriteTObject(self.Hists[h])
        nwrite = len(self.Hists)
        rf.Close()
        print 'second.endroot closed',rfn

        # next line deletes all hists/trees from memory according to Rene Brun 3Apr2002
        gDirectory.GetList().Delete()
        print 'second.endroot Wrote',nwrite,'objects to',rfn,'and deleted them from memory'

        
        return
if __name__ == '__main__' :
    maxevt = 9999999
    if len(sys.argv)>1: maxevt = int(sys.argv[1])

    fn = 'Results/20160208_163935/Output/run600-run602.h5'

    S = second()
    datadir = '/Users/djaffe/work/GIT/ONETON/ReconCalibDataFiles/'
    fnlist = S.P.getRawDataList(datadir)

    #fnlist = fnlist[:1] ##### TEMPORARY

    #fnlist = [ 'ReconCalibDataFiles/run600-run602.h5.gz' ] #### TEST GZIP
    #fnlist = [ 'ForTesting/run600-run602.h5' ]  ##### TEST SOME FIGURES


    fnlist.sort() # order filelist by runs
    
    rfn = S.parentDir + 'second.root'
    print 'second__main__',len(fnlist),'files in input list'


    S.book()
    for fn in fnlist:
        ok = S.open(fn)
        if ok: 
            S.loop(maxevt=maxevt)
            S.close()
    S.endroot(rfn)
