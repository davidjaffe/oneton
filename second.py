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
import reader
import math
import os
import Logger
import wfanal
import calibThermocouple
import writer
from optparse import OptionParser

class second():
    def __init__(self):
        self.f = None
        self.writer = writer.writer()

        self.CTtimes = {} # for WFD
        for c in ['S0','S1','S2','S3','S4','S5','S6','S7']:
            self.CTtimes[c] = [120.,200.]
            if c in ['S4','S5']: self.CTtimes[c] = [90.,160.]
        
        return
    def open(self,fn):
        self.f = h5py.File(fn,'r')
        print 'second.open Opened',self.f.filename
        return
    def close(self):
        n = self.f.filename
        self.f.close()
        print 'second.close Closed',n
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
                if evtnum%1000==0 : print 'second.loop event',evtnum
                if evtnum>maxevt: break
                data = Events[event]
                TDC = self.getTDC(data['TDC'])
                QDC = self.getQDC(data['QDC'])
                WFDtime = WFDarea = None
                if 'WFD' in data: WFDtime,WFDarea = self.getWFD(data['WFD'])
                Time = data['Time']
                Temp = data['Temp']
                self.analyze(TDC,QDC,WFDtime,WFDarea,Time,Temp)
        return
    def book(self):
        self.Hists = {}
        for trig in ['CT','M','LED']:
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

                name,nx,xmi,xma = trig + '_WFD_area_tcut_'+x,200,0.,2000.
                title = name.replace('_',' ')
                self.Hists[name] = TH1D(name,title,nx,xmi,xma)
            if 'LED'!=trig:
                Href = 'H0'
                for H in ['H0','H1','H2','H3','H4','H5']:
                    if H!=Href:
                        name,nx,xmi,xma = trig+'_dT_'+H+'wrt'+Href,100,-20.,20.
                        title = name.replace('_',' ')
                        self.Hists[name] = TH1D(name,title,nx,xmi,xma)

                    name,nx,xmi,xma = trig+'_TDC_'+H,100,0.,200.
                    title = name.replace('_',' ')
                    self.Hists[name] = TH1D(name,title,nx,xmi,xma)
                
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
            
    def analyze(self,TDC,QDC,WFDtime,WFDarea,Time,Temp):
        '''
        analysis
        '''
        triggers = self.getTriggers(TDC)
        if WFDtime is not None:
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
            if 'CT' in triggers:
                for x in WFDtime:
                    for ts,areas in zip(WFDtime[x],WFDarea[x]):
                        for t,a in zip(self.makeList(ts),self.makeList(areas)):
                            if self.CTtimes[x][0]<=t and t<=self.CTtimes[x][1]:
                                name = trig + '_WFD_area_tcut_'+x
                                self.Hists[name].Fill(abs(a))
        for trig in triggers:
            if 'LED'!=trig:
                Href = 'H0'
                tref = None
                if Href in TDC: tref = TDC[Href][1]
                for H in ['H0','H1','H2','H3','H4','H5']:
                    if H in TDC:
                        t = TDC[H][1]
                        #print 'H,t',H,t
                        name = trig+'_TDC_'+H
                        self.Hists[name].Fill(t)
                        if H!=Href and tref is not None:
                            name = trig+'_dT_'+H+'wrt'+Href
                            self.Hists[name].Fill(t-tref)
        return
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
        WFDtime,WFDarea = {},{}
        #print 'second.getWFD data',data,data.name
        #self.writer.show('xx',group=data.name,h5pyf=self.f)
        for x in data['time']:
            WFDtime[x] = data['time'][x]
            WFDarea[x] = data['area'][x]
        return WFDtime,WFDarea
    def endroot(self,rfn):
        rf = TFile(rfn,"RECREATE")
        for h in self.Hists:rf.WriteTObject(self.Hists[h])
        rf.Close()
        print 'second.endroot closed',rfn
        return
if __name__ == '__main__' :
    maxevt = 9999999
    if len(sys.argv)>1: maxevt = int(sys.argv[1])

    fn = 'Results/20160208_163935/Output/run600-run602.h5'
    rfn = 'second.root'
    S = second()
    S.book()
    S.open(fn)
    S.loop(maxevt=maxevt)
    S.close()
    S.endroot(rfn)
