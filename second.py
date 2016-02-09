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
    def loop(self):
        for run in self.f['Run']:
            runnum = int(run)
            print 'second.loop run',runnum
            Events = self.f['Run/'+run+'/Event']
            for event in Events:
                evtnum = int(event)
                if evtnum%1000==0 : print 'second.loop event',evtnum
                data = Events[event]
                TDC = self.getTDC(data['TDC'])
                QDC = self.getQDC(data['QDC'])
                WFDtime = WFDarea = None
                if 'WFD' in data: WFDtime,WFDarea = self.getWFD(data['WFD'])
                Time = data['Time']
                Temp = data['Temp']
        return
    def book(self):
        self.Hists = {}
        for x in ['S0','S1','S2','S3','S4','S5']:
            name = 'WFD_dt_'+x
            title = name.replace('_',' ')
            self.Hists[name] = TH1D(name,title,100,0.,400.)
        print 'second.book booked hists'
        return
    
            
    def analyze(self,TDC,QDC,WFDtime,WFDarea,Time,Temp):
        '''
        analysis
        '''
        if WFDtime is not None:
            for x in WFDtime:
                ts = WFDtime[x]
                name = 'WFD_dt_'+x
                if len(ts)>1:
                    for i,t in enumerate(ts):
                        if i+1<len(ts):
                            dt = ts[i+1]-t
                            self.Hists[name].Fill(dt)
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
        rf.close()
        print 'second.endroot closed',rfn
        return
if __name__ == '__main__' :
    fn = 'Results/20160208_163935/Output/run600-run602.h5'
    rfn = 'second.root'
    S = second()
    S.book()
    S.open(fn)
    S.loop()
    S.close()
    S.endroot(rfn)
