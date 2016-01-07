#!/usr/bin/env python
'''
read in and parse hdf5 file
20151231
'''
import h5py
import numpy
import graphUtils
from ROOT import TFile


class reader():
    def __init__(self):
        self.datasetNames = [ "Digitizer_1","Digitizer_2", "Event_Temp", "Event_Time", "QDC_1","QDC_2", "Scaler", "TDC"]
        self.scalerTypes =  ['float32','float64','uint32']

        # used in unpackConfigurationSettings
        self.iDNCAT = None
        self.hiQDCbinwidth = None
        self.hiQDCthresbin = None
        self.loQDCbinwidth = None
        self.loQDCthresbin = None
        self.TDCbinwidth   = None
        self.TDCthresbin   = None

        
        self.gU = graphUtils.graphUtils()
        print 'reader: initialized'
        return
    def first(self,fn=None):
        print '\n reader.first filename is',fn
        f = h5py.File(fn)
        for x in f: print x,
        print ''

        RunInfo = f['Run_Info']
        CalData = f['Calibration']
        EvtDir = f['Events']

        print 'RunInfo:'
        for x in RunInfo:
            print x,RunInfo[x]
            for y in RunInfo[x]:
                print y,RunInfo[x][y],RunInfo[x][y][()]
        
        self.unpackRunInfo(RunInfo)

        
        #print EvtDir
        Hists = []
        for evtkey in EvtDir.keys():
            print 'first event'
            print evtkey,EvtDir[evtkey],EvtDir[evtkey].name

            for thing in self.datasetNames: #EvtDir[evtkey]:
                d = EvtDir[evtkey][thing].dtype
                print thing,EvtDir[evtkey][thing],'shape:',EvtDir[evtkey][thing].shape,'dtype:',EvtDir[evtkey][thing].dtype,EvtDir[evtkey][thing].name
                if d in self.scalerTypes:
                    print thing,EvtDir[evtkey][thing][()]
                else:
                    for x in EvtDir[evtkey][thing]:
                        print x,
                    print ''
                if 'Scaler' in thing:
                    scaler = self.unpackScaler(EvtDir[evtkey][thing])

                    print ''
                        
                if 'Digitizer' in thing:
                    WFD = self.unpackDigitizer(EvtDir[evtkey][thing])
                    WFDps=self.unpackDigitizer(EvtDir[evtkey][thing],CalData[thing]['Pedestal'])
                    #print thing,'WFD',WFD
                    #print 'ped-subtracted',WFDps
                    H = self.displayWFD(WFD,pretitle=thing+'raw')
                    Hists.extend(H)
                    Hps = self.displayWFD(WFDps,pretitle=thing+'pedsub')
                    Hists.extend(Hps)
                    self.gU.drawMultiHists(H,fname=thing+'raw',figdir='Figures/WF/',dopt='HIST',statOpt=0)
                    self.gU.drawMultiHists(Hps,fname=thing+'pedsub',figdir='Figures/WF/',dopt='HIST',statOpt=0)

            break
        rfn = 'rdr.root'
        rf = TFile(rfn,"RECREATE")
        for h in Hists: rf.WriteTObject(h)
        rf.Close()
        print 'wrote',len(Hists),'hists to',rfn
        return
    def unpackRunInfo(self,raw):
        print 'Elements in RunInfo'
        for elm in raw:
            print '\n',elm
            if elm=='Configuration_Settings':
                self.unpackConfigurationSettings(raw[elm])
            elif elm=='Digitizer_Data':
                for x in raw[elm]: print x
                print ''
                self.unpackDigitizerData(raw[elm])
            elif elm=='QDC_Data':
                print ''
            elif elm== 'Run_Details':
                print ''
            elif elm== 'Scaler_Data':
                print ''
            elif elm== 'TDC_Data':
                print ''
            elif elm== 'Temperature':
                print ''
            elif elm== 'Time_Info':
                print ''
            elif elm== 'Trigger_Bits':
                print ''
            elif elm== 'VME_Configuration':
                print ''
            else:
                w = 'reader.unpackRunInfo: WARNING Unexpected element',elm
                print w
        return
    def unpackDigitizerData(self,raw):
        self.wfdChanDesc = [x for x in raw['Channel_Description']]
        self.wfdOHodo = raw['Only_Hodo']
        print 'wfdChanDesc:',self.wfdChanDesc

        return

    def unpackConfigurationSettings(self,raw):
        self.iDNCAT = raw['Digitizer_Num_Cols_After_Trig']
        self.hiQDCbinwidth = raw['QDC_High_Range_Bin_Width_fC']
        self.hiQDCthresbin = raw['QDC_High_Range_Threshold_bin']
        self.loQDCbinwidth = raw['QDC_Low_Range_Bin_Width_fC']
        self.loQDCthresbin = raw['QDC_Low_Range_Threshold_bin']
        self.TDCbinwidth   = raw['TDC_Bin_Width_ns']
        self.TDCthresbin   = raw['TDC_Threshold_bin']
        return
    def displayWFD(self,WFD,pretitle=''):
        '''
        return list of histograms, one per WFD channel
        use pretitle to make hist name, title
        '''
        H = []
        for chan in WFD:
            bins = numpy.array(range(len(WFD[chan])))
            title = pretitle + ' ' + str(chan)
            h = self.gU.makeTH1Dwtd(bins,WFD[chan],title)
            H.append(h)
        return H
    def unpackDigitizer(self,raw,ped=None):
        '''
        unpack 4 channel digitizer input
        return dict[chan] = [v0,v1,...]
        subtract pedestals if provided
        '''
        d = {}
        data = raw
        if ped is not None:
            data = numpy.subtract(raw,ped)
        for i,x in enumerate(data):
            d[i] = numpy.array(x)
        return d
    def unpackScaler(self,raw):
        '''
        unpack sl
        '''
        s = numpy.array(raw)
        return s
if __name__ == '__main__' :
    r = reader()
    fn = '/Users/djaffe/work/1TonData/run462.h5'
    r.first(fn=fn)
