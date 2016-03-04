#!/usr/bin/env python
'''
read in and parse hdf5 file
20151231
'''
import h5py
import numpy
import graphUtils
from ROOT import TFile
import sys
import datetime
import os,zipfile


class reader():
    def __init__(self):
        self.scalerTypes =  ['float32','float64','uint32']

        # used by unpackEvent, summary
        self.datasetNames = [ "Digitizer_1","Digitizer_2", "Event_Temp", "Event_Time", "QDC_1","QDC_2", "Scaler", "TDC"]
        
        self.QDCNames = [x for x in self.datasetNames if 'QDC' in x]
        self.TDCNames = [x for x in self.datasetNames if 'TDC' in x]
        self.WFDNames = [x for x in self.datasetNames if 'Digitizer' in x]
        #print 'reader.init QDCnames',self.QDCNames,'TDCNames',self.TDCNames,'WFDNames',self.WFDNames
        
        self.eventCalls = 0
        self.incompleteEvents = {}

        # used by unpackTrigger
        self.validTriggers = ['CT','M','LED','CT+M+LED'] # as of 20160106
        self.triggerOR     = ['CT+M+LED']
        self.uniqueTriggers = [x for x in self.validTriggers if x not in self.triggerOR]
        #print 'reader__init__ self.uniqueTriggers',self.uniqueTriggers
        
        self.gU = graphUtils.graphUtils()
        print 'reader: initialized'
        return
    def start(self,fn=None,RequireCalib=True):
        '''
        start of run initialization
        handle zip files using zipfile:
        identify as zipfile, check if namelist in zipped files agrees with expected name,
        then extract .h5 file. keep track of name of extracted file so file can be deleted
        after use
        If RequireCalib, then report it, go through file closing procedure and return False
        
        Return True for successful open satisfying all requirements.
        '''
        OK = True
        bn = os.path.basename(fn)
        bnsuf = bn.split('.')[-1]
        if bnsuf=='zip':
            print 'reader.start open,extract',fn
            zf = zipfile.ZipFile(fn)
            namel = zf.namelist()
            bnpre = bn.split('.')[0]
            if namel[0]!=bnpre+'.h5': # should not happen
                sys.exit('reader.start ERROR processing ' + fn + ' because namelist ' + namelist[0] + ' does not match prefix ' + bnpre)
            h5f = zf.extract(namel[0])
            self.localHDF5file = namel[0]
            self.f = h5py.File(h5f,'r')
        elif bnsuf=='h5':
            self.f = h5py.File(fn,'r')
            self.localHDF5file = None
        else:
            sys.exit('reader.start ERROR processing ' + fn +  ' UNKNOWN suffix ' + bnsuf)
            
        print 'reader.start run. File',fn
        
        self.RunInfo = self.f['Run_Info']
        self.RunNumber = self.getRunDetail('run')
        if 'Calibration' in self.f:
            self.CalData = self.f['Calibration']
        else:
            self.CalData = None
            print 'XXXXXXXXXXXXXX NO CALIBRATION DATA IN FILE XXXXXXXXXXX'
            if RequireCalib: OK = False
        self.EvtDir  = self.f['Events']

        if not OK: self.closeHDF5File()

        return OK
    def closeHDF5File(self):
        '''
        close the current HDF5
        check to see if a local file hdf5 was made by unzipping
        '''
        name = self.f.filename
        self.f.close()
        words = 'reader.closeHDF5File closed ' + name
        if self.localHDF5file is not None:
            if os.path.isfile(self.localHDF5file):
                os.remove(self.localHDF5file)
                words += ' and deleted it'
        print words
        return
    def summary(self):
        '''
        summary of stats accumulated by reader
        '''
        L = len(self.incompleteEvents)
        print ' ++++++ reader.summary',self.eventCalls,'events with',L,'events with a missing dataset'
        if L>0:
            print ' ++++++ reader.summary missing datasets',
            for x in self.incompleteEvents: print self.incompleteEvents[x],x,',',
            print ''
        return
    def getEvtDir(self):
        return self.EvtDir
    def getCalData(self):
        return self.CalData
    def testGet(self):
        '''
        exercise getXXX modules for RunInfo
        '''
        print 'reader.testGet....................................'
        for Q in ['Run_Number','Run_Type','Material','Comments']: print 'RunDetail',Q,self.getRunDetail(Q)
        print 'TrigBits',self.getTriggerBits('anything')
        print 'Temperature',self.getTemperature()
        for Q in ['Start_Time_str','Start_Time_UNIX','datetime']: print 'Time',Q,self.getTime(Q)
        for Q in ['Model_Name','Module_Description','VME_Address']: print 'VMEConfig',Q,self.getVMEConfig(Q)
        print 'Digitizer_Num_Cols_After_Trig',self.getDigitizerNumColAfterTrig()
        for Q in ['HighWidth','LowWidth','HighThres','LowThres']:
            print 'QDC'+Q,self.getQDCConfig(Q)
        for Q in ['width','thres']:
            print 'TDC'+Q,self.getTDCConfig(Q)
        for M in ['Digitizer_Data','TDC_Data','QDC_Data','Scaler_Data']:
            for Q in ['ChanDesc','Hodo']: print 'ModuleData',M,Q,self.getModuleData(M,Q)
        print '............................................reader.testGet'
        return
    def getTriggerBits(self,Q):
        '''
        from RunInfo return definition of trigger bits given input Q
        there is only one valid input Q
        '''
        TB = self.RunInfo['Trigger_Bits']
        if Q in TB: return [x for x in TB[Q]]
        return [x for x in TB['Description_of_Trigger_Bits']]
    def getVMEConfig(self,Q):
        '''
        from RunInfo return config of VME slots model, module descrip, VME address
        '''
        VC = self.RunInfo['VME_Configuration']
        if Q in VC: return [x for x in VC[Q]]
        sys.exit('reader.getVMEConfig ERROR Invalid input '+Q)
        return
    def getTime(self,Q):
        '''
        from RunInfo return start time as string, unix time or as datetime object
        '''
        TI = self.RunInfo['Time_Info']
        if Q in TI: return TI[Q][()]
        if Q.lower()=='datetime':
            S = TI['Start_Time_str'][()] # unused for now
            U = TI['Start_Time_UNIX'][()]
            return datetime.datetime.fromtimestamp(U)
        sys.exit('reader.getTime ERROR Invalid input '+Q)
        return
    def getTemperature(self):
        '''
        from RunInfo return position of temperature probes
        '''
        return self.RunInfo['Temperature']['Location'][()]
    def getRunDetail(self,Q):
        '''
        from RunInfo return details of run
        run number
        run type
        material in vessel
        comments
        '''
        RD = self.RunInfo['Run_Details']
        if Q in RD          : return RD[Q][()]
        if Q.lower()=='run' : return RD['Run_Number'][()]
        if Q.lower()=='type': return RD['Run_Type'][()]
        sys.exit('reader.getRunDetail ERROR Invalid input '+Q)
        return
    def getRunNum(self):
        return self.getRunDetail('run')
    def reportRunDetail(self):
        '''
        summarize run details
        '''
        runnum = self.getRunDetail('run')
        runtype= self.getRunDetail('type')
        mat    = self.getRunDetail('Material')
        st     = self.getTime('Start_Time_str')
        c      = self.getRunDetail('Comments')
        print runtype,'run',runnum,'Start',st,'Fill',mat
        print c
        return
    def getModuleData(self,module,Q):
        '''
        from RunInfo return per channel info for QDC, TDC, Scaler
        '''
        m = module 
        if module not in self.RunInfo:
            m = module + '_Data'
            if m not in self.RunInfo:
                sys.exit('self.getModuleData ERROR Invalid module '+module)
        DD = self.RunInfo[m]
        if Q in DD      : return [x for x in DD[Q]]
        if Q=='ChanDesc': return [x for x in DD['Channel_Description']]
        if Q=='Hodo'    : return [x for x in DD['Only_Hodo']]
        sys.exit('reader.getModuleData ERROR Module Invalid input quantity '+Q)
        return
    def getDigitizerNumColAfterTrig(self):
        return self.RunInfo['Configuration_Settings']['Digitizer_Num_Cols_After_Trig'][()]
    def getQDCConfig(self,Q):
        '''
        from RunInfo return QDC configuration setting given input string Q
        high or low range bin width in fC
        high or low range threshold bin
        fail if Q is invalid
        '''
        CS = self.RunInfo['Configuration_Settings']
        if Q in CS: return CS[Q][()]
        if Q=='HighWidth': return CS['QDC_High_Range_Bin_Width_fC'][()]
        if Q=='LowWidth' : return CS['QDC_Low_Range_Bin_Width_fC'][()]
        if Q=='HighThres': return CS['QDC_High_Range_Threshold_bin'][()]
        if Q=='LowThres' : return CS['QDC_Low_Range_Threshold_bin'][()]
        sys.exit('reader.getQDCConfig ERROR Invalid input quantity ' + Q)
        return 
    def getTDCConfig(self,Q):
        '''
        from RunInfo return TDC configuration setting given input string Q
        TDC bin width in ns
        TDC threshold bin
        '''
        CS = self.RunInfo['Configuration_Settings']
        if Q in CS: return CS[Q][()]
        if Q=='width' or Q=='binwidth': return CS['TDC_Bin_Width_ns'][()]
        if Q=='thres' : return CS['TDC_Threshold_bin'][()]
        sys.exit('reader.getTDCConfig ERROR Invalid input quantity ' + Q)
        return
    def displayWFDnew(self,WFD,detel,pretitle=None,NX=None,XMI=None,XMA=None):
        '''
        return a single histogram of input waveform
        hist name and title derived from pretitle and detel = detector element
        '''
        title = detel
        if pretitle is not None: title = pretitle + detel
        #print 'reader.displayWFDnew detel,WFD=',detel,WFD
        bins = numpy.array(range(len(WFD)))
        h = self.gU.makeTH1Dwtd(bins,WFD,title,NX=NX,XMI=XMI,XMA=XMA)
        return h
    def displayWFD(self,WFD,pretitle='',ipt=None):
        '''
        return list of histograms, one per WFD channel
        use pretitle to make hist name, title
        OR use ipt as title of each channel
        '''
        H = []
        pret = ipt
        if ipt is None:
            pret = []
            for chan in WFD: pret.append(pretitle + str(chan))
        elif len(ipt)<len(WFD):
            pret = []
            for chan in WFD: pret.append(pretitle + str(chan))
                

        #print 'reader.displayWFD len(ipt),len(WFD),pret',len(ipt),len(WFD),pret
        
        for p,chan in zip(pret,WFD):
            bins = numpy.array(range(len(WFD[chan])))
            title = p
            h = self.gU.makeTH1Dwtd(bins,WFD[chan],title)
            H.append(h)
        return H
    def addChanDesc(self,Unpack,module):
        '''
        given list of unpacked data for TDC, QDC or Digitizer and the module
        return ordered list of channel description

        Assumes single TDC, two QDCs, two WFDs
        '''
        CD = []

        s  = 'ChanDesc'
        if module in self.TDCNames:
            names = self.getModuleData(module,s)
            for b in Unpack:
                i = int(b[0])
                CD.append( names[i] )
        elif module in self.QDCNames:
            names = self.getModuleData('QDC',s)
            #print 'reader.addChanDesc:names',names,'module',module,'Unpack',Unpack
            i0 = 0
            if '2' in module: i0 = 8
            for b in Unpack:
                i = i0 + int(b[0]) 
                #print 'reader.addChanDesc:b,i0,i',b,i0,i
                CD.append( names[i] )
        elif module in self.WFDNames:
            names = self.getModuleData('Digitizer',s)
            i0 = 0
            if '2' in module: i0 = 4
            for j,b in enumerate(Unpack):
                i = i0 + j
                CD.append( names[i] )
            pass
        else:
            w = 'reader.addChanDesc ERROR Invalid module '+module
            sys.exit(w)
        return CD
    def unpackEvent(self,raw):
        '''
        check that all datasets are present in event
        return list of missing datasets and tally it
        '''
        self.eventCalls += 1
        missing = []
        for ds in self.datasetNames:
            if ds in raw:
                continue
            else:
                missing.append(ds)
                if ds not in self.incompleteEvents: self.incompleteEvents[ds] = 0.
                self.incompleteEvents[ds] += 1
        return missing
    def assocWFD(self,raw,pedsub=False):
        '''
        input is signal event, pedsub=True, then subtract pedestals, otherwise no.
        generate error if pedestals are not available
        return dict with key=detector element, value = array of WF data
        '''
        if pedsub and (self.CalData is None):
            sys.exit('\n reader.assocWFD ERROR Cannot subtract pedestals. No pedestal data available \n')
            
        U,WFD = {},{}
        for Q in self.WFDNames:
            if pedsub:
                U[Q] = self.unpackDigitizer(raw[Q],ped=self.CalData[Q]['Pedestal'])
            else:
                U[Q] = self.unpackDigitizer(raw[Q])
            cd  = self.addChanDesc(U[Q],Q)
            #print 'readder.assocWFD Q,U[Q],cd',Q,U[Q],cd
            for x,y in zip(cd,U[Q]): WFD[x] = U[Q][y]
        #print 'reader.assocWFD U=',U
        #print 'reader.assocWFD WFD=',WFD
        return WFD
    def unpackDigitizer(self,raw,ped=None):
        '''
        unpack 4 channel digitizer input
        return dict[chan] = [v0,v1,...]
        subtract pedestals if provided
        '''
        d = {}
        data = raw
        if ped is not None:
            if 493<self.RunNumber and self.RunNumber<751 :
                ped = self.fixPed(ped)
            data = numpy.subtract(raw,ped)
        for i,x in enumerate(data):
            d[i] = numpy.array(x)
        return d
    def fixPed(self,ped):
        '''
        Lindsey email of 19jan2016:
        OK, I've now fixed this bug.
    It was caused by overflowing an undersized variable.
    I got the error because I changed the number of averaged waveforms from 10 to 100.
    It is pretty easy to recover the true pedestal data from the values that were written:
    Pedestal_true = ((Pedestal_written)x100 + (3x65536 - 1))/100
    The pedestal data from run 751 onwards will not require that correction.
        '''
        offset = float((3*65536 - 1))
        q = numpy.multiply(ped,100.)
        r = numpy.add(q,offset)
        s = numpy.divide(r,100.)
        return s
        
    def assocTDC(self,raw,mode=None):
        '''
        input is a single event
        return dict with key=det element, value = TDC value
        mode is defined in unpackTDC
        '''
        TDC,U = {},{}
        for Q in self.TDCNames:
            U[Q] = self.unpackTDC(raw[Q],mode=mode)
            cd  = self.addChanDesc(U[Q],Q)
            for x,y in zip(cd,U[Q]): TDC[x] = y
        return TDC
    def unpackTDC(self,raw,mode=None):
        '''
        return array of pairs (channel, TDC(channel) ) given mode
        mode = None, 'raw' => raw data
        mode = 'inns','ns' => convert raw data to ns
        '''
        if mode is None or mode=='raw':
            TDC = numpy.array(raw)
        elif mode=='inns' or mode=='ns':
            c = self.getTDCConfig('width')
            TDC = []
            for pair in raw:
                TDC.append( [pair[0], c*float(pair[1])] )
        else:
            sys.exit('reader.unpackTDC ERROR Invalid mode '+mode)
        return TDC
    def assocQDC(self,raw,removeNC=True):
        '''
        input is a single event
        if removeNC is True, then output dict will not have detector elements 'N/C'
        return dict with key=detector element (i.e., 'S0'), value = list of QDC elements
        each QDC element is an array (chan#, value, lo/hi, overflow?)
        
        '''
        U = {}
        QDC = {}
        for Q in self.QDCNames:
            U[Q] = self.unpackQDC(raw[Q])
            cd  = self.addChanDesc(U[Q],Q)
            for x,y in zip(cd,U[Q]):
                if removeNC and x=='N/C':
                    pass
                else:
                    if x not in QDC: QDC[x] = []
                    QDC[x].append(y)
        return QDC

    def unpackQDC(self,raw):
        s = numpy.array(raw)
        return s
    def unpackScaler(self,raw):
        '''
        unpack sl
        '''
        s = numpy.array(raw)
        return s
    def unpackTrigger(self,raw):
        '''
        return list of triggers that fired based on TDC information
        excludin the OR of all triggers
        '''
        TDCmd = self.getModuleData('TDC','ChanDesc')
        trigs = []
        for pair in raw:
            #print 'reader.unpackTrigger pair',pair
            T = TDCmd[int(pair[0])]
            if  T in self.validTriggers:
                trigs.append(T)
        # special treatment since LED signal is outside TDC dynamic range
        if trigs==self.triggerOR : trigs.append('LED')
        # exclude OR of all triggers
        trigs.remove(self.triggerOR[0])
        return trigs
    def unpackTemperature(self,raw):
        return raw['Event_Temp'][()]
    def unpackTime(self,raw):
        return raw['Event_Time'][()]
    def printDataset(self,X):
        '''
        avoid annoying "TypeError: Can't iterate over a scalar dataset"
        when trying to print
        '''
        if X.size==1.:
            more = False
            print X[()],
        else:
            more = True
            print X,

        print ''
        return more
if __name__ == '__main__' :
    #r = reader()
    #fn = '/Users/djaffe/work/1TonData/run462.h5'
    #r.first(fn=fn)
    print("directly calling data file reader object")