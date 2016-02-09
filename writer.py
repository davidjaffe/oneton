#!/usr/bin/env python
'''
write out hdf5 file for 1ton

workflow
loop on runs
   open input file for each run
   openFile for output
   setRunNum
   loop on events
      writeEvent
   close input file
close output file
   
20160203
'''
import h5py
import numpy

import random # for testing 

class writer():
    def __init__(self):
        self.f = None
        self.r = None
        self.RunNum   = None
        self.EventNum = None
        return
    def openFile(self,fn):
        '''
        open hdf5 file for writing
        '''
        self.f = h5py.File(fn,'w')
        print 'writer.openFile opened',self.f.filename
        return
    def closeFile(self):
        '''
        close currently open file
        '''
        n = self.f.filename
        self.f.close()
        print 'writer.closeFile closed',n
        return
    def setRunNum(self,rn):
        '''
        set current run number
        reset event number
        '''
        self.RunNum = rn
        rnm = str(self.RunNum).zfill(6)
        self.r = self.f.create_group('Run/'+rnm)
        self.EventNum = None
        return
    def writeEvent(self,datalabel,data,evtNo=None):
        '''
        write array of data to hdf5 event group
        can be called multiple times for same event
        if no event number given, then event# is incremented
        returns current event number
        '''
        if evtNo is None:
            if self.EventNum is None:
                self.EventNum = 0
            else:
                self.EventNum += 1
        else:
            self.EventNum = evtNo
        eN = self.EventNum
        enm = str(eN).zfill(6)
        #E = self.r.create_group('Event/'+enm)
        Elab = 'Event/'+enm+'/'
        if type(datalabel) is list:
            DL,D = datalabel,data
        else:
            DL,D = [datalabel],[data]
        for l,d in zip(DL,D):
            #print 'writer.writeEvent Elab+l,d',Elab+l,d
            self.r.create_dataset(Elab+l,data=d)
            
        return eN
    def writeRunData(self,label,data):
        '''
        write run data
        add data to output with path given by label + path for this run
        data and label should be equal length lists
        '''
        if type(label) is list:
            DL,D = label,data
        else:
            DL,D = [label],[data]
        for l,d in zip(DL,D):
            self.r.create_dataset(l,data=d)
            print 'writer.writeRunData label',l,'data',d
        return
    def writeData(self,label,data):
        '''
        add data to output with path given by label.
        data and label should be equal length lists
        '''
        if type(label) is list:
            DL,D = label,data
        else:
            DL,D = [label],[data]
        for l,d in zip(DL,D):
            self.f.create_dataset(l,data=d)
            print 'writer.writeData label',l,'data',d
        return

    def find(self,name):
        if name is not None: self.Found.append(name)
        return None
    def uniqify(self,L):
        '''
        return list with only unique strings
        example: L = ['a','aa','ab','abc'] => ['aa','abc']
        '''
        more = True
        K = L
        while more:
            more = False
            takeOut = []
            for i,a in enumerate(K):
                for b in K[i+1:]:
                    if a in b:
                        takeOut.append(a)
                        break 
            for b in takeOut:
                more = True
                i = K.index(b)
                del( K[i] )
        return K
    def show(self,fn,group=None,h5pyf=None):
        '''
        print content of uniqified elements of group in hdf5 file named fn
        OR contents of elements in h5pyf
        depends on writer.find, write.uniqify
        '''
        if h5pyf is None:
            f = h5py.File(fn,'r')
        else:
            f = h5pyf
        self.Found = []
        if group is None:
            f.visit(self.find)
        elif group in f:
            f[group].visit(self.find)
            self.Found = ['{1}/{0}'.format(a,group) for a in self.Found]

        #print 'BEFORE self.Found',self.Found
        self.Found = self.uniqify(self.Found)
        #print 'AFTER self.Found',self.Found
        
        print '\nwrite.show group',group,'in file',fn
        if len(self.Found)==0:print 'write.show Nothing found for group',group
        for a in self.Found:
            print a,
            if f[a].shape==():
                print f[a][()],
            else:
                for b in f[a]: print b,
            print ''
        if h5pyf is None: f.close()
        return

if __name__ == '__main__' :
    w = writer()
    fn = 'bugshdf5.h5'
    w.openFile(fn)

    datalabel = ['time','temperature']
    dl2 = ['banana']
    for rn in range(387,389):
        w.setRunNum(rn)
        for ev in range(4):
            faketime = float(ev)+random.random()+0.5
            faketemp = random.gauss(float(ev),1.)
            data = [faketime, faketemp]
            evtn = w.writeEvent(datalabel, data)
            q = []
            for i in range(3): q.append(float(ev+1)*random.random())
            data = [ numpy.array(q) ]
            evtn = w.writeEvent(dl2,data,evtNo=evtn)

    datalabel = ['/Calibration/TDC', '/Calibration/ADC']
    data = [ numpy.array([random.random() for x in range(3)]), numpy.array([random.gauss(0.,.1) for x in range(3)]) ]
    w.writeData(datalabel,data)
    
    w.closeFile()

    w.show(fn)
    w.show(fn,'/Run/000387')
    w.show(fn,'/Run/000222')
    w.show(fn,'/Calibration')
        
        
            
    
