#!/usr/bin/env python
'''
write out hdf5 file for 1ton
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
    def writeEvent(self,data,datalabel,evtNo=None):
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
        for l,d in zip(datalabel,data):
            self.r.create_dataset(Elab+l,data=d)
        return eN
    def printname(self,name):
        print name
        return
    def testrd(self,fn):
        self.f = f = h5py.File(fn,'r')
        
        print 'testing',fn
        for a in f:
            print a,
            for b in f[a]:
                print b,
                for c in f[a][b]:
                    print c,
                    for d in f[a][b][c]:
                        print d,
                        for e in f[a][b][c][d]:
                            print e,f[a][b][c][d][e][()],
                print ''
        f.close()
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
            evtn = w.writeEvent(data,datalabel)
            q = []
            for i in range(3): q.append(float(ev+1)*random.random())
            data = [ numpy.array(q) ]
            evtn = w.writeEvent(data,dl2,evtNo=evtn)
            
    w.closeFile()
    w.testrd(fn)
        
        
            
    
