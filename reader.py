#!/usr/bin/env python
'''
read in and parse hdf5 file
20151231
'''
import h5py


class reader():
    def __init__(self):
        return
    def first(self,fn=None):
        f = h5py.File(fn)
        EvtDir = f['Events']
        for evtkey in EvtDir.keys():
            print 'first event'
            print EvtDir[evtkey]
            for thing in EvtDir[evtkey]:
                print thing,EvtDir[evtkey][thing],'shape:',EvtDir[evtkey][thing].shape,'dtype:',EvtDir[evtkey][thing].dtype
    #            for x in EvtDir[evtkey][thing]:
    #                print x,
                print ''
            break
        return
if __name__ == '__main__' :
    r = reader()
    fn = '/Users/djaffe/work/1TonData/run273.h5'
    r.first(fn=fn)
