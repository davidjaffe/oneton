#!/usr/bin/env python
'''
simple estimate of precision of cosmic telescope
units are cm
20160301
'''
import numpy
#import graphUtils
#import ROOT
#from ROOT import TFile,TH1D,TH2D,gDirectory
import math

class cosmictele():
    def __init__(self):
        return
    def doit(self,Nevt=100,zextrap=150.):
        Nmeas = 4
        fNmeas= float(Nmeas)
        sig   = 1./math.sqrt(12.)
        y0 = 0.
        z0 = 0.
        zv = zextrap+z0
        dz = 15./4. # separation between each measurement plane in z
        Z = numpy.array([z0+float(i)*dz for i in range(Nmeas)])
        Sz = sum(Z)
        Szz= sum(Z*Z)

        yave = yrms = 0.
        for e in range(Nevt):
            Y = numpy.random.normal(y0, sig, Nmeas)
            Sy= sum(Y)
            Syz= sum(Y*Z)
            m = (fNmeas*Syz - Sy*Sz) / (fNmeas*Szz - Sz*Sz)
            b = (Syz - m*Szz)/Sz
            y = m*zv + b
            yave += y
            yrms += y*y
        yave = yave/float(Nevt)
        yrms = math.sqrt( (yrms-yave*yave)/float(Nevt-1) )
        print 'Nevt,yave,yrms',Nevt,yave,yrms,'for zextrap',zextrap,'yrms/zextrap',yrms/zextrap,'z separation of planes',dz,'yres',sig
        return
if __name__ == '__main__' :

    CT = cosmictele()
    for z in [50., 100., 150., 200., 250.]:
        CT.doit(Nevt=100000,zextrap=z)
