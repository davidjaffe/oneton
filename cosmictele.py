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
    def doit(self,Nevt=100,zextrap=150.,Nplane=1,dzplane2=100.):
        Nlayers = 4
        
        Nmeas = Nlayers
        sig   = 1./math.sqrt(12.)
        y0 = 0.
        z0 = 0.
        zv = zextrap+z0 # Z-position to calculate Y based on measurements
        dz = 15./4. # separation between each measurement plane in z
        Z = [z0+float(i)*dz for i in range(Nmeas)]
        Z0plane = [Z[0]]
        if Nplane==2:

            Nmeas += Nlayers
            Z.extend([z0+dzplane2+float(i)*dz for i in range(Nlayers)])
            Z0plane.append(Z[Nlayers])
                           
        fNmeas= float(Nmeas)
        Z = numpy.array(Z)
        Sz = sum(Z)
        Szz= sum(Z*Z)
        #print 'Z',Z

        yave = yrms = 0.
        slope = slrms = 0.
        bave = brms  = 0.
        for e in range(Nevt):
            Y = numpy.random.normal(y0, sig, Nmeas)
            #if e==0: print 'Y',Y
            Sy= sum(Y)
            Syz= sum(Y*Z)
            m = (fNmeas*Syz - Sy*Sz) / (fNmeas*Szz - Sz*Sz)
            b = (Syz - m*Szz)/Sz
            y = m*zv + b
            yave += y
            yrms += y*y
            slope += m
            slrms += m*m
            bave += b
            brms += b*b
        yave = yave/float(Nevt)
        yrms = math.sqrt( (yrms-yave*yave)/float(Nevt-1) )
        slope = slope/float(Nevt)
        slrms = math.sqrt( (slrms-slope*slope)/float(Nevt-1) )
        bave = bave/float(Nevt)
        brms = math.sqrt( (brms-bave*bave)/float(Nevt-1) )
        words = 'Nevt '+ str(Nevt)+' #planes '+str(Nplane)+' Z0plane '
        words += ' '.join('%.1f'% x for x in Z0plane)
        words += ' yave,yrms,zextrap '
        words += ' '.join('%.3f'% x for x in [yave,yrms,zextrap])
        words += ' slope,sloperms '
        words += ' '.join('%.5f'% x for x in [slope,slrms])
        words += ' intercept,intrms '
        words += ' '.join('%.5f'% x for x in [b,bave])
        words += ' yres {0:.3f}'.format(sig)
        print words
        #print '#planes',Nplane,'Z0plane',Z0plane,'Nevt,yave,yrms',Nevt,yave,yrms,'for zv',zv,'yrms/zextrap',yrms/zextrap,'z separation of planes',dz,'yres',sig
        return
if __name__ == '__main__' :

    CT = cosmictele()

    for nplane in [1,2]:
        dzplaneList = [0.]
        if nplane==2: dzplaneList = [-100., 600., 1000.]
        for z in [50.,  150.,  250.]:
            print ''
            for dzplane in dzplaneList:
                CT.doit(Nevt=100000,zextrap=z,Nplane=nplane,dzplane2=dzplane)
