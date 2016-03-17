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
import sys

class cosmictele():
    def __init__(self):
        # see notes 20160314 (back of 20160204/3)
        self.zTop        = 184.48 # height of active volume above floor
        self.zHeight     = 125.0  # height of active volume
        self.dzFloor     = 350.   # assumed height of a floor
        self.zCSC1TozTop = 59.48
        self.zCSC2ceiling = self.zCSC1TozTop + 100.
        self.zCSC1first   = self.zCSC1TozTop - (self.zTop + self.dzFloor)
        self.zCSC1basement= self.zCSC1TozTop - (self.zTop + 2.*self.dzFloor)
        
        return
    def pathAngle(self,d=99.5,h=125.,tC=math.acos(1./1.33)):
        '''
        CSC is 15cm high
        compute position of track at top of each CSC
        assume CSC on top of dark box
        2d CSC either under ceiling, on 1st floor, in basement
        # assumes 350cm/floor
        '''
        #        dark
        #        box    under 
        #        top    ceiling      on 1st floor         in basement
        CSCdz = []
        CSCdz.append( self.zCSC1TozTop )   # top of CSC at top of dark box
        CSCdz.append( self.zCSC2ceiling )  # under 2d floor ceiling
        CSCdz.append( self.zCSC1first )    # on 1st floor 
        CSCdz.append( self.zCSC1basement ) # in basement
        words = ' '.join('%10.1f' % x for x in CSCdz)
        words += ' (cm) relative to top of active 1ton volume'
        print '{0:>45} {1:}'.format('',words)
        print '{0:5} {1:5} {2:10} {3:10} {4:>10} {5:>10} {6:>10} {6:>10} {6:>10}'.format('x','y','theta(rad)','theta(deg)','path(cm)','CSC1r0','CSC2r0')
        fractions = [1./4., 1./3., 0.4]
        
        for fx in fractions:
            x = fx*d
            for fy in fractions:
                y = fy*d
                theta = math.atan( h/(d-(x+y)) ) - (math.pi/2.-tC)
                l = x/math.sin(theta)
                print '{0:5.2f} {1:5.2f} {2:10.3f} {3:10.1f} {4:10.1f}'.format(x,y,theta,theta*180./math.pi,l),
                for c in CSCdz:
                    g = c*math.tan(theta)
                    dz = d/2. - (x+g)
                    print ' {0:>10.1f}'.format(dz),
                print ''
                #print 'fx,fy,theta(rad),theta(deg),mupath',fx,fy,theta,theta*180/math.pi,l
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
    mode = 0
    if len(sys.argv)>0 : mode=int(sys.argv[1])

    
    CT = cosmictele()

    if mode==0 or mode==2:

        for nplane in [1,2]:
            dzplaneList = [0.]
            if nplane==2: #dzplaneList = [-100., 600., 1000.]
                dzplaneList = [-(CT.zCSC2ceiling - CT.zCSC1TozTop),  -(CT.zCSC1first - CT.zCSC1TozTop), -(CT.zCSC1basement- CT.zCSC1TozTop)]
            for z in [CT.zCSC1TozTop, CT.zCSC1TozTop+CT.zHeight]:
                print ''
                for dzplane in dzplaneList:
                    CT.doit(Nevt=100000,zextrap=z,Nplane=nplane,dzplane2=dzplane)
    if mode==1 or mode==2:
        CT.pathAngle(h=CT.zHeight)
