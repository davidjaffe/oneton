#!/usr/bin/env python
'''
estimate delta_t distributions for PMT pairs for stopped muons
20221118 

'''
import numpy

import os

import datetime
import Logger
import pickle

import sys
import math
import matplotlib.pyplot as plt
import random

class deltat():
    def __init__(self,debug=-1,nToy=0):
        
### usual directory structure to preserve output(log,figures,pickle,...) of each job
        now = datetime.datetime.now()
        self.now = now.strftime('%Y%m%dT%H%M%S')

        self.figDir = 'DELTAT/'
        
        ### xyz of center of PMT faces. units are mm. from table IV of tech note
        self.pmtPosition = {}
        zbot, ztop = -506.4, 806.4
        self.pmtPosition['S0'] = numpy.array([381.35, -30.79, zbot])
        self.pmtPosition['S1'] = numpy.array([131.35, -30.79, zbot])
        self.pmtPosition['S2'] = numpy.array([-118.35, -30.79, zbot])
        self.pmtPosition['S3'] = numpy.array([-368.65, -30.79, zbot])
        self.pmtPosition['S4'] = numpy.array([-165.03, 280.38, ztop])
        self.pmtPosition['S5'] = numpy.array([-165.03,-128.43, ztop])
        self.pmtPosition['S6'] = numpy.array([381.35, 169.21, zbot])
        self.pmtPosition['S7'] = numpy.array([131.35, 169.21, zbot])

        ### inner dimensions of cylinder
        self.cylZ = [zbot-25.4, ztop+25.4]
        self.dz   = abs(self.cylZ[1]-self.cylZ[0])
        self.cylR = 995.0/2.

        self.twopi = 2.*math.acos(-1.)
        self.clight = 299792458. * 1000. * 1.e-9 # mm/ns
        print('deltat.__init__ clight(mm/ns)',self.clight)
        self.waterIOR = 4./3.
        self.cinwater = self.clight/self.waterIOR
        
        
        return
    def genPos(self,ZMI=None,ZMA=None,XLIM=None,YLIM=None):
        '''
        return random position within the cylinder or 
        within specified limits within cylinder
        '''
        keepGoing = True
        while (keepGoing):
            zmi = self.cylZ[0]
            zma = zmi + self.dz
            if ZMI is not None : zmi = max(ZMI,self.cylZ[0])
            if ZMA is not None : zma = min(ZMA,self.cylZ[1])
            z = random.uniform(zmi,zma)
            r = random.uniform(0.,self.cylR)
            phi = random.uniform(0.,self.twopi)
            y = r*math.sin(phi)
            x = r*math.cos(phi)
            keepGoing = False
            if XLIM is not None:
                if x<XLIM[0] or XLIM[1]<x : keepGoing = True
            if YLIM is not None:
                if y<YLIM[0] or XLIM[1]<y : keepGoing = True
        return numpy.array([x,y,z])
    def btwn(self,U,V):
        '''
        return distance between two three vectors U and V
        '''
        return numpy.linalg.norm(U-V)
    def distToPMTs(self,X):
        '''
        return distance from point X to each PMT as a dict
        '''
        d = {}
        for pmt in self.pmtPosition:
            d[pmt] = self.btwn(X,self.pmtPosition[pmt])
        return d
        
    def main(self):
        Nevt = 1000
        self.plotOneExpt(Nevt=Nevt,words='Uniform')
        self.plotOneExpt(Nevt=Nevt,XLIM=[100.,400.], YLIM=[-150., 200.],words='radial SM')
        zcenter = 0.5*(self.cylZ[0]+self.cylZ[1])
        self.plotOneExpt(Nevt=Nevt,ZMI=self.cylZ[0]+10.,ZMA=zcenter,XLIM=[100.,400.], YLIM=[-150., 200.],words='bottom SM')
        self.plotOneExpt(Nevt=Nevt,ZMI=zcenter,ZMA=self.cylZ[1]-10.,XLIM=[100.,400.], YLIM=[-150., 200.],words='top SM')

        return
    def plotOneExpt(self,Nevt=1000,ZMI=None,ZMA=None,XLIM=None,YLIM=None,words=''):
        '''
        generate, plot deltaT distributions
        '''
        
        ## flight times
        T = {x:[] for x in self.pmtPosition}
        for ievt in range(Nevt):
            X = self.genPos(ZMI=ZMI,ZMA=ZMA)
            d = self.distToPMTs(X)
            for pmt in T: T[pmt].append(d[pmt]/self.cinwater)
        ## convert to numpy array
        tmi,tma = 0.,8.
        for pmt in T:
            T[pmt] = numpy.array(T[pmt])
            tmi = min(tmi,min(T[pmt]))
            tma = max(tma,max(T[pmt]))
        dt = (tma-tmi)/20.
        tmi,tma = tmi-dt,tma+dt
        
        ## plot time distributions of each PMT
        bins = 40
        nrows, ncols = 4,2
        fig, ax = plt.subplots(ncols=ncols,nrows=nrows, sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0,wspace=0)
        for pmt in T:
            ipmt = int(pmt[1])
            irow, icol = ipmt//ncols, ipmt%ncols
            mean, std = numpy.mean(T[pmt]), numpy.std(T[pmt])
            htit = pmt + ' $\mu$ {:.1f} $\sigma$ {:.1f}'.format(mean,std)
            ax[irow,icol].hist( T[pmt], range=(tmi,tma), bins=bins )
            ax[irow,icol].set_title(htit,y=1.0,pad=-14,loc='left')
#            ax[irow,icol].grid()
            ax[irow,icol].tick_params(top=True,right=True,direction='in')
        cmi,cma, cx, cy = 'zmi default','zma default','x default','y default'
        if ZMI is not None: cmi = 'zmi {:.1f}'.format(ZMI)
        if ZMA is not None: cma = 'zma {:.1f}'.format(ZMA)
        if XLIM is not None: cx = '{:.1f}<x<{:.1f}'.format(XLIM[0],XLIM[1])
        if YLIM is not None: cy = '{:.1f}<y<{:.1f}'.format(YLIM[0],YLIM[1])
            
        title = cmi + ' '  + cma + ', ' + cx + ', ' + cy
        
        fig.suptitle(title)
        plt.show()
            
        ## plot time difference distributions of all pairs of PMTs
        bins = 40
        nrows,ncols = 7,4 
        fig, ax = plt.subplots(ncols=ncols,nrows=nrows, sharex=True, sharey=True, figsize=(5,7))
        fig.subplots_adjust(hspace=0,wspace=0)
        i = 0
        npmts = T.keys()
        for pmtA in T:
            iA = int(pmtA[1])
            tA = T[pmtA]
            for pmtB in T:
                iB = int(pmtB[1])
                if iB>iA :
                    tB = T[pmtB]
                    dt = tB-tA
                    irow, icol = i//ncols, i%ncols
                    i += 1
                    mean, std = numpy.mean(dt), numpy.std(dt)
                    htit = pmtB+'-'+pmtA+ ' $\mu$ {:.1f} $\sigma$ {:.1f}'.format(mean,std)
                    ax[irow,icol].hist( dt, bins=bins)
                    ax[irow,icol].set_title(htit,y=1.0,pad=-10,loc='left',fontsize=6)
                    ax[irow,icol].tick_params(top=True,right=True,direction='in')
        cmi,cma, cx, cy = 'zmi default','zma default','x default','y default'
        if ZMI is not None: cmi = 'zmi {:.1f}'.format(ZMI)
        if ZMA is not None: cma = 'zma {:.1f}'.format(ZMA)
        if XLIM is not None: cx = '{:.1f}<x<{:.1f}'.format(XLIM[0],XLIM[1])
        if YLIM is not None: cy = '{:.1f}<y<{:.1f}'.format(YLIM[0],YLIM[1])
            
        title = cmi + ' '  + cma + ', ' + cx + ', ' + cy 
        if words!='' : title += ', '+words
        
        fig.suptitle(title,fontsize=10)
        pdf = self.figDir+title.replace(' ','_')+'.pdf'
        plt.savefig(pdf)
        print('deltat.plotOneExpt Wrote',pdf)
        plt.show()
                    

        return
if __name__ == '__main__' :
    debug = -1

    if len(sys.argv)>1 :
        if 'help' in sys.argv[1].lower():
            print('USAGE: python deltat.py debug[{0}]'.format(debug))
            sys.exit()
    
    if len(sys.argv)>1 : debug = int(sys.argv[1])
    
    P = deltat(debug=debug)
    P.main()

