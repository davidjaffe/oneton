#!/usr/bin/env python
'''
optics of 1ton
try to estimate effect of reflectivity of black barrier
20220125

also calculate path-lengths on inner cylinder surface and from volume
20220208
'''
import numpy

import os
#import datetime
import sys
import math
import matplotlib.pyplot as plt
#import mpl_interface
import random

#import Logger

class optics():
    def __init__(self):
        self.ior_water = 1.331
        self.ior_acrylic = 1.506
        self.tankD = 995. # mm, tank ID
        self.tankH = 1250. # mm, tank, inner height
        
        return
    def rt(self,Ti,n1,n2):
        '''
        return reflectivity and transmission given
        input angle Ti 
        input, output index of refraction n1,n2
        '''
        sinTt = n1/n2*math.sin(Ti)
        if abs(sinTt)>1. :
            R,T = 1.,0.
        else:
            Tt = math.asin(sinTt)
            alpha = math.cos(Tt)/math.cos(Ti)
            beta = n2/n1
            den = 1. + alpha*beta
            R = (1. - alpha*beta)/den
            R = R*R
            T = alpha*beta*2.*2./den/den
        return R,T
    def main(self):
        '''

        '''
        n1 = self.ior_water
        Tc1 = math.acos(1./n1)
        n2 = self.ior_acrylic
        Tc2 = math.acos(1./n2)
        
        piby2 = math.pi/2.
        
        Ti = numpy.linspace(2./9.*piby2,piby2,100)
        R,T,H = [],[],[]
        x,z0 = 0.,0.
        xr = self.tankD/2.
        for q in Ti:
            r,t = self.rt(q,n1,n2)
            R.append(r)
            T.append(t)
            Tr = q
            h = (xr-x)*math.tan(piby2-Tr) + z0
#            print('q,r,t,h',q,r,t,h)
            H.append(h/self.tankH)
        R = numpy.array(R)
        T = numpy.array(T)
        H = numpy.array(H)
        yma1 = max(H)*1.1
        yma2 = max(R*H)*1.1
        ymi = 0.
        for ymax in [yma1,yma2]:
            plt.plot(Ti,R,'r-',label='reflectivity')
            plt.plot(Ti,T,'b-',label='transmission')
            plt.plot(Ti,H,'g-',label='normalized height')
            plt.plot(Ti,H*R,'c:',label='reflectivity*normed height')
            plt.grid()
            plt.legend()
            plt.axvline(Tc1,label='Water Cerenkov angle',color='black')
            title = 'water to acrylic boundary'
            plt.xlabel('Incident angle')
            plt.ylim(ymi,ymax)
            savefig = False
            if savefig : 
                pdf = 'reflectivity_transmission'+ title.replace(' ','_') + '.pdf'
                print('optics.main Saved figure to',pdf)
                plt.savefig(pdf)
            plt.show()


        return
    def paths(self):
        '''
        compare average path length from random points in a cylinder to random points on the surface of that cylinder

        center of cylinder is origin
        PMT is at phi=0, radius = R/2, top surface
        '''
        R = self.tankD/2. # mm, tank radius
        H = self.tankH  # mm, tank, inner height
        Atop = Abot = numpy.pi*R*R # area of top or bottom
        Aside= numpy.pi*R*2.*H # area of cylindrical wall
        ftop = Atop/(Atop+Abot+Aside)

        x = R/2.
        y = 0.
        z = H/2.
        PMT = numpy.array([x,y,z])
        Origin = numpy.array([0.,0.,0.])
        toOrigin = numpy.linalg.norm(PMT-Origin)
        print('distance of PMT to origin(mm)',toOrigin)

        
        npoints = 10000
        # points in cylindrical volume
        Vr = numpy.random.uniform(0.,R,npoints)
        Vct= numpy.random.uniform(-1.,1.,npoints)
        Vphi=numpy.random.uniform(0.,2.*numpy.pi,npoints)
        Vst= numpy.sqrt(1.-Vct*Vct)
        Vx = Vr*Vst*numpy.cos(Vphi)
        Vy = Vr*Vst*numpy.sin(Vphi)
        Vz = Vr*Vct
        dV = []
        for i in range(npoints):
            point = numpy.array( [Vx[i],Vy[i],Vz[i]] )
            dist = numpy.linalg.norm(PMT-point) 
            dV.append( dist*dist )
        dV = numpy.array(dV)
        print('Volume distance^2 mean',dV.mean(),'stddev',dV.std(),'points',npoints)
        volx = numpy.copy(Vx)
        voly = numpy.copy(Vy)
        volz = numpy.copy(Vz)

        npbot = nptop = int(npoints * ftop)
        npside = npoints - npbot - nptop
        nptop = 0
        Ntot = npbot + nptop + npside

        # top
        N = nptop
        Vz = numpy.array([H/2. for q in range(N)])
        Vr = numpy.random.uniform(0.,R,N)
        Vphi=numpy.random.uniform(0.,2.*numpy.pi,N)
        Vx = Vr*numpy.cos(Vphi)
        Vy = Vr*numpy.sin(Vphi)

        # bottom
        N = npbot
        Vz= numpy.append(Vz, numpy.array([-H/2. for q in range(N)]) )
        Vr = numpy.random.uniform(0.,R,N)
        Vphi=numpy.random.uniform(0.,2.*numpy.pi,N)
        Vx = numpy.append(Vx, Vr*numpy.cos(Vphi) )
        Vy = numpy.append(Vx, Vr*numpy.sin(Vphi) )

        # side
        N = npside
        Vr = numpy.array([R for q in range(N)])
        Vz = numpy.append(Vz, numpy.random.uniform(-H/2.,H/2,N) )
        Vphi = numpy.random.uniform(0.,2.*numpy.pi,N)
        Vx = numpy.append(Vx, Vr*numpy.cos(Vphi) )
        Vy = numpy.append(Vy, Vr*numpy.sin(Vphi) )
        sidex = numpy.copy( Vx )
        sidey = numpy.copy( Vy )
        sidez = numpy.copy( Vz )

        dS =[]
        for i in range(Ntot):
            point = numpy.array( [Vx[i],Vy[i],Vz[i]] )
            dist = numpy.linalg.norm(PMT-point) 
            dS.append( dist*dist )
        dS = numpy.array( dS )
        print('Side distance^2 mean',dS.mean(),'stddev',dV.std(),'points',Ntot,'top points',nptop)
        print('volume/side {:.2f}, corrected for misses {:.2f}'.format(dV.mean()/dS.mean(),(dV.mean()/dS.mean())/(float(Ntot)/float(npoints))))

        #print('dV',dV)
        #print('dS',dS)
        for L in [5.e3, 10.e3, 20.e3]:
            atten = []

            ncomb = 10000
            for j in range(ncomb):
                ivol = numpy.random.randint(0,npoints)
                iside= numpy.random.randint(0,Ntot)
                emit = numpy.array( [volx[ivol], voly[ivol], volz[ivol]] )
                refl = numpy.array( [sidex[iside], sidey[iside], sidez[iside]] )
                lre = numpy.linalg.norm(emit-refl)
                lde = numpy.linalg.norm(emit-PMT)
                lrd = numpy.linalg.norm(refl-PMT)
                att = numpy.exp(-lre/L)*numpy.exp(-lrd/L)*numpy.exp(lde/L)
                atten.append( att )
            atten = numpy.array(atten)
            print('average attenuation factor {:.3f} stddev {:.3f} for {} combinations assuming {:.1f} mm attenuation length'.format(atten.mean(),atten.std(),ncomb,L))
        
        
        
        
        
        return
if __name__ == '__main__' :
    ca = optics()
    if False:
        ca.main()
        sys.exit('done')
    ca.paths()
