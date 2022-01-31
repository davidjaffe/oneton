#!/usr/bin/env python
'''
optics of 1ton
try to estimate effect of reflectivity of black barrier
20220125
'''
import numpy

import os
#import datetime
import sys
import math
import matplotlib.pyplot as plt
#import mpl_interface

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
if __name__ == '__main__' :
    ca = optics()
    ca.main()
