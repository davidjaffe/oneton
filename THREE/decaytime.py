#!/usr/bin/env python
'''
plot and integration decay time distribution from 2003.10491
20221229
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import numpy

class decaytime():
    def __init__(self):

        # from table 2: (rise fraction, time(ns) ),3x (decay fraction, time)
        self.decaypar = {'1%': [ (1.0, 0.23), (0.87, 2.00), (.068, 12.0), (0.062, 110) ] }
        
        return
    def decayfunc(self,t,conc='1%',use=[1.,1.,1.,1.]):
        '''
        return decayfunc = sum(f*exp(-t/tau)) with t and tau in ns
        '''
        func = -1.*use[0]

        for i,pair in enumerate(self.decaypar[conc]):
            f,tau = pair
            if i==0:
                func += use[i]*(1.-math.exp(-t/tau))
            else:
                func += use[i]*f*math.exp(-t/tau)
        return func
    def main(self):

        tmi,tma,nt = 0.,500.,5000
        tlim = 60.
        T = numpy.linspace(tmi,tma,nt)
        
        f = []
        f1,f2,f3 = [],[],[] 
        facc = []
        s = 0.

        for i,t in enumerate(T):

            y = self.decayfunc(t)
            f.append(y)
            s += y
            facc.append(s)

            f1.append( self.decayfunc(t,use=[0.,1.,0.,0.]) )
            f2.append( self.decayfunc(t,use=[0.,0.,1.,0.]) )
            f3.append( self.decayfunc(t,use=[0.,0.,0.,1.]) )

        f = numpy.array(f)
        f1,f2,f3 = numpy.array(f1),numpy.array(f2),numpy.array(f3)
        m = facc[-1]
        facc = numpy.array(facc)/m
        for t in [25.,50.]:
            i = numpy.argmin(abs(T-t))
            print('Cumulative integral at',i,'t=',T[i],'ns is',facc[i])
        

        uselog = [0,0,1,1]
        xma    = [-1,tlim,-1,tlim]
        for i in range(4):
            
            plt.plot(T,f,'r')
            plt.plot(T,f1,linestyle='dotted',color='grey')
            plt.plot(T,f2,linestyle='dotted',color='grey')
            plt.plot(T,f3,linestyle='dotted',color='grey')
            plt.plot(T,facc,'b')
            plt.grid()
            plt.ylim(1.e-4,1.1)
            if uselog[i]==1: plt.yscale('log')
            if xma[i]>0    : plt.xlim(-1.,xma[i])    
            plt.show()

        
        return
if __name__ == '__main__' :
    P = decaytime()
    P.main()
    
