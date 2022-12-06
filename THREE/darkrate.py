#!/usr/bin/env python
'''
calculate dark rate from Darkcount_LED.pdf histograms.
Rong says this is LED triggers from randomly selected runs in water data
20220308
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt

class darkrate():
    def __init__(self):
        self.Ntrig = 677736 # take from # of entries
        self.hmean = [0.0005518, 0.0005547, 0.0005303, 0.0007668, 0.0003655, 0.0002035, 0.0007539, 0.0008614] # mean of hists. S0, S1...S7
        self.dt = 1500.e-9 # integration window per trigger

        integrated_time = float(self.Ntrig) * self.dt
        print('Total integrated time(s) {:.2f}'.format(integrated_time))
        return
    def poisson(self,mu,n):
        '''
        return f = mu^n * exp(-mu) / n!
        '''
        return math.pow(mu,n) * math.exp(-mu) / math.factorial(n)
    def main(self):
        rate = []
        nhit = []
        for i,x in enumerate(self.hmean):
            nhit.append( x * self.Ntrig )
            r = x/self.dt
            rate.append( r )
            dist = []
            f = 1.
            mu = x
            n = 0
            while f>0.1 :
                f = self.Ntrig * self.poisson(mu,n)
                dist.append(f)
                n += 1
                
            print('S{:}:'.format(i),'expected cts/bin ',' '.join(['{:.1f}'.format(q) for q in dist]))
        print('Rate(Hz) ',', '.join(['{:.1f}'.format(x) for x in rate]))
        print('Hits/PMT ',', '.join(['{:.0f}'.format(x) for x in nhit]))
        phit = []
        tprob = 40. / (2000.-600.) # prob of overlap of 40ns coinc in 1400ns selection
        for i,x in enumerate(self.hmean):
            p = 0.
            for j,y in enumerate(self.hmean):
                p += (1. - math.exp(-y))
            p = p * (1. - math.exp(-x))*tprob
            phit.append( p )
        print('2fold coinc prob ',' '.join(['{:.1g}'.format(q) for q in phit]),'sum ','{:.1g}'.format(sum(phit)))
        
        return
if __name__ == '__main__' :
    P = darkrate()
    P.main()
    
