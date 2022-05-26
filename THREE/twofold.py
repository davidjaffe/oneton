#!/usr/bin/env python
'''
determine efficiency of each signal PMT relative to MC 
using table of unique twofold coincidences from WbLS_0525_2022
20220525
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt

class twofold():
    def __init__(self,debug=-1):

        self.debug = debug

        self.sources   = ['DATA','MC']
        self.TwoFold = {}
        self.TwoFold['MC'] = numpy.array( (0,101836, 62216, 37352, 7560, 10556, 105196, 89348, # S0 in coinc S0-S7
                               101836, 0, 112644, 62832, 13356, 15372, 98756, 131880,
                               62216, 112644, 0, 89936, 16800, 19292, 66584, 101472,
                               37352, 62832, 89936, 0, 16828, 19320, 40768, 60452,
                               7560, 13356, 16800, 16828, 0, 66668, 11508, 15764,
                               10556, 15372, 19292, 19320, 66668, 0, 10304, 13216,
                               105196, 98756, 66584, 40768, 11508, 10304, 0, 115080,
                               89348, 131880, 101472, 60452, 15764, 13216, 115080,0)
                            ).reshape( (8,8) )
        self.TwoFold['DATA']=numpy.array( (0, 4258, 1892, 1940, 875, 720, 4156, 3739,
                               4258, 0, 2550, 2474, 988, 777, 3412, 4344,
                               1892, 2550, 0, 2098, 664, 522, 1538, 2185,
                               1940, 2474, 2098, 0, 999, 735, 1685, 2218,
                               875, 988, 664, 999, 0, 1355, 920, 1064,
                               720, 777, 522, 735, 1355, 0, 622, 723,
                               4156, 3412, 1538, 1685, 920, 622, 0, 4032,
                               3739, 4344, 2185, 2218, 1064, 723, 4032,0)
                               ).reshape( (8,8) )
        print('twofold.__init__ Consistency check')
        self.nPMT = 8
        OK = True
        maximum = {}
        for src in self.sources: 
            L = len(self.TwoFold[src][0])
            if L!=self.nPMT : print('inconsistent array length',L,'and nPMT',self.nPMT)
            m = 0
            for i in range(L):
                for j in range(L):
                    Nij = self.TwoFold[src][i][j]
                    Nji = self.TwoFold[src][j][i]
                    m = max(Nij, m)
                    if Nij!=Nji:
                        OK = False
                        print(src,'inconsistent i',i,'j',j,'Nij',Nij,'Nji',Nji)
            maximum[src] = m
        if OK : print('twofold.__init__ Consistency check OK')
        self.maximum = maximum

        X = self.TwoFold['MC']
        m = maximum[src]
        if self.debug > 0 : print('twofold.__init__ self.TwoFold[`MC`]',X)
        self.prob = X/m
        if self.debug > 0 : print('twofold.__init__ prob',self.prob)
        

        self.figDir = 'TWOFOLD_FIGURES/'
			
        return
    def chisqr(self,param):
        '''
        chisquared = sum_i=0,7 sum_j=i+1,7 (Cij - effi*effj*Pij*N)^2 / Cij
        C[i,j] = count rate in date for coincidence between PMTs i and j
        P[i,j] = probability of coincidence between PMTs i,j
        effi = efficiency of PMT i ( = param[0:8] )
        N = total number of data events ( = param[8] )
        '''
        eff = param[0:8]
        N   = param[8]
        cs = 0.
        C = self.TwoFold['DATA']
        P = self.prob
        for i in range(self.nPMT):
            for j in range(i+1,self.nPMT):
                num = (C[i,j] - eff[i]*eff[j]*P[i,j]*N)
                cs = num*num/C[i,j]
        return cs
    def main(self):

        return
if __name__ == '__main__' :
    P = twofold()
    P.main()
    
