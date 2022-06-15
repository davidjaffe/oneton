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
from scipy.optimize import minimize

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
        self.Npmt = self.nPMT = 8
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
        print('twofold.__init__ maximum',maximum)
        

        self.figDir = 'TWOFOLD_FIGURES/'
			
        return
    def oneExpt(self):
        '''
        return 'Data' for one generated experiment given input efficiencies and number of events
        '''
        NN = self.Npmt*self.Npmt
        eff = self.input[:self.Npmt]
        epp = []
        for i in range(self.Npmt):
            for j in range(self.Npmt):
                epp.append( eff[i]*eff[j]*self.prob[i,j] )
        epp = numpy.array(epp)
        X = numpy.zeros( NN )
        Ones = numpy.ones( NN )
        N   = self.input[self.Npmt]
        for evt in numpy.arange(N):
            R = numpy.random.random( NN )
            X += (R<=epp)*Ones
        X = numpy.fix( X )
        X = X.reshape(self.Npmt,self.Npmt)
        return X
    def oneTrueExpt(self):
        '''
        return 'Data' for one generated experiment given input efficiencies and number of events
        '''
        eff = self.input[:self.Npmt]
        N   = self.input[self.Npmt]
        X = numpy.zeros( self.Npmt*self.Npmt ).reshape( (self.Npmt,self.Npmt) )
        for i in range(self.Npmt):
            for j in range(self.Npmt):
                v = 0. 
                if i!=j : v = eff[i]*eff[j]*self.prob[i,j]*N
                X[i,j] += v
        X = numpy.fix( X )
        return X
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
        C = self.inputData # TwoFold['DATA']
        P = self.prob
        for i in range(self.nPMT):
            for j in range(i+1,self.nPMT):
                num = (C[i,j] - eff[i]*eff[j]*P[i,j]*N)
                cs += num*num/C[i,j]
        return cs
    def find(self,method='Powell',Nguess = 5000.):
        param = [1. for x in range(self.nPMT)]
        bounds = [(None,None) for x in range(self.nPMT)]
        bounds.append( (self.maximum['DATA'], 2.*Nguess) )
        param.append(Nguess)
        chi2 = self.chisqr(param)
        if self.debug > 0: 
            print('\ntwofold.find method',method)
            print('twofold.find chi2',chi2,'with input params',param)
        res = minimize(self.chisqr, param, method=method,bounds=bounds)
        if self.debug > 1 : print('twofold.find res',res)
        pout = res.get('x')
        hess_inv = res.get('direc') #???
        chi2 = self.chisqr(pout)
        fitpar = 'twofold.find {0} chi2 {1:.1f} fitpar '.format(method,chi2)
        if self.debug > 1:
            print('twofold.find res.keys()',res.keys())
            print('twofold.find hess_inv',hess_inv)
        for i,p in enumerate(pout):
            punc = 0.
            if hess_inv is not None: punc = math.sqrt(max(0.,hess_inv[i,i]))
            fitpar += ' {0:.3f}({1:.3f})'.format(p,punc)

        print(fitpar)
        return pout
    def analyzeResults(self,fitresults):
        '''
        compare fitresults with input for toyMC
        '''
        print('\ntwofold.analyzeResults')
        L = len(self.input)
        means,stddevs = {},{}
        for method in fitresults:
            means[method],stddevs[method] = [],[]
            results = numpy.array( fitresults[method] )
            #print('method,results',method,results)
            Mresult = method + ' means'
            Sresult = method + 'stddev'
            for i in range(L):
                #if i==0: print('i,results[:,i]',i,results[:,i])
                m,s = numpy.mean( results[:,i] ), numpy.std(results[:,i] )
                means[method].append( m )
                stddevs[method].append( s )
                Mresult += ' {:.3f}'.format(m)
                Sresult += ' {:.3f}'.format(s)
            print(Mresult)
            print(Sresult)
        Winput = 'Input pararameters'
        for i in range(L):
            Winput += ' {:.3f}'.format(self.input[i])
        print(Winput)

        return
    def main(self):

        self.debug = -1

        ## these are minimization methods for scipy 1.8.1
        ## Unconstrained: CG, BFGS 
        ## bound-constrained : 'Nelder-Mead'
        methods = []
        unconstrained_methods = ['CG','BFGS']
        boundconstrained_methods = ['Nelder-Mead','Powell','L-BFGS-B']
        methods.extend( unconstrained_methods )
        methods.extend( boundconstrained_methods )
        
        toyMC = True
        Nexpt = 1
        if toyMC:
            Nexpt = 1000
            Nguess = 5000.
            self.input = [1.2 - numpy.random.random()*0.4 for x in range(self.Npmt)]
            print('twofold.main toyMC input effy',' %.2f'*len(self.input)%tuple(self.input))
            self.input.append( Nguess )
            
        fitresults = {m:[] for m in methods}
        for expt in range(Nexpt):
            print('twofold.main expt',expt)
            if toyMC:
                self.inputData = self.oneExpt()
            else:

                self.inputData = self.TwoFold['DATA']
                Nguess = self.maximum['DATA']

            for method in methods:
                fitpar = self.find(method=method,Nguess=Nguess)
                fitresults[method].append( fitpar )
            #print('fitresults',fitresults)
        if toyMC : self.analyzeResults(fitresults)
        return
if __name__ == '__main__' :
    P = twofold()
    P.main()
    
