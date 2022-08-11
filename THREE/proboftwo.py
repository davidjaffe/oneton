#!/usr/bin/env python
'''
determine probably of a two-fold coincidence between PMTs given the mean PE per PMT for 3 values of simulated light yield (LY)
data taken from WbLS_0802_2022 for 10k p=0 muons generated in the center of the detector.
20220811
'''
import numpy

import os

import datetime
import Logger
import pickle

import sys
import math
import matplotlib.pyplot as plt

from scipy.stats import poisson

import iminuit

class proboftwo():
    def __init__(self,debug=-1,nToy=0):

        if nToy>0: 
### usual directory structure to preserve output(log,figures,pickle,...) of each job
            now = datetime.datetime.now()
            self.now = now.strftime('%Y%m%dT%H%M%S')

            self.rootDir = 'JOBS_PROBOFTWO/'
            parentDir = self.rootDir+self.now
            dirs = [parentDir]
            self.Figures = parentDir  + '/FIGURES/'
            self.logDir = parentDir
            dirs.append( self.Figures)
            dirs.append( self.logDir)

            for d in dirs:
                if not os.path.exists(d):
                    os.makedirs(d)
                    print('calibf.__init__ create directory',d)

            lf = self.logDir + '/logfile.log'
            sys.stdout = Logger.Logger(fn=lf)
            print('proboftwo.__init__ Output directed to stdout and',lf)
            self.figDir = self.Figures #'PROBOFTWO_FIGURES/'

## process, report inputs
        self.debug = debug
        self.nToy = nToy
        self.toyMC = nToy > 0
        print('proboftwo.__init__ debug',self.debug,'nToy',self.nToy,'toyMC',self.toyMC)


## table of mean PE per PMT for LY = 300, 500, 700 from WbLS_0802_2022
        self.meanPE = {300: [1.02, 0.65, 0.59, 0.50, 0.45, 0.17, 1.06, 0.79],
                       500: [1.09, 0.83, 0.66, 0.69, 0.53, 0.24, 1.08, 0.97],
                       700: [1.19, 0.98, 0.70, 0.89, 0.65, 0.33, 1.18, 1.19]}
        self.coincRate = {300: [1160, 1203, 1259],
                          500: [1927, 1994, 1989],
                          700: [2485, 2501, 2523]}  # EVEN, RED, INNER coincidences in 10k events from WbLS_0802_2022
        self.coincDef  = {'EVEN': [0,2,4,6],
                          'RED' : [0,2,5,7],
                          'INNER':[1,2,5,7]} # definitions of coincidence criteria
        self.nSim = 10000 # number of simulated zero momentum muon events per LY
        self.Npmt = self.nPMT = 8
        self.MEANS = {}

			
        return
    def readPickle(self,timestamp):
        '''
        return input & fitresults given job specified by timestamp
        '''
        pickle_fn = self.rootDir + timestamp + '/' + timestamp + '.pickle'
        
        f = open(pickle_fn,'rb')
        Input = pickle.load(f)
        fitresults = pickle.load(f)
        f.close()
        print('proboftwo.readPickle Read',pickle_fn)
        
        return Input, fitresults
    def writePickle(self,Input, fitresults):
        '''
        write inputData & fitresults to pickle file
        order of dump must be matched in load according to https://stackoverflow.com/questions/20716812/saving-and-loading-multiple-objects-in-pickle-file
        '''
        timestamp = self.now
        pickle_fn = self.rootDir + timestamp + '/' + timestamp + '.pickle'


        f = open(pickle_fn, 'wb')
        pickle.dump(Input, f)
        pickle.dump(fitresults, f)
        f.close()

        print('proboftwo.writePickle Wrote',pickle_fn)
        return
    def oneEvent(self,LY,coincType):
        '''
        return list of generated number of PE for each PMT given 
        LY and coincType which determines meanPE per PMT and PMTs to consider for coincidences
        '''
        means = self.REDUCE(LY,coincType)
        genPE = []
        for m in means:
            genPE.append( poisson.rvs(m) )
        return genPE
    def coinc(self,genPE,criterion='>1'):
        '''
        return True if genPE satisfies requirement of >1 hit in 2 or more PMTs
        '''
        if criterion=='>1':
            c = sum([x>0 for x in genPE])>1
        elif criterion=='>2':
            c = sum([x>0 for x in genPE])>2
        elif criterion=='>3':
            c = sum([x>0 for x in genPE])>3
        elif criterion=='==1':
            c = sum([x>0 for x in genPE])==1
        elif criterion=='==2':
            c = sum([x>0 for x in genPE])==2
        elif criterion=='==3':
            c = sum([x>0 for x in genPE])==3
        elif criterion=='==4':
            c = sum([x>0 for x in genPE])==4
        else:
            sys.exit('proboftwo.coinc ERROR Unknown criterion '+criterion)
        return c
    def oneExpt(self,LY,coincType,Nevt=-1,nCoinc=['>1']):
        '''
        return dict with the number of events with nCoinc coincidences among Nevt events 
        given meanPE determined by LY and coincType 
        using MC method
        '''
        if Nevt==-1: Nevt = self.nSim
        Ncoinc = {}
        for cr in nCoinc:
            Ncoinc[cr] = 0
        for evt in range(Nevt):
            genPE = self.oneEvent(LY,coincType)
            for cr in nCoinc:
                if self.coinc(genPE,criterion=cr) : Ncoinc[cr] += 1
        return Ncoinc
    def reduce(self,meanPE,cDef):
        means = []
        for i in cDef: means.append( meanPE[i] )
        return means
    def REDUCE(self, LY, coincType):
        key = str(LY)+coincType
        if key in self.MEANS : return self.MEANS[key]
        means = self.reduce(self.meanPE[LY], self.coincDef[coincType] )
        self.MEANS[key] = means
        return means
    def calculate(self,LY,coincType):
        '''
        determine probability of two-fold coincidence for coincType and LY
        '''
        if LY not in self.meanPE: sys.exit('proboftwo.calculate ERROR LY ' + str(LY) + ' not valid')
        if coincType not in self.coincDef: sys.exit('proboftwo.calculate ERROR coincType ' + coincType + ' not valid')
            
        means = self.REDUCE(LY,coincType)

        p2f = self.pof2(means)
        return p2f
    def pof2(self,means):
        '''
        return probability of two-fold coincidence given mean PE in means
        '''
        L = len(means)
        Pnohits = 1.
        for m in means:
            Pnohits *= poisson.pmf(0,m)
        P1hits  = 0.
        for j in range(L):
            pj = poisson.sf(0,means[j])
            for i in range(L):
                if i!=j : pj *= poisson.pmf(0,means[i])
            P1hits += pj
        if self.debug>0 : print('proboftwo.pof2 Pnohits,P1hits',Pnohits,P1hits)
        pof2 = 1. - (Pnohits + P1hits)
        return pof2
    def ebinom(self,N,Ntot):
        '''
        return binomial error for N successes out of Ntot tries
        '''
        if N<=0 or N>=Ntot: return 0.
        f = float(N)/float(Ntot)
        return math.sqrt(f*(1.-f)/float(Ntot))*float(Ntot)
    def main(self):
        '''
        calculation and toy simulation for two-fold coincidences
        '''
        coincCriteria = [ '>1', '>2', '>3', '==1', '==2', '==3','==4']
        favCC = coincCriteria[0]

        results = {}

        for LY in self.meanPE:
            results[LY] = []

            for coincType in self.coincDef:
                p2f = self.calculate(LY,coincType)
                N = int(p2f*self.nSim)
                results[LY].append(N)

                Ncoinc = self.oneExpt(LY,coincType,Nevt=self.nSim,nCoinc=coincCriteria)
                e = self.ebinom(Ncoinc[favCC],self.nSim)
                s = -999.
                if e>0: s = float(N-Ncoinc[favCC])/e
                print('proboftwo.main LY,coincType,N,Ncoinc,signif',LY,coincType,N,Ncoinc[favCC],'{:.2f}'.format(s),'coincidence criteria is',favCC)
                if len(coincCriteria)>1:
                    words = ''
                    for CC in coincCriteria: words += CC + ':{}, '.format(Ncoinc[CC])
                    print('proboftwo.main',words)

        h = '{:>3}'.format('LY')
        for coincType in self.coincDef: h+= ' {:>7}'.format(coincType)
        print(h)
        for LY in self.meanPE:
            h = '{:>3}'.format(LY)
            for i in range(len(self.coincDef)):
                h += ' {:>7}'.format(results[LY][i])

            print(h)
        return
if __name__ == '__main__' :
    debug = -1
    nToy  = 0
    timestamp = None
    if len(sys.argv)>1 :
        if 'help' in sys.argv[1].lower():
            print('USAGE: python proboftwo.py debug[{0}] nToy[{1}] timestamp[{2}]'.format(debug,nToy,timestamp))
            print('USAGE: if timestamp is not None, then analyze data corresponding to timestamp')
            print('USAGE: if timestamp==`MANY`, then analyze data from many timestamps')
            sys.exit()
    
    if len(sys.argv)>1 : debug = int(sys.argv[1])
    if len(sys.argv)>2 : nToy  = int(sys.argv[2])
    if len(sys.argv)>3 : timestamp = sys.argv[3]
    
    P = proboftwo(debug=debug,nToy=nToy)
    if timestamp is None:
        P.main()
    elif timestamp=='MANY':
        P.MANY()
    else:
        P.readAndAnalyze(timestamp)
    
