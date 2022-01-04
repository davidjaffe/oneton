#!/usr/bin/env python
'''
analyze cry input file
20220104
'''
import numpy
#import graphUtils
#import ROOT
#from ROOT import TFile,TH1D,TH2D,gDirectory
#import gfit
#import pipath
import os
import datetime
import sys
import math
import matplotlib.pyplot as plt
import mpl_interface

import Logger

class cry_anal():
    def __init__(self):
       #self.mpl_interface = mpl_interface.mpl_interface(internal=True)
       return
    def main(self):
        '''

        '''
        fn = 'CRY_DATA/out_20210823_1.dat'
        f = open(fn,'r')
        mom = {13:[], -13:[]}
        pmi,pma = 1.e6,-1.e6
        for line in f:
            if len(line)>10:
                X = line.split()
                pid = int(X[1])
                px  = float(X[4])
                py  = float(X[5])
                pz  = float(X[6])
                p = math.sqrt(px*px+py*py+pz*pz)
                
                if abs(pid)==13: 
                    mom[pid].append(p)
                    pmi = min(p,pmi)
                    pma = max(p,pma)
        f.close()
        r = float(len(mom[-13]))/float(len(mom[13]))
        print('len(mom[13]),len(mom[-13]),r',len(mom[13]),len(mom[-13]),r)
        print('pmi,pma',pmi,pma)
        nbin,xlo,xhi = 20,0.,40.
        bins = numpy.array((0.,1.,2.,4.,8.,16.,32.,64.))#,1024.))
        phist,pedges = numpy.histogram(mom[-13],bins=bins)
        mhist,medges = numpy.histogram(mom[13],bins=bins)
        rhist = phist/mhist
        drhist = rhist*numpy.sqrt(1./phist + 1./mhist)
        for a,b in zip(pedges,medges):
            print('pedge {0:.3f} medge {1:.3f}'.format(a,b))
        for a,b in zip(rhist,drhist):
            print('{0:.3f} {1:.3f}'.format(a,b))
        x = (bins[:-1]+bins[1:])/2.
        dx= (bins[1:]-bins[:-1])/2.
        plt.errorbar(x,rhist,yerr=drhist,xerr=dx,linestyle='')
        plt.xlabel('Muon momentum (GeV/c)')
        plt.ylabel('Positive to negative muon ratio')
        plt.grid()
        pdf = 'cry_muon_charge_ratio_vs_momentum.pdf'
        plt.savefig(pdf)
        print('cry_anal.main Saved figure to',pdf)
        plt.show()
#        self.mpl_interface.histo(mom[13],xlo,xhi,dx=2.,title='Positive muons',xlabel='P(GeV)',grid=True)
        return
if __name__ == '__main__' :
    ca = cry_anal()
    ca.main()
