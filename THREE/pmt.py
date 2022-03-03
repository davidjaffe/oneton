#!/usr/bin/env python
'''
look at PMTCHARGE in RAT implementation 
https://github.com/RongZhao-arlier/WbLS_One_Ton/blob/main/ratpac/data/20210322watermc/PMTCHARGE.ratdb

20220303
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt

class pmt():
    def __init__(self):
        self.ratPMTChargeFile = 'PMTCHARGE.ratdb' # downloaded 20220303
        self.Figures = 'PMT_FIGURES/'
        return
    def readCharge(self):
        '''
        read pmt spe charge distribution for each PMT
        '''
        debug = False
        f = open(self.ratPMTChargeFile,'r')
        Qpmt  = {} # {PMT: [ (charge), (charge_prob) ]
        contents = f.read().split('\n')
        f.close()
        words = ['index:', 'charge:', 'charge_prob:']
        locs = {w:[] for w in words}
        for word in words:
            for i,line in enumerate(contents):
                if word in line:
                    locs[word].append( i )
        if debug : print('pmt.readCharge locs',locs)
        for i,jdx in enumerate(locs['index:']):
            jch = locs['charge:'][i]
            jpb = locs['charge_prob:'][i]
            line = contents[jdx]
            wdx = line[:-1].split()[1].replace('"','')
            line = contents[jch]
            chtuple = self.fillList(line,'charge:')
            line = contents[jpb]
            pbtuple = self.fillList(line,'charge_prob:')
            Qpmt[wdx] = chtuple,pbtuple
        if debug : print('pmt.readCharge Qpmt',Qpmt)

        if debug :
            for wdx in Qpmt:
                ch,pb = Qpmt[wdx]
                for i in range(5):
                    print('pmt.readCharge',wdx,' i,ch,pb',i,ch[i],pb[i],)
        return Qpmt
    def fillList(self,line,words):
        debug = False
        if debug : print('pmt.fillList line[:50]',line[:50],'line[-20:]',line[-10:])
        L = line[:-1].replace(words,'').replace('[','').replace(']','')
        if debug : print('pmt.fillList L[:40]',L[:40],'L[-40:]',L[-40:])

        LIST = []
        for x in L.split(','):
            if x and not x.isspace()  : LIST.append(float(x))
        return LIST
    def plotCharge(self,Qpmt):
        '''
        plot the prob as a function of charge prob
        '''
        ls = ['solid','dashed']
        for i,pmtName in enumerate(Qpmt):
            if pmtName != 'r7723':
                pmt = pmtName.replace('r7723_','')
                x,y = Qpmt[pmtName]
                x = numpy.array( x )
                y = numpy.array( y )
                y = y/numpy.sum(y)
                ymax = y.max()
                imax = y.argmax()
                xmax = x[imax]
                il,ir = imax,imax
                xl,xr = None,None
                for i in range(len(x)):
                    il = imax - i
                    ir = imax + i
                    if xl is None and y[il]<ymax/2. : xl = il
                    if xr is None and y[ir]<ymax/2. : xr = ir
                fwhm = xr-xl
                cxmax = ' ({:.0f},{:.0f})'.format(xmax,fwhm)
                plt.plot(x,y,label=pmt+cxmax,linestyle=ls[i%2])
        plt.grid()
        plt.legend(loc='best',ncol=1)
        plt.xlabel('Charge')
        plt.ylabel('Probability')
        plt.title('Single photoelectron charge spectra (max,fwhm)')
        pdf = self.Figures + 'simulated_spe_charge_spectra.pdf'
        plt.savefig(pdf)
        print('pmt.plotCharge Created',pdf)
        plt.show()
        return
    def main(self):
        Qpmt = self.readCharge()
        self.plotCharge(Qpmt)
        return
if __name__ == '__main__' :
    P = pmt()
    P.main()

        
