#!/usr/bin/env python
'''
fit 1d hists produced by Rong to measure muon lifetime 
20211221
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt

class stability():
    def __init__(self):
        self.input = 'stability.txt'
        self.Figures = 'STABILITY_FIGURES/'
        return
    def readInput(self):
        '''
        return dict {PMT: [ (Nentries, mean, std.dev., period), ...] }
        '''
        f = open(self.input,'r')
        print('stability.readInput Opened',self.input)
        results = {}
        for line in f:
            if line[0]!='#':
                s = line.split()
                PMT = s[0]
                Nentries = int(s[1])
                mean = float(s[2])
                stddev=float(s[3])
                period=s[4]
                if PMT not in results : results[PMT] = []
                results[PMT].append( (Nentries, mean, stddev, period) )
        f.close()
        return results
    def plotResults(self,results,liquid='L'):
        '''
        plot mean/PMT vs period
        liquid = 'L' ==> WbLS
        liquid = 'W' ==> Water
        liquid = 'A' ==> All
        '''
        fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2)
        ax0.set_title('mean NPE vs period')
        period = results['S0'][0][3]
        ax1.set_title('Relative to '+period)
        ax0.grid()

        ax1.grid()
        for PMT in results:
            x,y,dy = [],[],[]
            for tup in results[PMT]:
                N,mean,stddev,period = tup
                if liquid=='A' or liquid in period:
                    dmean = stddev/math.sqrt(N)
                    x.append( period )
                    y.append( mean )
                    dy.append( dmean )
            Y = numpy.array(y)
            DY= numpy.array(dy)
            wtav = numpy.average(Y,weights=1./DY/DY)
            nY = Y/wtav
            nDY= DY/wtav
            fY = Y/Y[0]
            fDY= DY/Y[0]
            ax0.errorbar(x,Y,yerr=DY,label=PMT,marker='o')
#            ax1.errorbar(x,nY,yerr=nDY,label=PMT,marker='o')
            ax1.errorbar(x,fY,yerr=fDY,label=PMT,marker='o')
        y1,y2 = ax0.get_ylim()
        ax0.set_ylim(0.,y2)
        ax1.legend(loc='best',ncol=2,title='PMT')
        
        ax1.set_ylim(0.5,1.15)
        x1,x2 = ax1.get_xlim()
        for d,ls in zip([0.1,0.05],['dashed','dotted']):
            for h in [1.-d,1.+d]:
                ax1.plot([x1,x2], [h,h], linestyle=ls,color='black')

        zoomY = False
        if zoomY : ax1.set_ylim(0.89,1.06) 
        # plt.show()
        pdf = self.Figures + 'NPE_vs_period_' + liquid + '.pdf'
        plt.savefig(pdf)
        print('stability.plotResults Wrote',pdf)
        return
    def main(self):
        results = self.readInput()
        # print('stability.main results',results)
        self.plotResults(results)
        return
if __name__ == '__main__' :
    S = stability()
    S.main()

