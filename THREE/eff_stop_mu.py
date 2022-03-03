#!/usr/bin/env python
'''
estimate stopped muon selection efficiency and bkgd with total NPE cut
uses figures decaymuon_npe.png and tgmuon_npe

20220303
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt

class eff_stop_mu():
    def __init__(self):

        self.Figures = 'EFF_STOP_MU_FIGURES/'
        if not os.path.exists(self.Figures): os.makedirs(self.Figures)

        self.upedge = [  4,   8,  12, 16, 20, 24, 28, 32, 36, 40]
        self.sm_cts = [860, 330, 265,215, 147, 128, 92, 125, 98, 100]
        self.sm_tot = 2621

        self.tg_cts = [270, 580, 945, 1050, 1100, 1360, 1560, 1670, 1780, 1890]
        self.tg_tot = 40869
        return
    def convert(self):
        '''
        convert counts/bin to cumulative counts for signal and background efficiencies

        plot results
        '''
        sm_cts = numpy.array( self.sm_cts )
        sm_cum = numpy.cumsum( sm_cts )
        sm_eff = sm_cum/self.sm_tot

        tg_cts = numpy.array( self.tg_cts )
        tg_cum = numpy.cumsum( tg_cts )
        tg_eff = tg_cum/self.tg_tot

        x = numpy.array( self.upedge )

        ## calculate effys if stop-muon effy is doubled
        sm_eff8 = numpy.interp(8., x, sm_eff)
        effnew = min(0.90, sm_eff8*2.)
        cut2 = numpy.interp(effnew,sm_eff,x)
        tg_eff8 = numpy.interp(8., x, tg_eff)
        tg_eff2 = numpy.interp(cut2, x, tg_eff)
        print('stopped muon effy for NPEsum>8 {:.3f} and throughgoing muon effy {:.3f}'.format(sm_eff8,tg_eff8))
        print('NPEsum> {:.1f} gives stopped muon effy of {:.3f} and throughgoing muon effy of {:.3f}'.format(cut2,effnew,tg_eff2))
        
        

        fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1)
        
        ax0.plot(x,sm_eff,'o-', label='Stopping muon efficiency')
        ax0.plot(x,tg_eff,'o-', label='Through-going muon efficiency')

        ax0.set_title('Efficiency vs total NPE cut')
        ax0.set_xlabel('Total NPE cut')
        ax0.set_ylabel('Efficiency')
        ax0.set_ylim(-0.02,1.0)
        ax0.legend(loc='best')
        ax0.grid()

        ax1.plot(tg_eff/tg_eff8, sm_eff/sm_eff8,'o-')
        ax1.set_title('Stopping muon vs through-going muon rates relative to NPE>8 cut')
        ax1.set_xlabel('Through-going muon relative rate')
        ax1.set_ylabel('Stopping muon relative rate')
        ax1.xaxis.set_ticks(numpy.arange(0.,14.+1.,1.))
        ax1.set_ylim(0.,2.2)
        ax1.grid()

        plt.tight_layout()

        pdf = self.Figures + 'sm_tg_selection_effy.pdf'

        plt.savefig(pdf)
        print('eff_stop_mu.convert Wrote',pdf)
        
        
        plt.show()

        
        
        
        return 
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
        self.convert()
        return
if __name__ == '__main__' :
    P = eff_stop_mu()
    P.main()

        
