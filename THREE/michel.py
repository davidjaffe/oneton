#!/usr/bin/env python
'''
plot distribution of michel electron energy
units are MeV
20220125
'''
import numpy

import os
#import datetime
import sys
import math
import matplotlib.pyplot as plt
import mpl_interface

import Logger

class michel():
    def __init__(self):
        em = self.emass = 0.510999
        mm = self.mumass= 105.658
        pmax = (mm*mm - em*em)/2/mm
        self.Eemax = math.sqrt(pmax*pmax + em*em)
        print('michel.__init__ Eemax {:.2f} MeV'.format(self.Eemax))
        return
    def edist(self,Ee):
        '''
        wikipedia sez positron energy distribution, integrated over polar angle is dG/dx = x^2*(3-2x) where x = Ee/Eemax
        same expression from RPP2021 equation 57.3
        this is for muon decay in vacuum
        '''
        if Ee > self.Eemax : return 0.
        x = Ee/self.Eemax
        return x*x*(3.-2.*x)
    def gdist(self,Ee):
        '''
        muon decay width 
        I(x) = 1. - 8x - 12x^2*ln(x) + 8x^3 - x^4 for
        x = 2*Ee/(m_muon*c^2)
        '''
        x = 2*Ee/self.mumass
        if x==0. or x>1. : return 0.
        f = 1. - 8.*x - 12*x*x*math.log(x) + 8*x*x*x - x*x*x*x
        return f
    def main(self):
        '''

        '''
        fn = 'CRY_DATA/out_20210823_1.dat'


        Ee = numpy.linspace(0.,self.Eemax+1.,1000)
        y = numpy.array( [self.edist(x) for x in Ee] )
        z = numpy.array( [self.gdist(x) for x in Ee] )

        for name,Y in zip(['edist'],[y]):
        
            plt.plot(Ee,Y,'-r')
            mean = numpy.average(Ee,weights=Y)
            title = 'mean={:.2f} MeV'.format(mean)
            plt.grid()
            plt.xlabel('Electron energy (MeV)')
            plt.ylabel('Michel electron energy distribution')
            plt.title(title)
            savefig = True
            if savefig : 
                pdf = 'michel_electron_energy_spectrum.pdf'
                print('michel.main Saved figure to',pdf)
                plt.savefig(pdf)
            plt.show()


        return
if __name__ == '__main__' :
    ca = michel()
    ca.main()
