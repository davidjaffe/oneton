#!/usr/bin/env python
'''
plot results of linear fit from arXiv:2006.00173 

20221230
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit

class linearLY():
    def __init__(self):

        # contents of table 1 absLY[concentration] = (central value, unc)
        self.absLY = {'1%': (234., 30.),
                      '5%': (770., 72.),
                      '10%':(1357,125.),
                      'LABPPO':(11076, 1004)}
        self.conc = {'1%': 1., 
                      '5%':5., 
                      '10%':10., 
                      'LABPPO':100.}
        
        # section 3.1.3 reported results of linear fit to 3 wbls points
        self.lfit = {'slope': (127.9,17.0),
                     'intercept': (108.3, 51.0)}
        return

    def toNPA(self,a):
        return numpy.array(a)
    def linear(self, x, m, b):
        return m*x + b
    def main(self):

        LY,dLY,x = [],[],[]
        for c in self.absLY:
            y,dy = self.absLY[c]
            LY.append(y)
            dLY.append(dy)
            x.append( self.conc[c] )
        LY, dLY, x = self.toNPA(LY), self.toNPA(dLY), self.toNPA(x)

        ## do linear fit to the first 3 points
        param, cov = curve_fit(self.linear, x[:3], LY[:3], sigma=dLY[:3])
        print('linearLY.main results of linear fit, param',param,'cov',cov)

        ## generate a family of fit parameters consistent with
        ## the fit results
        size = 500
        p = numpy.random.multivariate_normal(param,cov,size=size)

            
        
        X = numpy.linspace(min(x)/10.,max(x),1000)
        m, b = self.lfit['slope'][0],self.lfit['intercept'][0]
        Y = m*X+b

        for xma in [101.,10.1]:
            for yscale in ['linear','log']:
                for xscale in ['linear','log']:

                    ### plot data
                    plt.errorbar(x,LY,fmt='o',yerr=dLY,color='black')

                    ### draw family of fit results as grey band
                    for pair in p:
                        m,b = pair
                        YY = m*X+b
                        plt.plot(X,YY,color='grey',linestyle='solid',alpha=0.1)
                    ### main fit result
                    plt.plot(X,Y,'r-')

                    plt.title('absolute LY arXiv:2006.00173 table 1')
                    plt.grid()

                    plt.xlabel('% Concentration')
                    plt.ylabel('Light yield (ph/MeV/% conc)')
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    if xscale=='linear': plt.xlim(-1.,xma)
                    if yscale=='linear' and xma<100. : plt.ylim(0.,2000.) 
    #        png = 'LINEARLY_FIGURES/linearLY_arXiv2110.13222.png'
    #        plt.savefig(png)
    #        print('linearLY.main Wrote',png)
                    plt.show()

        
        return
if __name__ == '__main__' :
    P = linearLY()
    P.main()
    
