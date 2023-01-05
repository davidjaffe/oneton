#!/usr/bin/env python
'''
plot results of linear fit from arXiv:2006.00173 
20221230
20230103 add results from 1508.07029
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

        # 1508.07029 contents of 4
        self.BNLabsLY = {'0.4%': (19.9, 2.3),
                         '1%'  : (108.9, 10.9),
                        'LABPPO':(9156,917)}
        
        
        # contents of table 1 absLY[concentration] = (central value, unc)
        self.absLY = {'1%': (234., 30.),
                      '5%': (770., 72.),
                      '10%':(1357,125.),
                      'LABPPO':(11076, 1004)}
        self.conc = {'0.4%':0.4,
                     '1%': 1., 
                     '5%':5., 
                     '10%':10., 
                     'LABPPO':100.}
        # calculate LY relative to LAPPPO assuming uncertainties are independent (add in quad)
        self.relLY = {}
        keyden = 'LABPPO'
        den,eden = self.absLY[keyden]
        for key in self.absLY:
            if key!=keyden:
                num,enum = self.absLY[key]
                r = num/den
                er = r*math.sqrt(enum*enum/num/num + eden/eden/den/den)
                self.relLY[key] = (r,er)
                print('linearLY.__init__ LY('+key+') relative to keyden {:.4f}({:.4f})'.format(r,er))

        self.relBNLLY = {}
        for key in self.BNLabsLY:
            if key!=keyden:
                num,enum = self.BNLabsLY[key]
                r = num/den
                er = r*math.sqrt(enum*enum/num/num + eden/eden/den/den)
                self.relBNLLY[key] = (r,er)
                print('linearLY.__init__ BNL LY('+key+') relative to keyden {:.4f}({:.4f})'.format(r,er))

                
        
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

        # BNL data
        BNLLY, dBNLLY, BNLx = [],[],[]
        for c in self.BNLabsLY:
            y,dy = self.BNLabsLY[c]
            BNLLY.append(y)
            dBNLLY.append(dy)
            BNLx.append( self.conc[c] )
        BNLLY, dBNLLY, BNLx = self.toNPA(BNLLY), self.toNPA(dBNLLY), self.toNPA(BNLx)

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

        for xma,rng in zip([101.,10.1],['Xfull','Xpartial']):
            for yscale in ['linear','log']:
                yw = 'Y'+yscale
                for xscale in ['linear','log']:
                    words = rng+yw+'X'+xscale 

                    ### plot data
                    plt.errorbar(x,LY,fmt='o',yerr=dLY,color='black',label='LBL')

                    ### plot BNL data
                    plt.errorbar(BNLx,BNLLY,fmt='o',yerr=dBNLLY,color='blue',label='BNL')

                    ### draw family of fit results as grey band
                    for pair in p:
                        m,b = pair
                        YY = m*X+b
                        plt.plot(X,YY,color='grey',linestyle='solid',alpha=0.1)
                    ### main fit result
                    plt.plot(X,Y,'r-')

                    plt.title('absolute LY arXiv:2006.00173 table 1')
                    plt.grid()
                    plt.legend(loc='best')
                    
                    plt.xlabel('% Concentration')
                    plt.ylabel('Light yield (ph/MeV/% conc)')
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    if xscale=='linear': plt.xlim(-1.,xma)
                    if yscale=='linear' and xma<100. : plt.ylim(0.,2000.) 
                    png = 'LINEARLY_FIGURES/linearLY_arXiv2006.00173_'+words+'.png'
                    plt.savefig(png)
                    print('linearLY.main Wrote',png)
                    plt.show()

        
        return
if __name__ == '__main__' :
    P = linearLY()
    P.main()
    
