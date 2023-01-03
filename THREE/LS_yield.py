#!/usr/bin/env python
'''
20230103 compare, plot, etc. a number of results on relative and absolute light yield

Ref1 - 2006.00173 : WbLS aLY (LAB-based) and LABPPO(2 g/L) aLY
Ref2 - 1508.07029 :  WbLS aLY (PC-based) and  LABPPO(2 g/L) aLY
Ref3 - NIM A 701 (2013) 133–144 : Relative LY of LAB+PPO, DIN+PPO to PC+PPO
Ref4 - 2205.15046 : LAB+4g/L PPO, DIN+4g/L PPO, p-Xylene LY, EJ301 relative to anthracene
Ref5 - 2003.10491 : LAB + g/L PPO relative to EJ301

'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit

class LS_yield():
    def __init__(self):

        # 1508.07029 WbLS LY contents of table 4
        self.BNLabsLY = {'0.4%': (19.9, 2.3),
                         '1%'  : (108.9, 10.9),
                        'LABPPO':(9156,917)}
        
        
        # 2006.00173 contents of table 1 absLY[concentration] = (central value, unc)
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
                print('LS_yield.__init__ LY('+key+') relative to keyden {:.4f}({:.4f})'.format(r,er))

        self.relBNLLY = {}
        for key in self.BNLabsLY:
            if key!=keyden:
                num,enum = self.BNLabsLY[key]
                r = num/den
                er = r*math.sqrt(enum*enum/num/num + eden/eden/den/den)
                self.relBNLLY[key] = (r,er)
                print('LS_yield.__init__ BNL LY('+key+') relative to keyden {:.4f}({:.4f})'.format(r,er))

                
        
        # section 3.1.3 reported results of linear fit to 3 wbls points
        self.lfit = {'slope': (127.9,17.0),
                     'intercept': (108.3, 51.0)}


        # 2205.15046 = Ref4 dict[solvent] = [WLS, (rLY,drLY), LY]
        # results from table 2
        self.LY4 = {'EJ301'     : [ 'NA', (78.,0.), 13570],
                    'p-Xylene'  : [ '4 g/L PPO', (77.2,1.9), 13430],
                    'LAB'       : [ '4 g/L PPO', (60.2,1.4), 10470],
                    'DIN'       : [ '4 g/L PPO', (70.0,1.7), 12190],
                    'relativeTo': ['anthracene']}

        # NIM A701 (2013) 133
        # results from tables 1, 6 with uncertainties taken from
        #  statements in the text
        self.LY3 = {'LAB15'  : ['1.5 g/L PPO', (1.16, 1.16*.03)],
                    'LAB20'  : ['2.0 g/L PPO', (1.22, 1.22*.03)],
                    'LAB30'  : ['3.0 g/L PPO', (1.26, 1.26*.03)],
                    'LAB60'  : ['6.0 g/L PPO', (1.33, 1.33*.03)],
                    'DIN15'  : ['1.5 g/L PPO', (0.96, 0.96*.03)],
                    'DIN20'  : ['2.0 g/L PPO', (0.99, 0.99*.03)],
                    'DIN30'  : ['3.0 g/L PPO', (1.07, 1.07*.03)],
                    'DIN60'  : ['6.0 g/L PPO', (1.12, 1.12*.03)],
                    'relativeTo': ['PC 1.5 g/L PPO']}

        # 2003.10491
        # values and uncertainies estimated from figure 6
        self.LY5 = {'LAB05'  : ['0.5 g/L PPO', (0.37, 0.01)],
                    'LAB10'  : ['1.0 g/L PPO', (0.46, 0.01)],
                    'LAB20'  : ['2.0 g/L PPO', (0.64, 0.01)],
                    'LAB30'  : ['3.0 g/L PPO', (0.67, 0.01)],
                    'LAB40est':['4.0 g/L PPO', (0.70, 0.01)],
                    'LAB50'  : ['5.0 g/L PPO', (0.72, 0.01)],
                    'relativeTo': ['EJ301']}

        self.LY = {'Ref3' : self.LY3,
                   'Ref4' : self.LY4,
                   'Ref5' : self.LY5}

        
        
        return
    def getConc(self,words):
        '''
        return concentration in percent by parsing input words
        '''
        s = words.split(' ')[0]
        if s=='NA' : return 0.
        return float(s)
        
    def unpackRef(self,ref):
        '''
        return LS names, WLS conc, relative LY and reference LS
        given ref
        '''
        if ref not in self.LY :
            sys.exit('LS_yield.unpackRef ERROR Invalid ref '+ref)

        x,y,dy,z = [],[],[],[]
        relTo = None
        block = self.LY[ref]
        for LS in block:
            if LS=='relativeTo' :
                relTo = block[LS][0]
            else:
                conc = block[LS][0]
                v,dv = block[LS][1]
                if len(block[LS])> 2: z.append( block[LS][2] )
                x.append( self.getConc(conc) )
                y.append( v )
                dy.append( dv )
        if relTo is None: sys.exit('LS_yield.unpackRef ERROR No relativeTo for ref '+ ref)
        X = self.toNPA(x)
        Y = self.toNPA(y)
        dY= self.toNPA(dy)
        Z = self.toNPA(z)
        return relTo, X,Y,dY,Z 
    
    def toNPA(self,a):
        return numpy.array(a)
    def linear(self, x, m, b):
        return m*x + b
    def main(self):

        ref = 'Ref4'
        relTo, X,Y,dY,Z = self.unpackRef(ref)
        return

    def OLDmain(self):

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
        print('LS_yield.main results of linear fit, param',param,'cov',cov)

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
                    plt.errorbar(x,LY,fmt='o',yerr=dLY,color='black')

                    ### plot BNL data
                    plt.errorbar(BNLx,BNLLY,fmt='o',yerr=dBNLLY,color='blue')

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
                    png = 'LY_YIELD_FIGURES/LS_yield_arXiv2006.00173_'+words+'.png'
                    plt.savefig(png)
                    print('LS_yield.main Wrote',png)
                    plt.show()

        
        return
if __name__ == '__main__' :
    P = LS_yield()
    P.main()
    
