#!/usr/bin/env python
'''
analyze table of data and CM coincidence rates for 
EVEN/ODD, RED/BLUE, INNER/OUTER groupings
20220519
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt

class evenodd():
    def __init__(self):
        self.sources   = ['DATA','MC']
        self.groupings = ['RED-BLUE','EVEN-ODD',  'INNER-OUTER']
 # A.B, A.!B, !A.B
        self.coinc = {}
        self.coinc['MC'] =  {'EVEN-ODD': [2487, 5678, 6706],  
                    'RED-BLUE':[2627,6290,5087],
                    'INNER-OUTER':[2319, 8411, 4809]}
			
        self.coinc['DATA'] =  {'EVEN-ODD': [4176, 4498, 4960],
                    'RED-BLUE':[800, 4468, 5196],
                        'INNER-OUTER':[796, 5366, 5507] }

        print('evenodd.__init__ Totals in each grouping')
        for src in self.sources:
            for group in self.groupings:
                print(src,group,sum(self.coinc[src][group]))

        self.figDir = 'EVENODD_FIGURES/'
			
        return
    def fraction(self,L):
        '''
        return  fractions and uncertainties for input list L of unique selections
        '''
        tot = sum(L)
        f,df = [],[]
        for x in L:
            F = x/tot
            f.append( F ) 
            df.append( math.sqrt(F*(1.-F)/tot) )
        return f,df
    def ratios(self,L):
        '''
        return ratios and uncertainties
        '''
        L0 = L[0]
        r,dr = [],[]
        for x in L[1:]:
            R = x/L0
            dR = R*math.sqrt(1./x + 1./L0)
            r.append( R )
            dr.append( dR )
        return r,dr
    def plot(self,Fracs,dFracs):
        '''
        first plot data and MC fractions for all groupings as bar chart
        second plot as pie chart
        '''
        xlab = []
        align = ['center','center']
        width = [.5, .5]
        color = ['red','blue']
        dx = 0.0
        xstep = 1. #2*dx
        xtext = []
        for i in range(3):
            xlab.extend(  [ r'$A\cdot B$', r'$A\cdot\bar{B}$', r'$\bar{A}\cdot B$'] )
        for I,src in enumerate(self.sources):
            x = float(I)*dx
            xticks = []
            for group in self.groupings:
                f = Fracs[src][group]
                df=dFracs[src][group]
                X,Y,dY = [],[],[]
                for pair in zip(f,df):
                    X.append( x )
                    x += xstep
                    Y.append(pair[0])
                    dY.append(pair[1])
                plt.errorbar(X,Y,dY,label=src+' '+group,marker='.',drawstyle='steps-mid')
                xticks.extend( X )
                xtext.append( 0.5*(min(X)+max(X)) )
        plt.xticks(xticks,labels=xlab)
        for xt,group in zip(xtext,self.groupings):
            plt.text(xt,0.01,group,horizontalalignment='center')
        plt.legend(ncol=2)
        plt.ylim([0.,.7])
        plt.grid(linestyle=':')
        plt.title('Compare fractions within each group')
        pdf = self.figDir + 'fractions_within_groups.pdf'
        plt.savefig(pdf)
        print('evenodd.plot Wrote',pdf)
        plt.show()

        fig, ax = plt.subplots(nrows=2,ncols=3)
        
        for I,src in enumerate(self.sources):
            x = float(I)*dx
            xticks = []
            for J,group in enumerate(self.groupings):
                f = Fracs[src][group]
                df=dFracs[src][group]
                ax[I,J].pie(f,labels= [ r'$A\cdot B$', r'$A\cdot\bar{B}$', r'$\bar{A}\cdot B$'],labeldistance=0.5)
                ax[I,J].set_title(group)
                if J==0: ax[I,J].set_ylabel(src)
        pdf = self.figDir + 'fractions_pie.pdf'
        plt.savefig(pdf)
        print('evenodd.plot Wrote',pdf)
        plt.show()
        

        
        return
        
    def main(self):
        Fracs,dFracs = {},{}
        Ratios,dRatios = {},{}
        for src in self.sources:
            Fracs[src],dFracs[src] = {},{}
            Ratios[src],dRatios[src] = {},{}
            for group in self.coinc[src]:
                f,df = self.fraction(self.coinc[src][group])
                Fracs[src][group],dFracs[src][group] = f,df
                #print('Fractions',src,group,[p for p in zip(f,df)])

                r,dr = self.ratios(self.coinc[src][group])
                Ratios[src][group],dRatios[src][group] = r,dr
                #print('Ratios',src,group,[p for p in zip(r,dr)])
                    
        #print('Fracs',Fracs,'\ndFracs',dFracs)
        self.plot(Fracs,dFracs)
        return
if __name__ == '__main__' :
    P = evenodd()
    P.main()
    
