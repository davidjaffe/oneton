#!/usr/bin/env python
'''
analyze table of data and CM coincidence rates for 
EVEN/ODD, RED/BLUE, INNER/OUTER groupings
20220519

add analysis of CF factors for different muon generators
from WbLS_0614_2022.pptx page 4
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


### added 20220614
        self.CF = {}
        self.CF['CRY'] = {'R/B':[0.89, 0.85, 0.61, 0.81, 1.85, 1.52, 0.87, 0.76],
                          'E/O':[1.20, 0.95, 0.69, 1.03, 1.97, 1.68, 1.02, 0.95],
                          'I/O':[1.09, 0.75, 0.60, 0.91, 1.77, 1.34, 0.88, 0.79]}
        self.CF['FLAT']= {'R/B':[0.83, 0.81, 0.60, 0.82, 1.97, 1.63, 0.84, 0.72],
                          'E/O':[1.14, 0.91, 0.67, 1.04, 2.13, 1.89, 0.97, 0.90],
                          'I/O':[1.04, 0.71, 0.60, 0.93, 1.92, 1.48, 0.87, 0.74]}
        self.CF['TGM'] = {'xxx':[1.15, 0.74, 0.42, 0.95, 1.63, 1.43, 1.00, 0.78]}
        self.stopmuModes = ['CRY','FLAT']
        self.PMT = [x for x in range(8)]

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
    def plotCF(self):
        '''
        plot contents of self.CF

        First create bar chart similar to that in original pptx

        Second histograms of CF distributions per PMT

        Third bar chart of ratios of CRY to FLAT for stopping muons

        '''
        fig, ax = plt.subplots()
        N = 0
        for mode in self.CF:
            for sel in self.CF[mode]: N += 1
        print('evenodd.plotCF N',N)
        
        binsize = 1.0
        abit = 0.15
        barsize = binsize - abit
        width = barsize/float(N)
        space = abit/4.
        i,delta = 0,0.
        bar = {}
        for mode in self.CF:
            for sel in self.CF[mode]:
                y = self.CF[mode][sel]
                x = [q+i*width -barsize/2.+delta for q in self.PMT]
                label = sel + '_' + mode
                bar[i] = ax.bar(x, y, width, label=label)
                i += 1
            delta += space
        ax.set_ylabel('CF')
        ax.set_xlabel('Signal PMT')
        ax.set_title('CF for stopped and through-going muons')
        ax.set_xticks(self.PMT)
        ax.legend()
        ax.grid(axis='y')
        pdf = self.figDir + 'CF_all.pdf'
        plt.savefig(pdf)
        print('evenodd.plotCF Wrote',pdf)
        plt.show()

        ########## distributions of CF for each PMT for stopped muons, compared to thru-going muons
        cfPMT = {i:[] for i in self.PMT}
        nbins = 8
        Title = 'CF for '
        for mode in self.stopmuModes:
            Title += mode + ' '
            for sel in self.CF[mode]:
                for i,v in enumerate(self.CF[mode][sel]):
                    cfPMT[i].append(v)
        Title += 'compared to thru-going muon(red line)'
        fig, ax = plt.subplots(nrows=2,ncols=4,sharey='all')
        chisqr = 0.
        for i in self.PMT:
            j,k = i//4,i%4
            ax[j,k].hist(cfPMT[i],bins=nbins,linewidth=1.0)
            cfTGM = self.CF['TGM']['xxx'][i]
            xmi, xma = min( min(cfPMT[i]),cfTGM ), max( max(cfPMT[i]),cfTGM ) 
            delta = (xma - xmi)/float(nbins-1)
            minmax = [xmi-1.*delta, xma+1.*delta]
            ax[j,k].set(xlim=minmax)
            Q = numpy.array(cfPMT[i])
            mean,std = numpy.mean(Q),numpy.std(Q)
            con = (cfTGM-mean)/std
            chisqr += con*con
            title = 'mean {:.2f}\nstd.dev. {:.2f} n$\sigma$ {:.2f}'.format(mean,std,con)
            ax[j,k].set_title(title,fontsize=7,y=1.0-0.11)
            ax[j,k].set_xlabel('S'+str(i),labelpad=-125,loc='right')
            ax[j,k].axline( (cfTGM,0.), (cfTGM,+1.), color='red' )
        Title += ' $\chi^2/$ndf {:.0f}/{:}'.format(chisqr,(len(self.PMT)))
        fig.suptitle(Title)
        pdf = self.figDir + 'CF_dists_stopmu.pdf'
        plt.savefig(pdf)
        print('evenodd.plotCF Wrote',pdf)
        plt.show()

        ################
        fig, ax = plt.subplots()
        flat = self.CF['FLAT']
        cry  = self.CF['CRY']
        binsize, abit = 1,0.15
        barsize = binsize-abit
        i = 0
        width = barsize/3.
        bar = {}
        for i,sel in enumerate(flat):
            f = numpy.array(flat[sel])
            c = numpy.array(cry[sel])
            y = f/c
            x = [q+i*width-barsize/2. for q in self.PMT]
            label = sel
            bar[i] = ax.bar(x,y,width,label=label)
        ax.set_ylim( (0.8,1.2) )
        ax.set_ylabel('CF(FLAT)/CF(CRY)')
        ax.set_xlabel('Signal PMT')
        ax.set_title('CF ratios for stopped muon selections')
        ax.legend()
        ax.grid(axis='y')
        pdf = self.figDir + 'CF_FLAT_over_CF_CRY.pdf'
        plt.savefig(pdf)
        print('evenodd.plotCF Wrote',pdf)
        plt.show()
            
        
        
        
        return
    def main(self):
        EvenOddCoincRates = False
        if EvenOddCoincRates : 
        
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

        CFtable = True
        if CFtable :
            self.plotCF()
            
        return
if __name__ == '__main__' :
    P = evenodd()
    P.main()
    
