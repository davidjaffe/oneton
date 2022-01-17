#!/usr/bin/env python
'''
analyze PMT calibration file
20220117
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

class pmtcal_anal():
    def __init__(self):
        self.figDir = 'PMTCAL_FIGURES/'
        
        return
    def main(self):
        '''

        '''
        fn = 'PMTcalibration.txt'
        # format is
        # firstrun last run
        # c0 c1 ... c7
        # repeats
        f = open(fn,'r')
        cal = {} # {firstrun : [c0,c1,...]}
        rr  = {} # {firstrun : lastrun}
        x,dx = [],[]
        y = [[] for i in range(8)] #[ [c0_r0, c1_r1,...], [c1_r0, c1_r1, ...], ... [c7_r0, c7_r1, ...] ]
        for iline,line in enumerate(f):
            s = line.split()
            if iline%2==0 :
                r1,r2 = int(s[0]),int(s[1])
                cal[r1] = []
                rr[r1]  = r2
                x.append( 0.5*float(r1+r2) )
                dx.append( 0.5*float(r2-r1) )
            else:
                c = [float(q) for q in s]
                cal[r1] = c
                for i,z in enumerate(c):
                    y[i].append( z )
        f.close()
        print('pmtcal_anal.main Read',len(cal),'sets of PMT constants from',fn)
        rmi = min(rr.keys())
        rma = max(rr.values())
        print('pmtcal_anal.main First,last run',rmi,rma)

        colors = ['b','g','purple','c','m','k', 'r', 'gold']
        
        dy = [0. for q in x]
        ydelta = 5. 
        ymi,yma = 130., 210.+10.
        dataSets = self.runRanges()
        runLimits = {'All': [], 'Water': [], 'WbLS': []}
        rLkeys = ['All','Water','WbLS'] 
        for selection in rLkeys:
            for name in dataSets:
                r1,r2 = dataSets[name][:2]
                if ('W' in name and 'Water' in selection) or ('L' in name and 'WbLS' in selection) or selection=='All':
                    if len(runLimits[selection])==0 : runLimits[selection] = [r1,r2]
                    rmi = min(r1,runLimits[selection][0])
                    rma = max(r2,runLimits[selection][1])
                    runLimits[selection] = [rmi,rma]
        print('pmtcal_anal.main runLimits',runLimits)
        for selection in runLimits:
            rmi,rma = runLimits[selection]
            for iz,z in enumerate(y):
                Y = numpy.array(z)
                color = colors[iz]
                plt.errorbar(x,Y,xerr=dx,yerr=dy,label='S'+str(iz),color=color)
            plt.xlabel('Run (error bars show run range)')
            plt.ylabel('Gain calibration constant')
            plt.title(selection)
            plt.ylim(ymi,yma)
            plt.xlim(rmi-10.,rma+10.)
            iname = 0
            shades = ['green','orange']
            for name in dataSets:   ### indicate different datasets by alternate shading
                iname += 1
                x1,x2 = dataSets[name][:2]
                xave = 0.5*float(x1+x2)

                shade = shades[iname%2]
                plt.axvspan(x1,x2,color=shade,alpha=0.2)
                plt.text(xave, yma+ydelta*float(iname%2), name, rotation='vertical', color=shade, horizontalalignment='left')
                
            plt.grid()
            plt.legend(ncol=2)
            pdf = self.figDir + 'PMT_calibration_'+ selection + '.pdf'
            plt.savefig(pdf)
            print('pmtcal_anal.main wrote',pdf)
            plt.show()
            
        for selection in runLimits:
            calVal = {i:[] for i in range(8)}
            xmi,xma = 1.e6,-1e6
            for firstRun in cal:
                if self.validRun(selection,firstRun,dataSets) and self.validRun(selection,rr[firstRun],dataSets):
                    for i,z in enumerate(cal[firstRun]):
                        calVal[i].append(z)
                        xmi = min(xmi,z)
                        xma = max(xma,z)
            ymax = 0.
            for i in sorted(calVal.keys()):
                H = numpy.array(calVal[i])
                mean,sd = H.mean(),H.std()
                if mean!=0.: label = 'S{0} {1:.1f} {2:.2f}%'.format(i,mean,sd/
                                                                        mean*100.)
                histo,bins,patches = plt.hist(calVal[i],color=colors[i],label=label,histtype='step')
                ymax = max(ymax,max(histo))
            plt.title(selection)
            plt.legend()
            plt.grid()
            plt.ylim(0.,1.1*ymax)
            pdf = self.figDir + 'PMT_gains_' + selection + '.pdf'
            plt.savefig(pdf)
            print('pmtcal_anal.main wrote',pdf)
            plt.show()
    
        return
    def validRun(self,selection,run,dataSets):
        '''
        return True if run is in dataSets based on selection
        '''
        if selection=='All' : return True
        if selection=='Water' or selection=='WbLS':
            for name in dataSets:
                r1,r2 = dataSets[name][:2]
                if (selection=='Water' and name[0]=='W') or (selection=='WbLS' and name[0]=='L') :
                    if r1<=run<=r2 : return True
        return False            
    def runRanges(self):
        '''
        return definition of dataSets
        dataSets = {} # = {setName : [firstrun, lastrun, firstdate, lastdate, comment]
        '''
        dataSets = {}
        dataSets['W00']=[20043,20696, '1/2/2018', '1/10/2018', 'Start at run 20043']
        dataSets['W01']=[20710,22061, '1/13/2018','1/19/2018','']
        dataSets['W02']=[22062,24468, '2/6/2018','3/25/2018', 'Optical coupling changed!']
        dataSets['W03']=[24470,28656, '3/26/2018','6/12/2018', '']
        dataSets['W04']=[28659,28789, '6/12/2018','6/13/2018', 'Circulation off']
        dataSets['W05']=[28796,28846, '6/14/2018','6/14/2018', 'Circulation on with different speed']
        dataSets['W06']=[28850,29431, '6/15/2018','6/15/2018', 'Circulation on with different speed']
        dataSets['W07']=[29432,29550, '6/25/2018','6/25/2018', 'Circulation on with different speed']

        dataSets['L00']=[32494,35741, '8/6/2018','8/27/2018', 'Start WbLS data','taking, with circulation on']
        dataSets['L01']=[35742,36276, '8/27/2018','7/7/2018', ' circulation off' ]
        dataSets['L02']=[36277,36297, '10/10/2018','10/12/2018', 'Water added' ]  
        dataSets['L03']=[36447,36815, '10/13/2018','10/26/2018', 'New dataset with hodo trigger only']
        dataSets['L04']=[36881,38614, '11/2/2018','1/11/2019', 'Circulation on, hodo trigger only']  
        dataSets['L05']=[38868,38971, '1/19/2019','1/28/2019', 'Circulation off, hodo trigger only']
        return dataSets
if __name__ == '__main__' :
    ca = pmtcal_anal()
    ca.main()
