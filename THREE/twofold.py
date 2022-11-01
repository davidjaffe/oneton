#!/usr/bin/env python
'''
determine efficiency of each signal PMT relative to MC 
using table of unique twofold coincidences from WbLS_0525_2022
20220525
'''
import numpy

import os

import datetime
import Logger
import pickle

import sys
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import pearsonr

import iminuit

class twofold():
    def __init__(self,debug=-1,nToy=0):
        
### usual directory structure to preserve output(log,figures,pickle,...) of each job
        now = datetime.datetime.now()
        self.now = now.strftime('%Y%m%dT%H%M%S')

        self.rootDir = 'JOBS_TWOFOLD/'
        parentDir = self.rootDir+self.now
        dirs = [parentDir]
        self.Figures = parentDir  + '/FIGURES/'
        self.logDir = parentDir
        dirs.append( self.Figures)
        dirs.append( self.logDir)

        for d in dirs:
            if not os.path.exists(d):
                os.makedirs(d)
                print('calibf.__init__ create directory',d)

        lf = self.logDir + '/logfile.log'
        sys.stdout = Logger.Logger(fn=lf)
        print('twofold.__init__ Output directed to stdout and',lf)

## process, report inputs
        self.debug = debug
        self.nToy = nToy
        self.toyMC = nToy > 0
        print('twofold.__init__ debug',self.debug,'nToy',self.nToy,'toyMC',self.toyMC)

## 20221011 definitions of SM selections
        self.pmtGroup = {'EVEN': [0,2,4,6], 'ODD':  [1,3,5,7],
                         'RED' : [0,2,5,7], 'BLUE': [1,3,4,6],
                         'INNER':[1,2,5,7], 'OUTER':[0,3,4,6] }
        self.pmtPairs = [ ['EVEN','ODD'], ['RED','BLUE'], ['INNER','OUTER'] ]
        self.grpPairs = []
        for p in self.pmtPairs: self.grpPairs.append( p[0]+'_'+p[1] )

        self.pmtCoinc = {}
        for group in self.pmtGroup:
            pmts = self.pmtGroup[group]
            self.pmtCoinc[group] = []
            for i,I in enumerate(pmts):
                for j in range(i+1,len(pmts)):
                    self.pmtCoinc[group].append( [pmts[i],pmts[j]] )

        self.uniqPairs = {}

        for pair in self.pmtPairs:
            A,B = pair
            print('twofold.__init__','pairs',A,self.pmtCoinc[A],'  ',B,self.pmtCoinc[B])
            aCoinc,bCoinc = [x for x in self.pmtCoinc[A]],[x for x in self.pmtCoinc[B]]
            for group in self.pmtGroup:
                if group not in pair:
                    for c in self.pmtCoinc[group]:
                        if self.debug > 1 : print(group,'c',c,'aCoinc',aCoinc,'bCoinc',bCoinc)
                        if c in aCoinc: aCoinc.remove(c)
                        if c in bCoinc: bCoinc.remove(c)
            self.uniqPairs[A] = aCoinc
            self.uniqPairs[B] = bCoinc
            print('twofold.__init__','unique pairs',A,self.uniqPairs[A],'  ',B,self.uniqPairs[B])


        
## table of unique twofold coincidences from WbLS_0525_2022
        self.sources   = ['DATA','MC']
        self.TwoFold = {}
        self.TwoFold['MC'] = numpy.array( (0,101836, 62216, 37352, 7560, 10556, 105196, 89348, # S0 in coinc S0-S7
                               101836, 0, 112644, 62832, 13356, 15372, 98756, 131880,
                               62216, 112644, 0, 89936, 16800, 19292, 66584, 101472,
                               37352, 62832, 89936, 0, 16828, 19320, 40768, 60452,
                               7560, 13356, 16800, 16828, 0, 66668, 11508, 15764,
                               10556, 15372, 19292, 19320, 66668, 0, 10304, 13216,
                               105196, 98756, 66584, 40768, 11508, 10304, 0, 115080,
                               89348, 131880, 101472, 60452, 15764, 13216, 115080,0)
                            ).reshape( (8,8) )
        self.TwoFold['DATA']=numpy.array( (0, 4258, 1892, 1940, 875, 720, 4156, 3739,
                               4258, 0, 2550, 2474, 988, 777, 3412, 4344,
                               1892, 2550, 0, 2098, 664, 522, 1538, 2185,
                               1940, 2474, 2098, 0, 999, 735, 1685, 2218,
                               875, 988, 664, 999, 0, 1355, 920, 1064,
                               720, 777, 522, 735, 1355, 0, 622, 723,
                               4156, 3412, 1538, 1685, 920, 622, 0, 4032,
                               3739, 4344, 2185, 2218, 1064, 723, 4032,0)
                               ).reshape( (8,8) )
        print('twofold.__init__ Consistency check')
        self.Npmt = self.nPMT = 8
        OK = True
        maximum = {}
        self.probTwoFold = {}
        for src in self.sources: 
            L = len(self.TwoFold[src][0])
            if L!=self.nPMT : print('inconsistent array length',L,'and nPMT',self.nPMT)
            m = 0
            for i in range(L):
                for j in range(L):
                    Nij = self.TwoFold[src][i][j]
                    Nji = self.TwoFold[src][j][i]
                    m = max(Nij, m)
                    if Nij!=Nji:
                        OK = False
                        print(src,'inconsistent i',i,'j',j,'Nij',Nij,'Nji',Nji)
            maximum[src] = m
            self.probTwoFold[src] = self.TwoFold[src]/(10.*m)
            if self.debug > 0 : print('twofold.__init__ self.probTwoFold[',src,']',self.probTwoFold[src])
        if OK : print('twofold.__init__ Consistency check OK')
        self.maximum = maximum

        X = self.TwoFold['MC']
        m = maximum[src]
        if self.debug > 0 : print('twofold.__init__ self.TwoFold[`MC`]',X)
        self.prob = X/m
        if self.debug > 0 : print('twofold.__init__ prob',self.prob)
        print('twofold.__init__ maximum',maximum)
## table of coincidences from tech note (and which prez?)                            
        self.cName = ['A&B ','A&~B', '~A&B']
        self.measCoinc = {'MC': {'EVEN_ODD': [2847, 5678, 6706], 'RED_BLUE': [2627,6290,5087], 'INNER_OUTER': [2319, 8411, 4809]},
                        'DATA': {'EVEN_ODD': [4176, 4498, 4960], 'RED_BLUE': [ 800,4468,5196], 'INNER_OUTER': [ 796, 5366,5507]} }

## calculate measured fractions
        self.fracCoinc = {src: {c:[] for c in self.grpPairs} for src in self.sources}
        for dmc in self.sources:
            for gname in self.grpPairs:
                v = numpy.array(self.measCoinc[dmc][gname])
                self.fracCoinc[dmc][gname] = v/sum(v)
        if self.debug > 0 :
            for dmc in self.sources:
                print('twofold.__init__ dmc,self.measCoinc',dmc,self.measCoinc)
                print('twofold.__inti__ dmc,self.fracCoinc',dmc,self.fracCoinc)
                
## report totals in all pairs in a group and for uniq pairs in each group
        for dmc in self.sources:
            for group in self.pmtCoinc:
                pmts = self.pmtCoinc[group]
                upmts= self.uniqPairs[group]
                tot,utot = 0,0
                for pair in pmts:
                    #print('pair',pair)
                    i,j = pair
                    cts = self.TwoFold[dmc][i][j]
                    tot += cts
                    if pair in upmts : utot += cts
                print('twofold.__init__ dmc,group,tot,utot,utot/max',dmc,group,tot,utot,'{:.3f}'.format(float(utot)/self.maximum[dmc]),'uniq pairs',upmts)        

        self.figDir = self.Figures #'TWOFOLD_FIGURES/'
			
        return
    def readPickle(self,timestamp):
        '''
        return input & fitresults given job specified by timestamp
        '''
        pickle_fn = self.rootDir + timestamp + '/' + timestamp + '.pickle'
        
        f = open(pickle_fn,'rb')
        Input = pickle.load(f)
        fitresults = pickle.load(f)
        f.close()
        print('twofold.readPickle Read',pickle_fn)
        
        return Input, fitresults
    def writePickle(self,Input, fitresults):
        '''
        write inputData & fitresults to pickle file
        order of dump must be matched in load according to https://stackoverflow.com/questions/20716812/saving-and-loading-multiple-objects-in-pickle-file
        '''
        timestamp = self.now
        pickle_fn = self.rootDir + timestamp + '/' + timestamp + '.pickle'


        f = open(pickle_fn, 'wb')
        pickle.dump(Input, f)
        pickle.dump(fitresults, f)
        f.close()

        print('twofold.writePickle Wrote',pickle_fn)
        return
    def oneExptOLD(self):
        '''
        return Random 'Data' for one generated experiment given input efficiencies and number of events

        OLD VERSION: TOO CLEVER, OUTPUT DATA IS NOT SYMMETRIC ABOUT DIAGONAL
        '''
        NN = self.Npmt*self.Npmt
        eff = self.input[:self.Npmt]
        epp = []
        for i in range(self.Npmt):
            for j in range(self.Npmt):
                epp.append( eff[i]*eff[j]*self.prob[i,j] )
        epp = numpy.array(epp)
        X = numpy.zeros( NN )
        Ones = numpy.ones( NN )
        N   = self.input[self.Npmt]
        for evt in numpy.arange(N):
            R = numpy.random.random( NN )
            X += (R<=epp)*Ones
        X = numpy.fix( X )
        X = X.reshape(self.Npmt,self.Npmt)
        return X
    def oneExpt(self):
        '''
        return Random 'Data' for one generated experiment given input efficiencies and number of events
        '''
        NN = self.Npmt*self.Npmt
        eff = self.input[:self.Npmt]
        epp = []
        for i in range(self.Npmt):
            for j in range(self.Npmt):
                epp.append( eff[i]*eff[j]*self.prob[i,j] )
        epp = numpy.array(epp).reshape(self.Npmt,self.Npmt)
        X = numpy.zeros( NN ).reshape(self.Npmt,self.Npmt)
        N   = self.input[self.Npmt]
        for evt in numpy.arange(N):
            for i in range(self.Npmt):
                for j in range(i+1,self.Npmt):
                    if numpy.random.random()<=epp[i,j]:
                        X[i,j] += 1.
                        X[j,i] += 1.
        X = numpy.fix( X )
        return X
    def oneTrueExpt(self):
        '''
        return TRUE 'Data' for one generated experiment given input efficiencies and number of events
        '''
        eff = self.input[:self.Npmt]
        N   = self.input[self.Npmt]
        X = numpy.zeros( self.Npmt*self.Npmt ).reshape( (self.Npmt,self.Npmt) )
        for i in range(self.Npmt):
            for j in range(self.Npmt):
                v = 0. 
                if i!=j : v = eff[i]*eff[j]*self.prob[i,j]*N
                X[i,j] += v
        X = numpy.fix( X )
        return X
    def chisqr(self,param):
        '''
        chisquared = sum_i=0,7 sum_j=i+1,7 (Cij - effi*effj*Pij*N)^2 / Cij
        C[i,j] = count rate in date for coincidence between PMTs i and j
        P[i,j] = probability of coincidence between PMTs i,j
        effi = efficiency of PMT i ( = param[0:8] )
        N = total number of data events = FIXED

        add array self.csTerms to pass individual terms for debugging
        '''
        eff = param[:self.Npmt]
        N   = self.FixedN
        cs = 0.
        self.csTerms = numpy.zeros( self.Npmt*self.Npmt ).reshape( (self.Npmt,self.Npmt) )
        C = self.inputData # TwoFold['DATA']
        P = self.prob
        for i in range(self.nPMT):
            for j in range(i+1,self.nPMT):
                num = (C[i,j] - eff[i]*eff[j]*P[i,j]*N)
                self.csTerms[i,j] = num
                cs += num*num/C[i,j]
        return cs
    def chisqrNINE(self,param):
        '''
        OLDER VERSION OF CHISQ WITH NINE PARAMETERS : 8 Efficiencies and N=total events

        chisquared = sum_i=0,7 sum_j=i+1,7 (Cij - effi*effj*Pij*N)^2 / Cij
        C[i,j] = count rate in date for coincidence between PMTs i and j
        P[i,j] = probability of coincidence between PMTs i,j
        effi = efficiency of PMT i ( = param[0:8] )
        N = total number of data events ( = param[8] )
        '''
        eff = param[0:8]
        N   = param[8]
        cs = 0.
        C = self.inputData # TwoFold['DATA']
        P = self.prob
        for i in range(self.nPMT):
            for j in range(i+1,self.nPMT):
                num = (C[i,j] - eff[i]*eff[j]*P[i,j]*N)
                cs += num*num/C[i,j]
        return cs
    def find(self,method='Powell',Nguess = 5000.):
        '''
        return best fit parameters from fit

        20220622 if method='migrad' or 'simplex', use iminuit interface
        '''
        ## initial guesses and bounds for parameters
        param = [1. for x in range(self.nPMT)]
        bounds = [(None,None) for x in range(self.nPMT)]

        self.FixedN = Nguess
#        bounds.append( (0.999*Nguess,1.001*Nguess) )  #(self.maximum['DATA'], 2.*Nguess) )
#        param.append(Nguess)
        chi2 = self.chisqr(param)
        if self.debug > 2 :
            with numpy.printoptions(precision=1,linewidth=200,suppress=True):
                print('\ntwofold.find Initial-values chisquare terms\n',numpy.array(self.csTerms))
        
        if self.debug > 0: 
            print('\ntwofold.find method',method)
            print('twofold.find chi2',chi2,'with input params',param)
        ## after 20220616T165320 don't give bounds for unbounded methods
        disp = self.debug > 1
        if method=='migrad' or method=='simplex':
            iminuit.Minuit.strategy = 2
            res = iminuit.minimize,minimize(self.chisqr, param, options={'disp':disp}) # migrad is default
        else:
            if self.bounded[method]: 
                res = minimize(self.chisqr, param, method=method,bounds=bounds)
            else:
                res = minimize(self.chisqr, param, method=method, options={'disp':disp})
            
        if self.debug > 1 :
            print('twofold.find res',res)
        if method=='migrad' or method=='simplex': ### deal with idiosyncracy of iminuit.minimize.minimize()
            res = res[-1]
        pout = res.get('x')        
        hess_inv = res.get('hess_inv') # 
        if type(hess_inv)!=numpy.ndarray : # then it is a useless object, so replace it
            if self.debug > 1 : print('twofold.find hess_inv is a useless object, so replace it zeros')
            hess_inv = numpy.zeros( self.nPMT*self.nPMT).reshape( (self.nPMT,self.nPMT) )
        success = res.get('success')
        message = res.get('message')
        chi2 = self.chisqr(pout)
        
        if self.debug > 2 :
            with numpy.printoptions(precision=1,linewidth=200,suppress=True):
                print('\ntwofold.find Best-fit chisquare terms\n',numpy.array(self.csTerms))
            with numpy.printoptions(precision=6,linewidth=300,suppress=True):
                print('\ntwofold.find hess_inv\n',hess_inv)
            if method=='migrad' or method=='simplex':
                print('twofold.find iminuit.Minuit.accurate',str(iminuit.Minuit.accurate))

        fitpar = 'twofold.find {0} chi2 {1:.1f} fitpar '.format(method,chi2)
        if self.debug > 1:
            print('twofold.find hess_inv\n',hess_inv)
        for i,p in enumerate(pout):
            punc = 0.
            if hess_inv is not None: punc = math.sqrt(max(0.,hess_inv[i,i]))
            fitpar += ' {0:.3f}({1:.3f})'.format(p,punc)

        if not success : fitpar += ' FIT FAILED. ' + message
        if self.debug > 0 or not success: print(fitpar)
        return pout
    def analyzeResults(self,fitresults):
        '''
        compare fitresults with input for toyMC
        return Chi2All and Chi2PMT
        '''
        print('\ntwofold.analyzeResults')
#        L = len(self.input)
        
### store analysis results
### chi2All = sum_all [(best fit - input)/standard deviation]^2
### chi2PMT = sum_PMT [(best fit - input)/standard deviation]^2
### all includes number of input events

        self.analysisResults = {} # [method: [chi2All, chi2PMT]]
        
        means,stddevs = {},{}
        CHI2ALL, CHI2PMT = {},{}
        for method in fitresults:
            L = len(fitresults[method][0])
            means[method],stddevs[method] = [],[]
            results = numpy.array( fitresults[method] )
            #print('method,results',method,results)
            Mresult = method + ' means'
            Sresult = method + 'stddev'
            chi2All, chi2PMT = 0., 0.
            for i in range(L):
                #if i==0: print('i,results[:,i]',i,results[:,i])
                m,s = numpy.mean( results[:,i] ), numpy.std(results[:,i] )
                c = (m-self.input[i])/s
                chi2All += c*c
                if i<self.nPMT : chi2PMT += c*c
                means[method].append( m )
                stddevs[method].append( s )
                Mresult += ' {:.3f}'.format(m)
                Sresult += ' {:.3f}'.format(s)
            print(Mresult)
            print(Sresult)
            CHI2ALL[method], CHI2PMT[method] = chi2All,chi2PMT
            print(method,'Chi2All {:.2f} Chi2PMT {:.2f}'.format(chi2All,chi2PMT))
            self.analysisResults[method] = [chi2All, chi2PMT]
        Winput = 'Input pararameters'
        for i in range(L):
            Winput += ' {:.3f}'.format(self.input[i])
        AA = self.input[:self.Npmt]
        m,p,s,span = numpy.mean(AA), numpy.prod(AA), numpy.std(AA), numpy.max(AA)-numpy.min(AA)
        Winput += ' mean {:.3f} prod {:.3f} stddev {:.3f} span {:3f}'.format(m,p,s,span)
        print(Winput)

        self.rankMethods()
        
        return CHI2ALL, CHI2PMT
    def rankMethods(self):
        '''
        produce list of methods ranked by smallest chi2PMT and chi2All
        '''
        A = self.analysisResults
        iAll,iPMT = 0,1
        s_PMT = sorted(A.items(), key=lambda A : A[1][iPMT])
        s_All = sorted(A.items(), key=lambda A : A[1][iAll])
        for i,word in zip( [iAll, iPMT], ['All', 'PMT']):
            print('\ntwofold.rankMethods Ranked by Chi2'+word)
            for X in s_PMT:
                method = X[0]
                Chi2   = X[1][i]
                print(method,'{:.3f}'.format(Chi2))
        return
    def plotResults(self,fitresults,method,timestamp):
        '''
        plot fitresults for minimization method for job specified by timestamp
        '''
        results = numpy.array( fitresults[method] )
        if method in self.analysisResults:
            chi2All, chi2PMT = self.analysisResults[method]
        else:
            chi2All, chi2PMT = -1., -1.
        Nevts = len(results)
        Input   = self.input
        #print('twofold.plotResults Input',Input)

        nrows, ncols = 2,4
        fig, ax = plt.subplots(nrows=nrows,ncols=ncols, sharey='all')
        for iPMT in range(self.nPMT):
            irow, icol = iPMT//ncols, iPMT%ncols
            label = 'S'+str(iPMT)
            ax[irow,icol].hist( results[:,iPMT])
            mean,std = numpy.mean(results[:,iPMT]),numpy.std(results[:,iPMT])
            title = 'mean {:.3f} input {:.3f}\nstd.dev. {:.3f}'.format(mean,Input[iPMT],std)
            ax[irow,icol].set_title(title,fontsize=7,y=1.0-0.11)
            ax[irow,icol].set_xlabel(label,labelpad=-125,loc='right')
            ax[irow,icol].axline( (Input[iPMT],0), (Input[iPMT],+1),color= 'red')
        wChi2 = ' $\chi^2=$ {:.2f}'.format(chi2PMT)
        Title = 'Results with method ' + method + wChi2 + '\n'+str(Nevts) + ' evts gen`d ' + timestamp
        fig.suptitle(Title)
        pdf = self.figDir + method + '.pdf'
        plt.savefig(pdf)
        print('twofold.plotResults Wrote',pdf)
        plt.show()
        return
    def readAndAnalyze(self,timestamp):
        '''
        read and analyze data from pickle file specified by timestamp
        '''
        Input, fitresults = self.readPickle(timestamp)
        self.input = Input
        CHI2ALL, CHI2PMT = self.analyzeResults(fitresults)
        methods = fitresults.keys()
        for method in methods:
            self.plotResults(fitresults,method,timestamp)
        return
    def analyzeMany(self,timestamps,method = 'L-BFGS-B'):
        '''
        read and analyze many toy mc given  pickle files specified by list of timestamps
        for a single minimization method
        '''
        bigDict = {}
        method = 'L-BFGS-B'
        for timestamp in timestamps:
            Input, fitresults = self.readPickle(timestamp)
            self.input = Input
            CHI2ALL, CHI2PMT = self.analyzeResults(fitresults)
            if method in CHI2ALL:
                bigDict[timestamp] = [Input,CHI2ALL, CHI2PMT]
        return bigDict
    def plotMany(self, method, bigDict):
        '''
        plot results from many toyMC
        '''
        print('twofold.plotMany method',method,'# points',len(bigDict))
        t_sort = sorted(bigDict) # sorted list of timestamps
        chi2pmts,means,prods,stds,spans = [],[],[],[],[]
        x = []
        for ix,timestamp in enumerate(t_sort):
            x.append(ix)
            Input, CHI2ALL, CHI2PMT = bigDict[timestamp]
            chi2All, chi2PMT = CHI2ALL[method], CHI2PMT[method]
            chi2pmts.append(chi2PMT)
            AA = Input[:self.Npmt]
            m,p,s,span = numpy.mean(AA), numpy.prod(AA), numpy.std(AA), numpy.max(AA)-numpy.min(AA)
            means.append(m)
            prods.append(p)
            stds.append(s)
            spans.append(span)
        
        bigY = [chi2pmts, means, prods, stds, spans]
        Ynames= ['$\chi^2(PMT)$','Mean','Prod','StdDev','Span']
        nrows,ncols = 5,1
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='all')
        for ia,Y in enumerate(bigY):
            name = Ynames[ia]
            ax[ia].plot(x,Y,'o-',label=name)
            ax[ia].set_ylabel(name)
            ax[ia].set_xticks(x)
            ax[ia].set_xticklabels(t_sort,rotation=90.)
            ax[ia].grid()
            rcorr, pval = pearsonr(chi2pmts, Y)
            text = 'r {:.3f} p {:.4f}'.format(rcorr,pval)
            print('twofold.plotMany corr between $\chi^2(PMT)$ and',name,text)
            if name=='$\chi^2(PMT)$' :
                ax[ia].set_ylim( (0.8*min(Y), 1.1*max(Y)) )
                ax[ia].set_yscale("log")
        plt.show()
        return
    def dPLOT(self,X1, X2):
        '''
        two-d plots
        '''
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2,ncols=2)
        im = ax1.imshow( numpy.array(X1), origin='lower' )
        fig.colorbar(im, ax=ax1)
        im = ax2.imshow( numpy.array(X2), origin='lower' )
        fig.colorbar(im, ax=ax2)
        im = ax3.imshow( numpy.array(X1)-numpy.array(X2), origin='lower')
        fig.colorbar(im, ax=ax3)
        ax1.set_ylabel('input')
        ax2.set_ylabel('true')
        ax3.set_ylabel('input-true')
        plt.show()
        return
    def getTimestamps(self):
        '''
        get list of timestamps that have directories with pickle files
        '''
        maindir = self.rootDir
        timestamps = []
        for d in os.listdir(maindir):
            md = maindir + d
            for root, dirs, files in os.walk(md):
                for file in files:
                    if file.endswith('pickle'): timestamps.append( d )
        print('twofold.getTimestamps Found',len(timestamps),'directories with pickle files')
        return sorted(timestamps)
    def MANY(self):
        '''
        main routine for analyzing many toy MC
        '''
        timestamps = self.getTimestamps()
        method = 'L-BFGS-B'
        bigDict = self.analyzeMany(timestamps,method=method)
        if len(bigDict)>0  : self.plotMany(method,bigDict)
        print('twofold.MANY log in',self.logDir)
        return
    def main(self):
        '''
        Fit real or fake data.

        These are minimization methods for scipy 1.8.1
        Unconstrained: CG, BFGS -- These methods frequently fail with message "Desired error not necessarily achieved due to precision loss." so omit them
        bound-constrained : 'Nelder-Mead','Powell','L-BFGS-B' --- best performance is achieved with the later
        '''

        debugPLOT = False ## activates use of dPLOT
        
        methods = []
        unconstrained_methods = [] # ['CG','BFGS']
        boundconstrained_methods = ['L-BFGS-B']
        self.bounded = {}
        for method in unconstrained_methods: self.bounded[method] = False
        for method in boundconstrained_methods: self.bounded[method] = True
        methods.extend( unconstrained_methods )
        methods.extend( boundconstrained_methods )

        methods = ["migrad"]
        
        toyMC = self.toyMC
        Nexpt = 1
        if toyMC:
            Nexpt = self.nToy
            Nguess = 5000.
            span = 0.8 # 0.4
            self.input = [1.+span/2. - numpy.random.random()*span for x in range(self.Npmt)]
            p = numpy.prod( self.input )
            print('twofold.main toyMC input effy',' %.2f'*len(self.input)%tuple(self.input),'span {:.2f} product {:.3f}'.format(span,p))
            self.FixedN = Nguess
            self.input.append( Nguess )
            
        fitresults = {m:[] for m in methods}
        for expt in range(Nexpt):
            words = 'twofold.main expt '+str(expt)
            
            if toyMC:   ##### Toy MC generation
                self.inputData = self.oneExpt()
                if debugPLOT :
                    Xtrue = self.oneTrueExpt()
                    self.dPLOT( self.inputData, Xtrue )
                chi2 = self.chisqr( self.input )
                words += ' {:} {:.3f}'.format('TruePar',chi2)
            else:       ##### data processing
                self.inputData = self.TwoFold['DATA']
                Nguess = self.maximum['DATA']
                self.FixedN = Nguess

            
            for method in methods:
                fitpar = self.find(method=method,Nguess=Nguess)
                fitresults[method].append( fitpar )
                chi2 = self.chisqr(fitpar)
                words += ' {:} {:.3f}'.format(method,chi2)
            print(words)
            #print('fitresults',fitresults)
        if toyMC :
            self.writePickle(self.input, fitresults)
            CHI2ALL, CHI2PMT = self.analyzeResults(fitresults)
            for method in methods:
                self.plotResults(fitresults,method,self.now)
        return
    def groupProb(self,group,dmc='DATA'):
        '''
        return calculated prob for group based on twofold coincidences
        '''
        rate = 0.
        notAlready = 1.

        for coinc in self.pmtCoinc[group]:
            if self.debug > 1 : print('twofold.groupRate group,coinc',group,coinc)
            if self.debug > 1 : print('twofold.groupRate dmc,coinc[0],self.probTwoFold[dmc][coinc[0]]',dmc,coinc[0],self.probTwoFold[dmc][coinc[0]])
            p = self.probTwoFold[dmc][coinc[0]][coinc[1]]
            if self.debug > 0: print('twofold.groupRate group,coinc,p,notAlready,rate',group,coinc,p,notAlready,rate)
            rate += p*notAlready
            notAlready *= (1. - p)
        return rate
    def getProbTF(self,scaleF,dmc=None):
        '''
        recalculate self.probTwoFold given scaleF and data or MC
        dmc = 'DATA' or 'MC'
        '''
        if dmc is None : sys.exit('twofold.getProbTF ERROR must specify dmc')
        sF = max(1.,scaleF)
        if scaleF<1. :
            print('twofold.getProbTF input scaleF',scaleF,'rescaled to',sF)
        self.probTwoFold[dmc] = self.TwoFold[dmc]/(sF*self.maximum[dmc])
        return
    def pairProb(self,gProb,gPair):
        '''
        return pairProb for the each of the group pairs in gPair
        '''
        pP = {}
        for pair in gPair:
            p1,p2 = gProb[pair[0]], gProb[pair[1]]
            pname = pair[0] + '_' + pair[1]
            pP[pname] = [p1*p2, p1*(1.-p2), p2*(1.-p1)]
        return pP
    def plotCResults(self,sF,pNames,results,Fractions=True,showMeas=False):
        '''
        plot pResults[dmc][j][pname][i] for j,x in enumerate(sF)
        dmc = 'DATA','MC'
        j = scale factor bin in sF
        pname = 'EVEN_ODD','RED_BLUE' or 'INNER_OUTER' in pNames
        i = 0,1,2 = A&B, A&~B, ~A&B
        if Fractions = True then calculate total = A&B + A&~B + ~A&B 
             and plot A&B/total, A&~B/total, ~A&B/total
        '''
        pts = ['-.',':','-']
        colpt = ['r','b','g']
        altcol= ['m','c','y']
        X = numpy.array(sF)
        x1,x2 = X[0],X[-1]
        for dmc in self.sources:
            for ip,pname in enumerate(pNames):
                ydict = {c:[] for c in self.cName}
                total = None

                for ic,c in enumerate(self.cName):
                    ydict[c] = [results[dmc][j][pname][ic] for j,jx in enumerate(sF)]
                    if total is None:
                        total = numpy.array(ydict[c])
                    else:
                        total+= numpy.array(ydict[c])
                for ic,c in enumerate(self.cName):
                    Y = numpy.array(ydict[c])
                    if Fractions : Y = Y/total
                    plt.plot(X,Y,pts[ip]+colpt[ic],label=c+' '+pname)
                    if showMeas:
                        y1 = y2 = self.measCoinc[dmc][pname][ic]
                        if Fractions: y1 = y2 = self.fracCoinc[dmc][pname][ic]
                        plt.plot((x1,x2),(y1,y2),altcol[ic]+pts[ip],label='meas '+c+' '+pname)
                        if self.debug > 2 : print('twofold.plotCResults dmc,pname,c,x1,y1,x2,y2',dmc,pname,c,x1,y1,x2,y2)
            plt.grid()
            plt.xlim(0,2.*x2)
            plt.legend(ncol=1,loc='right',fontsize='small')
            plt.title(dmc)
            plt.show()
        return
            
        
    def newMain(self):
        '''
        try to figure out distributions of EVEN/ODD, RED/BLUE, INNER/OUTER from twofold rates 20221011
        '''
        cmap = 'rainbow'
        im = [None,None]
        fig, (ax, axn) = plt.subplots(nrows=2,ncols=2)
        for i,x in enumerate(self.sources):
            print('twofold.newMain',x)
            print(self.TwoFold[x])
            im[i] = ax[i].matshow(self.TwoFold[x]/self.maximum[x],cmap=cmap)
            ax[i].set_title(x + ' normed')
            fig.colorbar(im[i],ax=ax[i])
        print('twofold.newMain DATA/MC')
        ''' from last answer to https://stackoverflow.com/questions/17870612/printing-a-two-dimensional-array-in-python '''
        numpy.set_printoptions(precision=3)
        print(numpy.matrix(self.TwoFold['DATA']/self.TwoFold['MC']))
        in0 = axn[0].matshow(self.TwoFold['DATA']/self.TwoFold['MC'],cmap=cmap)
        axn[0].set_title('Data/MC')
        fig.colorbar(in0, ax=axn[0])
        in1 = axn[1].matshow((self.TwoFold['DATA']/self.TwoFold['MC'])*(self.maximum['MC']/self.maximum['DATA']),cmap=cmap)
        axn[1].set_title('Data normed/MC normed')
        #axn[1].set_zlim(0,3.6)
        fig.colorbar(in1, ax=axn[1])
        plt.tight_layout()
        plt.show()

        fmt0 = '{:.3f} '
        fmt = fmt0+fmt0+fmt0

        sF = [float(x)/10. for x in range(10,20,2)]
        sF.extend( [float(x) for x in range(2,21)] )
        fmt2 = ''
        fmtN = ''
        for i in sF :
            fmt2 += fmt0
            fmtN += '{} '


        pResults = {}  # probabilities
        nResults = {}  # numbers of events
        
        for dmc in self.sources:
            pResults[dmc] = []
            nResults[dmc] = []
            for scaleF in sF:
                if self.debug > 0 : print('twofold.newMain dmc,scaleF',dmc,scaleF)
                self.getProbTF(scaleF,dmc=dmc)
                gP  = {}
                for group in self.pmtGroup:
                    gP[group] = pb = self.groupProb(group,dmc=dmc)
                    if self.debug > 0 : print('twofold.newMain',dmc,'group',group,'pb','{:.3f}'.format(pb))
                pProb = self.pairProb(gP,self.pmtPairs)
                pNames = list(pProb.keys())
                pResults[dmc].append( pProb )
                nProb = {}
                if self.debug > 1: print('twofold.newMain pProb',pProb)
                for pname in pProb:
                    if self.debug > 1 : print('twofold.newMain pname,pProb[pname]',pname,pProb[pname])
                    y = []
                    for x in pProb[pname]: y.append( int(scaleF*self.maximum[dmc]*x) )
                    nProb[pname] = y
                nResults[dmc].append( nProb )
                for pair in pProb:
                    if self.debug > 0 : print('twofold.newMain',pair,'pairProb',fmt.format(*pProb[pair]))
                    if self.debug > 0 : print('twofold.newMain',pair,'pairNum',nProb[pair])
        print(pNames)
        for dmc in self.sources:
            if self.debug > 0 : print('>>>',pResults[dmc])
            pname = pNames[0]
            i = 0
            j = 0
            if self.debug > 1 :
                print('>>>>',dmc,j,pResults[dmc][j])
                print('>>>>>',dmc,j,pname,pResults[dmc][j][pname])
                print('>>>>>',dmc,j,pname,i,pResults[dmc][j][pname][i])
            if self.debug > 0 : 
                for pname in pNames:
                    for i in range(3):
                        print('twofold.newMain',dmc,pname,self.cName[i],fmt2.format(*[pResults[dmc][j][pname][i] for j,x in enumerate(sF)]))
                    for i in range(3):
                        print('twofold.newMain',dmc,pname,self.cName[i],fmtN.format(*[nResults[dmc][j][pname][i] for j,x in enumerate(sF)]))
        plotCR = False
        if plotCR : 
            self.plotCResults(sF,pNames,nResults,Fractions=False,showMeas=True)
        return
if __name__ == '__main__' :
    debug = -1
    nToy  = 0 - 1
    timestamp = None
    if len(sys.argv)>1 :
        if 'help' in sys.argv[1].lower():
            print('USAGE: python twofold.py debug[{0}] nToy[{1}] timestamp[{2}]'.format(debug,nToy,timestamp))
            print('USAGE: if nToy < 0, then use newMain')
            print('USAGE: if timestamp is not None, then analyze data corresponding to timestamp')
            print('USAGE: if timestamp==`MANY`, then analyze data from many timestamps')
            sys.exit()
    
    if len(sys.argv)>1 : debug = int(sys.argv[1])
    if len(sys.argv)>2 : nToy  = int(sys.argv[2])
    if len(sys.argv)>3 : timestamp = sys.argv[3]
    
    P = twofold(debug=debug,nToy=nToy)
    if nToy < 0:
        P.newMain()
    else:
        if timestamp is None:
            P.main()
        elif timestamp=='MANY':
            P.MANY()
        else:
            P.readAndAnalyze(timestamp)

