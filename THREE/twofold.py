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
        if OK : print('twofold.__init__ Consistency check OK')
        self.maximum = maximum

        X = self.TwoFold['MC']
        m = maximum[src]
        if self.debug > 0 : print('twofold.__init__ self.TwoFold[`MC`]',X)
        self.prob = X/m
        if self.debug > 0 : print('twofold.__init__ prob',self.prob)
        print('twofold.__init__ maximum',maximum)

### store analysis results
### chi2All = sum_all [(best fit - input)/standard deviation]^2
### chi2PMT = sum_PMT [(best fit - input)/standard deviation]^2
### all includes number of input events

        self.analysisResults = {} # [method: [chi2All, chi2PMT]]
        

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
    def oneExpt(self):
        '''
        return 'Data' for one generated experiment given input efficiencies and number of events
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
    def oneTrueExpt(self):
        '''
        return 'Data' for one generated experiment given input efficiencies and number of events
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
        '''
        
        param = [1. for x in range(self.nPMT)]
        bounds = [(None,None) for x in range(self.nPMT)]
        bounds.append( (self.maximum['DATA'], 2.*Nguess) )
        param.append(Nguess)
        chi2 = self.chisqr(param)
        
        if self.debug > 0: 
            print('\ntwofold.find method',method)
            print('twofold.find chi2',chi2,'with input params',param)
        ## after 20220616T165320 don't give bounds for unbounded methods
        disp = self.debug > 1
        if self.bounded[method]: 
            res = minimize(self.chisqr, param, method=method,bounds=bounds)
        else:
            res = minimize(self.chisqr, param, method=method, options={'disp':disp})
            
        if self.debug > 1 : print('twofold.find res',res)
        pout = res.get('x')
        hess_inv = res.get('direc') #???
        success = res.get('success')
        chi2 = self.chisqr(pout)
        fitpar = 'twofold.find {0} chi2 {1:.1f} fitpar '.format(method,chi2)
        if self.debug > 1:
            print('twofold.find res.keys()',res.keys())
            print('twofold.find hess_inv',hess_inv)
        for i,p in enumerate(pout):
            punc = 0.
            if hess_inv is not None: punc = math.sqrt(max(0.,hess_inv[i,i]))
            fitpar += ' {0:.3f}({1:.3f})'.format(p,punc)

        if not success : fitpar += ' FIT FAILED. ' + res.get('message')
        if self.debug > 0 or not success: print(fitpar)
        return pout
    def analyzeResults(self,fitresults):
        '''
        compare fitresults with input for toyMC
        return Chi2All and Chi2PMT
        '''
        print('\ntwofold.analyzeResults')
        L = len(self.input)
        means,stddevs = {},{}
        CHI2ALL, CHI2PMT = {},{}
        for method in fitresults:
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
        if len(bigDict)>0 : self.plotMany(method,bigDict)
        return
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
            ax[ia].grid()
#            if name=='$\chi^2(PMT)$' : ax[ia].set_yscale("log")
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
        return timestamps
    def MANY(self):
        '''
        main routine for analyzing many toy MC
        '''
        timestamps = self.getTimestamps()
        method = 'L-BFGS-B'
        bigDict = self.analyzeMany(timestamps,method=method)
        if len(bigDict)>0  : self.plotMany(method,bigDict)
        return
    def main(self):

        ## these are minimization methods for scipy 1.8.1
        ## Unconstrained: CG, BFGS -- These methods frequently fail with message "Desired error not necessarily achieved due to precision loss." so omit them
        ## bound-constrained : 'Nelder-Mead','Powell','L-BFGS-B' --- best performance is achieved with the later
        methods = []
        unconstrained_methods = [] # ['CG','BFGS']
        boundconstrained_methods = ['L-BFGS-B']
        self.bounded = {}
        for method in unconstrained_methods: self.bounded[method] = False
        for method in boundconstrained_methods: self.bounded[method] = True
        methods.extend( unconstrained_methods )
        methods.extend( boundconstrained_methods )
        
        toyMC = self.toyMC
        Nexpt = 1
        if toyMC:
            Nexpt = self.nToy
            Nguess = 5000.
            span = 0.6 # 0.4
            self.input = [1.+span/2. - numpy.random.random()*span for x in range(self.Npmt)]
            p = numpy.prod( self.input )
            print('twofold.main toyMC input effy',' %.2f'*len(self.input)%tuple(self.input),'span {:.2f} product {:.3f}'.format(span,p))
            self.input.append( Nguess )
            
        fitresults = {m:[] for m in methods}
        for expt in range(Nexpt):
            words = 'twofold.main expt '+str(expt)
            
            if toyMC:   ##### Toy MC generation
                self.inputData = self.oneExpt()
                chi2 = self.chisqr( self.input )
                words += ' {:} {:.3f}'.format('TruePar',chi2)
            else:       ##### data processing
                self.inputData = self.TwoFold['DATA']
                Nguess = self.maximum['DATA']

            
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
if __name__ == '__main__' :
    debug = -1
    nToy  = 0
    timestamp = None
    if len(sys.argv)>1 :
        if 'help' in sys.argv[1].lower():
            print('USAGE: python twofold.py debug[{0}] nToy[{1}] timestamp[{2}]'.format(debug,nToy,timestamp))
            print('USAGE: if timestamp is not None, then analyze data corresponding to timestamp')
            print('USAGE: if timestamp==`MANY`, then analyze data from many timestamps')
            sys.exit()
    
    if len(sys.argv)>1 : debug = int(sys.argv[1])
    if len(sys.argv)>2 : nToy  = int(sys.argv[2])
    if len(sys.argv)>3 : timestamp = sys.argv[3]
    
    P = twofold(debug=debug,nToy=nToy)
    if timestamp is None:
        P.main()
    elif timestamp=='MANY':
        P.MANY()
    else:
        P.readAndAnalyze(timestamp)
    
