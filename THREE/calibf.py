#!/usr/bin/env python
'''
MC study of determination of calibration factors which match
observed NPE distribution to scaled simulated NPE distribution

20220315
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker # https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib/19972993#19972993

import datetime
import Logger

import pickle
import copy

class calibf():
    def __init__(self):


        now = datetime.datetime.now()
        self.now = now.strftime('%Y%m%dT%H%M%S')

        self.rootDir = 'JOBS_CALIBF/'
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
        print('calibf.__init__ Output directed to stdout and',lf)
        
        self.speModels = {'DEFAULT': (1., 1./2.35)} # (mean, sigma)

        self.binwidth  = 1.0 # for NPE distributions with electronics response
        
        print('calibf.__init__ Completed')
        return
    def SPE(self,Ninput,model='DEFAULT'):
        '''
        return electronics response for Ninput single photoelectrons based on input model

        model = 'DEFAULT' is normal distribution with no negative elements allowed
        
        '''
        debug = -1

        
        if model not in self.speModels : 
            sys.exit('calibf.SPE ERROR Unknown model '+ model)

        N = Ninput
        if model=='DEFAULT' :
            mean,sigma = self.speModels[model]
            Nneg = N
            outputSPE = []
            while Nneg > 0 :
                if debug > 0 : print('calibf.SPE Ninput,Nneg',Ninput,Nneg)
                response = numpy.random.normal(loc=mean,scale=sigma,size=Nneg)
                Nneg = sum(response<0.)
                if Nneg>0 : # negative elements
                    sr = sorted(response,reverse=True)
                    outputSPE.extend( sr[:-Nneg] )
                else:
                    outputSPE.extend( response )
        return outputSPE
    def genEvent(self,mu,model='DEFAULT'):
        '''
        return 
        NPE = integer number of photoelectrons assuming poisson distribution with mean mu, and
        fNPE= float of number of photoelectrons with electronics response included
        '''
        NPE = numpy.random.poisson(lam=mu)
        fNPE = 0.
        if NPE > 0 :
            response = self.SPE(NPE,model=model)
            fNPE = sum(response)
        return NPE, fNPE
    def genExpt(self,Nevt,mu,model='DEFAULT'):
        '''
        return NPE_gen, fNPE_gen for Nevt events with mean NPE mu and model as specified.
        where NPE_gen = array of integer number of photoelectrons, and
        fNPE_gen = array of float number of photoelectrons (with electronics response included)
        '''
        NPE_gen, fNPE_gen = [],[]
        for i in range(Nevt):
            NPE, fNPE = self.genEvent(mu,model=model)
            NPE_gen.append( NPE )
            fNPE_gen.append( fNPE )
        
        return numpy.array(NPE_gen), numpy.array(fNPE_gen)
    def NEWgenExpt(self,Nevt,mu,model='DEFAULT',Nbkg=None, beta=None):
        '''
        return NPE_gen, fNPE_gen for Nevt events with mean NPE mu and model as specified.
        where NPE_gen = array of integer number of photoelectrons, and
        fNPE_gen = array of float number of photoelectrons (with electronics response included)

        model = 'DEFAULT' = poisson distribution with mean mu
        model = 'BKGD1'   = (Nevt-Nbkg) events from poisson distribution with mean mu, Nbkg events from exponential exp(-x/beta)/beta
        '''

        N = Nevt
        model1 = model
        if model=='BKGD1' :
            model1 = 'DEFAULT'
            if Nbkg is None or beta is None :
                sys.exit('calibf.NEWgenExpt ERROR model='+model,' but required parameters Nbkg='+str(Nbkg)+' and beta='+str(beta))
            N = Nevt - Nbkg
        
        NPE_gen, fNPE_gen = [],[]
        for i in range(N):
            NPE, fNPE = self.genEvent(mu,model=model1)
            NPE_gen.append( NPE )
            fNPE_gen.append( fNPE )

        if model=='BKGD1' :
            fNPE_gen.extend( list( numpy.random.exponential(beta,Nbkg) ) )
        
        return numpy.array(NPE_gen), numpy.array(fNPE_gen)
    def fitExpt(self,fNPE,model='DEFAULT',muMC=None,PLOT=True,allSteps=False,useLogLike=False,debug=-1,comment=''):
        '''
        return bfNPE, cf given input NPE distribution fNPE where
        bfNPE = best fit NPE distribution
        cf  = calibration factor applied to MC to give best fit

        inputs
        fNPE = one experiment's of events. each event is NPE with electronics response applied
        model = model to be used to generate 'MC' distribution. 
        muMC = input mean PE for 'MC' distribution, if None, use mean(fNPE)
        PLOT = produce plots showing scan of chi2 vs CF
        allSteps = if True, take fixed steps between limits, if False, use iterative procedure
        useLogLike = if False, use chisquare method, if True, use logLikelihood 
        debug > 0 generates output
        comment is text to be used for output file (if PLOT=True)

        '''

        ## define scan method
        minCF = 0.2
        maxCF = 2.2
        if allSteps :
            nstep = 2000
            maxiter = 1
        else:
            nstep = 20
            maxiter = 3
        minstepsize = 0.001
        stepCF = max(minstepsize, (maxCF-minCF)/float(nstep))
        if allSteps : print('calibf.fitExpt allSteps',allSteps,'initial stepCF',stepCF)

        Ndata = len(fNPE)
        rmcd  = 10.
        NMC   = int(rmcd)*Ndata
        
        ## use bins >> maximum of input data to ensure that full range of calibration can be explored
        ## add large overflow bin instead of increasing overflow bins
        binwidth = self.binwidth
        bins = numpy.arange(0.,max(fNPE)+binwidth,binwidth)
        lastedge = bins[-1]
        bins = numpy.append(bins, lastedge + 200.*binwidth )
        dataHist,dataEdges = numpy.histogram(fNPE,bins=bins)
        if debug > 0 : print('calibf.fitExpt bins',bins)

        ## generate 'MC' experiment
        mu = numpy.mean(fNPE)
        if muMC is not None : mu = muMC
        NPE, mcNPE = self.genExpt(NMC,mu,model=model)
        
        
        meanmu = numpy.mean(fNPE)
        words = 'fitExpt_meanData_{:}_meanMC_{:}permille'.format(int(1000*meanmu),int(1000*mu)) + comment
        if comment!='' : print('calibf.fitExpt',words)


        chimin = 1.e20
        bfCF     = None
        bfNPE  = None
        CFval = []
        chival= []
        for iter in range(maxiter):
            for CF in numpy.arange(minCF,maxCF,stepCF):
                MChist,MCedges = numpy.histogram(CF*mcNPE,bins=bins)
                if useLogLike :
                    chisqr = 0.
                    nterm = 0
                    for d,m in zip(dataHist,MChist):
                        chisqr += m-d
                        if d>0. and m>0.:
                            nterm += 1
                            chisqr += d*numpy.log(d/m)
                    chisqr = 2.*chisqr
                else: 
                    chisqr = 0.
                    nterm = 0
                    for d,m in zip(dataHist,MChist):
                        sigsqr = d + m/rmcd/rmcd
                        if sigsqr>0 :
                            chisqr += (d-m)*(d-m)/sigsqr
                            nterm += 1
                if debug > 1 : print('calibf.fitExpt iter, CF,chisqr,nterm {} {:.2f} {:.0f} {}'.format(iter,CF,chisqr,nterm))
                CFval.append( CF )
                chival.append( chisqr )
                if chisqr<chimin or bfCF is None:
                    chimin = chisqr
                    bfCF = CF
            minCF = bfCF - stepCF
            maxCF = bfCF + stepCF
            stepCF = max(minstepsize, (maxCF-minCF)/float(nstep))
        if debug > 0 : print('calibf.fitExpt chimin,bfCF {:.0f} {:.3f}'.format(chimin,bfCF))

        bfNPE = bfCF * mcNPE

        if PLOT : 
            ### plot results. Show full scan and then two progressively finer zoom-ins around minimum
            X,Y = numpy.array(CFval),numpy.array(chival)
            npts = len(CFval)

            ylabel = r'$\chi^2$'
            if useLogLike : ylabel = '2*loglike'
            
            fig, (ax0,ax1,ax2) = plt.subplots(nrows=3,ncols=1)
            ax0.plot(X,Y,'.')
            ax0.plot((bfCF),(chimin),'x',color='red')
            method = r' $\chi^2$ method'
            if useLogLike : method = ' logLikelihood method'
            ax0.set_xlabel('Calibration factor Npoints {}'.format(npts) + method)
            
            ax0.set_ylabel(ylabel)
            ax0.grid()

            ax1.plot(X,Y,'.')
            ax1.plot((bfCF),(chimin),'x',color='red')
            ax1.set_ylim((chimin-20.,chimin+1000.))
            xlim = self.defineXlim(X,Y,chimin+1000.)
            ax1.set_xlim(xlim) 
            ax1.set_xlabel('Calibration factor')
            ax1.set_ylabel(ylabel)
            ax1.grid()

            ax2.plot(X,Y,'.')
            ax2.plot((bfCF),(chimin),'x',color='red')
            ax2.set_ylim((chimin-5.,chimin+200. ))
            xlim = self.defineXlim(X,Y,chimin+200. )
            ax2.set_xlim(xlim) 
            ax2.set_xlabel('Calibration factor')
            ax2.set_ylabel(ylabel)
            ax2.grid()

            plt.tight_layout()

            pdf = self.Figures + words + '.pdf'
            plt.savefig(pdf)
            print('calibf.fitResult Wrote',pdf)
            
            plt.show()
        
        return bfNPE, bfCF
    def defineXlim(self,Xin,Yin,ymax):
        '''
        return abscissa limits such that they contain the most extreme values of Yin<ymax. 
        input Xin,Yin = arrays; ymax = maximum y value

        '''
        X,Y = numpy.array(Xin), numpy.array(Yin)
        igood = []
        for i,y in enumerate(Y-ymax):
            if y<0: igood.append(i)

        if len(igood)>0: 
            i1 = max(0,igood[0]-1)
            i2 = min(len(X)-1,igood[-1]+1)
            xlo = X[i1]
            xhi = X[i2]
        else:
            xlo = X[0]
            xhi = X[-1]
        xlim = [xlo,xhi]
        return xlim
    def stats(self,a):
        '''
        return mean, standard deviation, median of input array a
        '''
        return numpy.mean(a), numpy.std(a), numpy.median(a)
    def writeStats(self,a,ax,xfactor=0.6, yfactor=0.8, words=''):
        '''
        write stats for array a on axis object ax 
        '''
        xlim,ylim = ax.get_xlim(),ax.get_ylim()
        mean,std,median = self.stats(a)
        x1 = xfactor*(xlim[1]-xlim[0])+xlim[0]
        y1 = yfactor*(ylim[1]-ylim[0])+ylim[0]
        ax.text(x1,y1,words+r' $\mu$ {:.2f} $\sigma$ {:.2f} $\tilde x$ {:.2f}'.format(mean,std,median))
        return
    def oneExpt(self,muData=0.63,muMC=0.5,Nevt=1000,PLOT=True,allSteps=False,useLogLike=False,debug=-1,comment='',model='DEFAULT'):
        '''
        generate, analyze and plot one experiment

        return bfCF, muGen, medianGen
        bfCF = best fit calibration factor
        muGen, medianGen = mean, median of MC generated NPE distribution without application of calibration factor 20220406

        Given input parameters
        muData = mean of poisson distribution for fake Data
        muMC   = mean of poisson dist for 'MC'
        Nevt = events to generate for this experiment
        PLOT = if True, show generated and fitted fake data distribution
        allSteps, useLogLike = Inputs to fitExpt
        debug > 0 generates output
        comment = text to append to file name if PLOT=True
        model = model used for generation of NPE distributions
        '''

        trueCF = muData/muMC

        if comment!='' : print('calibf.oneExpt muData',muData,'muMC',muMC,'trueCF',trueCF,comment)
                
        mu = muData
        NPE, fNPE = self.genExpt(Nevt,mu,model=model)

        bfNPE, bfCF = self.fitExpt(fNPE, model=model, muMC=muMC, PLOT=PLOT, allSteps=allSteps, useLogLike=useLogLike,debug=debug,comment=comment)
        words = 'oneExpt_muData_{:}_muMC_{:}_Nevt_{:}'.format(int(1000*muData),int(1000*muMC),Nevt) + comment
        
        muGen,medianGen = -1.,-1. # mean,median of MC generated NPE distribution
        if bfCF>0. :
            mcNPE = bfNPE/bfCF
            muGen = numpy.mean( mcNPE )
            medianGen = numpy.median( mcNPE )

        if PLOT :
            ## two panels:
            ## top panel is fake data NPE distribution w/o electronics effects
            ## bottom panel is fake data NPE dist with electronics effects applied
            ##      overlaid with best fit NPE distribution from 'MC'
            binwidth = 1.
            bins = numpy.arange(0.,max(max(NPE),max(fNPE))+binwidth,binwidth)

            fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1)
            ax0.set_title('NPE no electronics response')
            ax0.hist(NPE,bins=bins)
            self.writeStats(NPE,ax0)

            binwidth = self.binwidth
            bins = numpy.arange(0.,max(max(NPE),max(fNPE))+binwidth,binwidth)

            method = r' $\chi^2$ method'
            if useLogLike : method = ' loglikelihood method'
            ax1.set_title('NPE true CF {:.3f} bestfit CF {:.3f}'.format(trueCF,bfCF) + method)
            Y,binedges = numpy.histogram(fNPE,bins)
            Yerr = numpy.sqrt(Y)
            X = 0.5*(binedges[1:]+binedges[:-1])
            ax1.errorbar(X,Y,fmt='o',color='black',yerr=Yerr,label='NPE mock data')
            self.writeStats(fNPE,ax1,yfactor=0.6,words='Data')

            weight = float(len(fNPE))/float(len(bfNPE))
            weights = numpy.full(len(bfNPE),weight)
            ax1.hist(bfNPE,bins,weights=weights,color='blue',label='Best fit',alpha=0.5)
            self.writeStats(bfNPE,ax1,yfactor=0.3,words='Best fit')

            plt.legend(loc='best')        
            plt.tight_layout()
            pdf = self.Figures + words + '.pdf'
            plt.savefig(pdf)
            print('calibf.oneExpt Wrote',pdf)
            plt.show()
                
        return bfCF, muGen, medianGen
    def plotExptData(self,exptData):
        '''
        plots results given dict exptData
        two additional lists were added to exptData on 20220406. Backward compatibility should be enforced...
##OLD        exptData = {} # {iConfig : [ (Nexpt, Nevt, muData, trueCF), [bfCF0, bfCF1, ...] ] }
        exptData = {} # {iConfig : [ (Nexpt, Nevt, muData, trueCF), [bfCF0, bfCF1, ...], [muGen0, muGen1, ... ], [mediaGen0, medianGen1,...] }

        iConfig is key that distinguishes between different configuations specified by parameters
        (Nexpt, Nevt, muData, trueCF). 
        [bfCF0, bfCF1, ..., bfCFi, ...] are the best-fit CF from experiments 0, 1, ..., i, ... 
        [muGen0, muGen1, ..., muGeni, ..] is the mean of the MC generated experiments 0, 1, ..., i, ..
        [medianGen0, medianGen1, ...]     idem for median

        '''
        
        print('\ncalibf.plotExptData Plot some stuff')
        print('calibf.plotExptData sorted(exptData.keys())',sorted(exptData.keys()))

        # repopulate some lists
        tCF0 = None
        trueCFlist = []
        Expts,Events,muDataList = [],[],[]
        biasLimits = [1.e20,-1.e20]
        for iConfig in exptData:
            Nexpt, Nevt, muData, trueCF = exptData[iConfig][0]
            CFs = exptData[iConfig][1]
            bias = numpy.array(CFs)-trueCF
            biasLimits[0] = min(biasLimits[0],min(bias))
            biasLimits[1] = max(biasLimits[1],max(bias))
            if trueCF not in trueCFlist : trueCFlist.append( trueCF )
            if tCF0 is None : tCF0 = trueCF
            if trueCF==tCF0 :
                Expts.append( Nexpt )
                Events.append( Nevt )
                muDataList.append( muData )


        binwidth = 0.01/2
        bins = numpy.arange(biasLimits[0],biasLimits[1]+binwidth,binwidth)

        muDataList = sorted(muDataList)
        trueCFlist = sorted(trueCFlist)
        ncols = len(set(trueCFlist))
        nrows = len(set(muDataList))
        print('calibf.plotExptData ncols=# trueCF=',ncols,'nrows=# muData=',nrows)
        fontsize = 6
        if nrows*ncols>12 : fontsize = 5.5
                
        # plot distribution of measured CF - true CF for each set of experiments
        # also fill dicts used later to plot bias, sigma, etc.
        ### # figsize=(6.4,4.8) =default figure size in inches width, height matplot3.5.1
        FIG, AX = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,constrained_layout=True, figsize=(1.1*6.4,1.3*4.8)) 
        
        xlimits = [1.e20,-1.e20]
        X,Ycf,eYcf,Ycfn,Ystd = {},{},{},{},{}
        Ycf1,eYcf1 = {},{} ## cf1 = muData/muGen
        tCF0 = trueCFlist[0]
        for iCF,trueCF in enumerate(trueCFlist):
            X[trueCF] = []
            Ycf[trueCF] = []
            eYcf[trueCF] = []
            Ycfn[trueCF] = []
            Ystd[trueCF] = []
            Ycf1[trueCF],eYcf1[trueCF] = [],[]
            for iConfig in exptData:
                Nexpt, Nevt, muData, tCF = params = exptData[iConfig][0]
                
                if trueCF==tCF :
                    CFs = exptData[iConfig][1]
                    CFs = numpy.array( CFs )
                    bias = CFs - tCF
                    xlimits[0] = min(xlimits[0],min(bias))
                    xlimits[1] = max(xlimits[1],max(bias))
                    iax = muDataList.index(muData)
                    label = r'X {:} v {:} $\mu$ {:.2f} tCF {:.1f}'.format(*params)
                    AX[iax][iCF].hist(bias,bins=bins)
                    AX[iax][iCF].set_title(label,fontsize=fontsize)
                    AX[iax][iCF].grid()
                    if len(exptData[iConfig])>2 :
                        muGens, medianGens = exptData[iConfig][2:]
                        CF1s = muData/numpy.array( muGens )
                        bias1 = CF1s - tCF
                        m1,s1 = numpy.mean(CF1s), numpy.std(CF1s)
                        Ycf1[trueCF].append( m1-trueCF )
                        eYcf1[trueCF].append( s1/numpy.sqrt(Nexpt) )
                    
                    mean,std = numpy.mean(CFs), numpy.std(CFs)
                    X[trueCF].append( muData )
                    Ycf[trueCF].append( mean-trueCF )
                    eYcf[trueCF].append( std/numpy.sqrt(Nexpt) )
                    Ystd[trueCF].append( std )
                    pull = 0.
                    if std>0.001 :
                        pull = (mean-trueCF)/std
                    else:
                        pull = 0
                        print('calibf.plotExptData ODD std=',std,'CFs',CFs)
                        print('calibf.plotExptData ODD',label)
                    Ycfn[trueCF].append( pull )

        for iax in range(len(AX)):
            for iCF in range(len(AX[0])):
                AX[iax][iCF].set_xlim( xlimits )
                AX[iax][iCF].tick_params(axis='both',which='major', labelsize=8)
                AX[iax][iCF].tick_params(axis='both',which='minor', labelsize=8)
                if iax==len(AX)-1 : AX[iax][iCF].set_xlabel('bias',fontsize=8)


        pdf = self.Figures + 'toy_CF_bias_table.pdf'
        plt.savefig(pdf)
        print('calibf.plotExptData Wrote',pdf)
        plt.show()
        
                        
        ### plot CF bias and standard deviation vs muData for different values of trueCF
        fig, (ax0,ax2) = plt.subplots(nrows=2,ncols=1, sharex=True)
        ax = (ax0,ax2)

        for trueCF in trueCFlist:
            ax0.errorbar(X[trueCF],Ycf[trueCF],fmt='o-',yerr=eYcf[trueCF], label='trueCF {:.2f}'.format(trueCF))
            if len(Ycf1[trueCF])>0 : ax0.errorbar(X[trueCF],Ycf1[trueCF],fmt='o:',yerr=eYcf1[trueCF],label=r'trueCF {:.2f} $\mu$'.format(trueCF))
            ax0.set_ylabel('CF bias')

            ax2.plot(X[trueCF],Ystd[trueCF],'o-',label='trueCF {:.2f}'.format(trueCF))
            ax2.set_xlabel('mean NPE')
            ax2.set_ylabel(r'$\sigma(CF)$')
            
        for a in [ax0]:
            ylim = a.get_ylim()
            yma = 1.1*max(abs(ylim[0]),abs(ylim[1]))
            ylim = [-yma,yma]
            a.set_ylim( ylim )
            
        # specifying my favorite tick spacing for bias axis (ordinate)
        # taken from https://stackoverflow.com/questions/12608788/changing-the-tick-frequency-on-x-or-y-axis-in-matplotlib/19972993#19972993
        tick_spacing = 0.1
        ax0.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        ylim = ax2.get_ylim()
        ylim = [0, 1.1*ylim[1]]
        ax2.set_ylim( ylim )
            
        ax0.legend(loc='best',ncol=2)
        ax0.grid()
        ax2.grid()
        plt.tight_layout()

        pdf = self.Figures + 'toy_CF_bias_sigma_vs_NPE.pdf'
        plt.savefig(pdf)
        print('calibf.plotExptData Wrote',pdf)
        plt.show()

                
        return
    def groupExpts(self, timestamps):
        '''
        return exptData collated from a number of jobs specified by the input list timestamps
        
        exptData = {} # {iConfig : [ (Nexpt, Nevt, muData, trueCF), [bfCF0, bfCF1, ...] ] }
        '''
        exptData = None
        keymax = 0
        for iloop,timestamp in enumerate(timestamps):
            if exptData is not None : keymax = max( exptData.keys() ) 
            eData = self.readPickle(timestamp)
            if exptData is None:
                exptData = copy.copy(eData)
            else:
                for jConfig in eData:
                    iConfig = jConfig + keymax
                    if iConfig in exptData:
                        print('calibf.groupExpts ERROR iConfig',iConfig,'already exists for jConfig',jConfig,'and keymax',keymax,'for timestamp',timestamp)
                        sys.exit('calibf.groupExpts ERROR Overwrite of dict')
                    exptData[iConfig] = eData[jConfig]
        return exptData
    def readPickle(self,timestamp):
        '''
        return exptData given job specified by timestamp
        '''
        pickle_fn = self.rootDir + timestamp + '/' + timestamp + '.pickle'

        f = open(pickle_fn, 'rb')
        exptData = pickle.load( f)
        f.close()
        print('calibf.readPickle Read',pickle_fn)

        return exptData
    def readAndPlot(self,timestamps):
        '''
        read in pickle files specified by timestamps, collate info and plot contents
        '''

        exptData = self.groupExpts(timestamps)
        self.plotExptData( exptData )
        return
        
    def main(self):
        testSPE = False  ############# TEST THE GENERATION OF SPE DISTRIBUTIONS
        if testSPE : 
            for Ninput in [10,200,5000,50000]:
                response = self.SPE(Ninput)
                response = numpy.array( response )
                print('calibf.main Ninput',Ninput,'len(response)',len(response),'sum(response<0.)',sum(response<0.),'sum(response)/float(len(response))',sum(response)/float(len(response)))
                sys.exit('calibf.main Completed testSPE')

        testBkgd = False
        if testBkgd :
            Nevt = 10000
            Nbkg =   int(.02*Nevt)
            mu = 12.
            beta = 25.
            model = 'BKGD1'

            NPE, fNPE = self.NEWgenExpt(Nevt,mu,model=model,Nbkg=Nbkg,beta=beta)
            binwidth = 1.
            bins = numpy.arange(0.,max(NPE)+2.*binwidth,binwidth)
            ovflow = bins[-2]+0.01*binwidth
            fNPE_trimmed = numpy.clip(fNPE,0.,ovflow)
            fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1)
            ax0.hist(NPE,bins=bins,label='NPE w/o bkgd or electronics response')
            ax0.grid()
            ax0.legend(loc='best')
            
            ax1.hist(fNPE_trimmed,bins=bins,label='NPE w bkgd and electronics response')
            for bar in ax1.containers[0]:
                if bar.get_x() <= ovflow <= bar.get_x()+bar.get_width() : bar.set_color('red')
            ax1.grid()
            ax1.legend(loc='best')
            
            plt.tight_layout()
            plt.show()
            sys.exit('calibf.main Completed testBkgd ---------------------------------------')
            

        showOne = False   ############# SHOW SINGLE EXPERIMENT RESULT (chi2 scan and NPE distributions)
        if showOne :
            print('\ncalibf.main showOne Generate a single experiment, or a series of single experiments with different fit methods and/or input configurations')
            Nevt = 10000
            trueCF = 0.5
            #trueCF = 1.0
            njob = 0
            for useLogLike in [False] : 
                for allSteps in [False,True]:
                    for muData in [ 0.5, 5.0, 15.]:
                        comment = 'Njob_{:}'.format(njob)
                        njob += 1
                        print('calibf.main showOne Nevt,muData',Nevt,muData,'useLogLike',useLogLike)
                        bfCF, muGen, medianGen = self.oneExpt(muData=muData,muMC=muData/trueCF,Nevt=Nevt,PLOT=True,allSteps=allSteps,debug=1,useLogLike=useLogLike,comment=comment)
                        CFmean, CFmedian = -1., -1.
                        if muGen > 0. : CFmean = muData/muGen 
                        if medianGen > 0. : CFmedian = muData/medianGen
                        print('calibf.main ShowOne bfCF {:.3f} muGen {:.2f} medianGen {:.2f}'.format(bfCF, muGen, medianGen))
                        print('calibf.main showOne trueCF {:.3f} CFmean {:.3f} CFmedian {:.3f}'.format(trueCF, CFmean, CFmedian) )
            sys.exit('calibf.main Completed showOne')
        

        print('\ncalibf.main Generate and fit multiple experiments. Pickle results')
        exptData = {} # {iConfig : [ (Nexpt, Nevt, muData, trueCF), [bfCF0, bfCF1, ...], [muGen0, muGen1, ... ], [mediaGen0, medianGen1,...] }
        

        iConfig = 0
        trueCFlist = [0.7, 1.4]

        trueCFlist = [1.0]

        trueCFlist = [0.4, 1.7]
        muDataList = [0.5, 1.5, 5.0, 10., 15.]
        Expts      = [100, 100, 50,  20,  20]
        nEvent  = 1000
        Events     = [nEvent for x in range(len(Expts))] 

        statStudy = False
        if statStudy :
        # study of statistics on CF determination 20220922
            muDataList = [1.5, 1.5, 1.5]
            Expts      = [100, 100, 50]
            Events     = [1000, 10000, 40000]
            nEvent     = None # not used
            print('calibf.main Study of statistics on CF determination')
        for jExpt,muData in enumerate(muDataList):
            Nevt = Events[jExpt]
            Nexpt= Expts[jExpt]
            for trueCF in trueCFlist:
                muMC = muData / trueCF
                iConfig += 1
                exptData[iConfig] =  [ (Nexpt, Nevt, muData, trueCF), [],[],[] ] 
                params = exptData[iConfig][0]
                print('calibf.main iConfig',iConfig,'Nexpt {:} Nevt {:} muData {:.2f} trueCF {:.1f}'.format(*params))
                for iExpt in range(Nexpt):
                    bfCF, muGen, medianGen = self.oneExpt(muData=muData, muMC=muMC, Nevt=Nevt, PLOT=False)
                    exptData[iConfig][1].append( bfCF )
                    exptData[iConfig][2].append( muGen )
                    exptData[iConfig][3].append( medianGen )
                CFs = exptData[iConfig][1]
                CFs = numpy.array( CFs )
                mean,std = numpy.mean(CFs), numpy.std(CFs)
                delta = mean - trueCF
                pull  = 99.
                if std>0. : pull = delta/std
                print('calibf.main iConfig',iConfig,'best fit CF mean,stddev,mean-true,pull {:.3f} {:.3f} {:.3f} {:.3f}'.format(mean,std,delta,pull))

        pickle_fn = self.logDir + '/' + self.now + '.pickle'
        f = open(pickle_fn, 'wb')
        pickle.dump( exptData, f)
        f.close()
        print('calibf.main Wrote',pickle_fn)

        self.plotExptData(exptData)
        
                
        return
if __name__ == '__main__' :
    P = calibf()

#20220316T153733 20220317T134917 20220317T090549
# [0.7, 1.4]      [0.4, 1.7]      [1.0]  = trueCF
    
    timestamps = None
    if len(sys.argv)>1 :
        timestamps = [a for a in sys.argv[1:]]

    if timestamps is not None :
        P.readAndPlot(timestamps)
    else:
        P.main()
