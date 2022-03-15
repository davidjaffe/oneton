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

class calibf():
    def __init__(self):

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
    def fitExpt(self,fNPE,model='DEFAULT',muMC=None,PLOT=True,allSteps=False,useLogLike=False):
        '''
        return bfNPE, cf given input NPE distribution fNPE where
        bfNPE = best fit NPE distribution
        cf  = calibration factor applied to MC to give best fit
        '''
        debug = 0
        
        meanmu = numpy.mean(fNPE)

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
        binwidth = self.binwidth
        bins = numpy.arange(0.,3.*max(fNPE)+binwidth,binwidth)  
        dataHist,dataEdges = numpy.histogram(fNPE,bins=bins)

        mu = numpy.mean(fNPE)
        if muMC is not None : mu = muMC
        NPE, mcNPE = self.genExpt(NMC,mu,model=model)
        
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
                    #sigsqr = dataHist + MChist/rmcd/rmcd
                    #chisqr = numpy.sum( (dataHist - MChist)**2/sigsqr )
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
            ### plot results
            X,Y = numpy.array(CFval),numpy.array(chival)

            ylabel = r'$\chi^2$'
            if useLogLike : ylabel = '2*loglike'
            
            fig, (ax0,ax1,ax2) = plt.subplots(nrows=3,ncols=1)
            ax0.plot(X,Y,'.')
            ax0.plot((bfCF),(chimin),'x',color='red')
            ax0.set_xlabel('Calibration factor')
            
            ax0.set_ylabel(ylabel)

            ax1.plot(X,Y,'.')
            ax1.plot((bfCF),(chimin),'x',color='red')
            ax1.set_ylim((chimin-20.,chimin+1000.))
            ax1.set_xlim((max(0.,bfCF-.2),bfCF+.2))
            ax1.set_xlabel('Calibration factor')
            ax1.set_ylabel(ylabel)

            ax2.plot(X,Y,'.')
            ax2.plot((bfCF),(chimin),'x',color='red')
            ax2.set_ylim((chimin-5.,chimin+1000./5.))
            ax2.set_xlim((max(0.,bfCF-.1),bfCF+.1))
            ax2.set_xlabel('Calibration factor')
            ax2.set_ylabel(ylabel)

            plt.tight_layout()
            plt.show()
        
        return bfNPE, bfCF
            
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
    def oneExpt(self,muData=0.63,muMC=0.5,Nevt=1000,PLOT=True,allSteps=False,useLogLike=False):
        '''
        generate, analyze and plot one experiment
        given input parameters
        '''

        trueCF = muData/muMC

                
        mu = muData
        NPE, fNPE = self.genExpt(Nevt,mu,model='DEFAULT')

        bfNPE, bfCF = self.fitExpt(fNPE, model='DEFAULT', muMC=muMC, PLOT=PLOT, allSteps=allSteps, useLogLike=useLogLike)


        if PLOT : 
            binwidth = 1.
            bins = numpy.arange(0.,max(max(NPE),max(fNPE))+binwidth,binwidth)

            fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1)
            ax0.set_title('NPE no electronics response')
            ax0.hist(NPE,bins=bins)
            self.writeStats(NPE,ax0)

            binwidth = self.binwidth
            bins = numpy.arange(0.,max(max(NPE),max(fNPE))+binwidth,binwidth)

            ax1.set_title('NPE true CF {:.3f} bestfit CF {:.3f}'.format(trueCF,bfCF))
            Y,binedges = numpy.histogram(fNPE,bins)
            Yerr = numpy.sqrt(Y)
            X = 0.5*(binedges[1:]+binedges[:-1])
            ax1.errorbar(X,Y,fmt='o',color='black',yerr=Yerr,label='NPE mock data')
            self.writeStats(fNPE,ax1,yfactor=0.55,words='Data')


            weight = float(len(fNPE))/float(len(bfNPE))
            weights = numpy.full(len(bfNPE),weight)
            ax1.hist(bfNPE,bins,weights=weights,color='blue',label='Best fit',alpha=0.5)
            self.writeStats(bfNPE,ax1,yfactor=0.4,words='Best fit')

            plt.legend(loc='best')        
            plt.tight_layout()        
            plt.show()
                
        return bfCF
    def main(self):
        testSPE = False
        if testSPE : 
            for Ninput in [10,200,5000,50000]:
                response = self.SPE(Ninput)
                response = numpy.array( response )
                print('calibf.main Ninput',Ninput,'len(response)',len(response),'sum(response<0.)',sum(response<0.),'sum(response)/float(len(response))',sum(response)/float(len(response)))


        showOne = False
        if showOne :
            Nevt = 1000
            for useLogLike in [False,True]:
                for allSteps in [False,True]:
                    for muData in [ 0.5, 1.5, 15.]:
                        print('calibf.main showOne Nevt,muData',Nevt,muData,'useLogLike',useLogLike)
                        bfCF = self.oneExpt(muData=muData,muMC=1.3*muData,Nevt=Nevt,PLOT=True,allSteps=allSteps)
            sys.exit('calibf.main Completed showOne')
        

        
        exptData = {} # {iConfig : [ (Nexpt, Nevt, muData, trueCF), [bfCF0, bfCF1, ...] ] }
        
        Nevt = 1000
        Nexpt = 50
        iConfig = 0
        trueCFlist = [0.7, 1.4]
        muDataList = [0.5, 1.5, 5.0, 10., 15.]
        for muData in muDataList: 
            for trueCF in trueCFlist:
                muMC = muData / trueCF
                iConfig += 1
                exptData[iConfig] =  [ (Nexpt, Nevt, muData, trueCF), [] ] 
                params = exptData[iConfig][0]
                print('calibf.main iConfig',iConfig,'Nexpt {:} Nevt {:} muData {:.2f} trueCF {:.1f}'.format(*params))
                for iExpt in range(Nexpt):
                    bfCF = self.oneExpt(muData=muData, muMC=muMC, Nevt=Nevt, PLOT=False)
                    exptData[iConfig][1].append( bfCF )
                CFs = exptData[iConfig][1]
                CFs = numpy.array( CFs )
                mean,std = numpy.mean(CFs), numpy.std(CFs)
                delta = mean - trueCF
                pull  = 99.
                if std>0. : pull = delta/std
                print('calibf.main iConfig',iConfig,'best fit CF mean,stddev,mean-true,pull {:.3f} {:.3f} {:.3f} {:.3f}'.format(mean,std,delta,pull))
        # plot bias on CF vs mean NPE
        print('\ncalibf.main Plot some stuff')
        X,Ycf,Ycfn = {},{},{}
        for trueCF in trueCFlist:
            X[trueCF] = []
            Ycf[trueCF] = []
            Ycfn[trueCF] = []
            for iConfig in exptData:
                Nexpt, Nevt, muData, tCF = params = exptData[iConfig][0]
                #print('calibf.main trueCF',trueCF,'iConfig',iConfig,'Nexpt {:} Nevt {:} muData {:.2f} tCF {:.1f}'.format(*params))
                
                if trueCF==tCF :
                    CFs    = exptData[iConfig][1]
                    mean,std = numpy.mean(CFs), numpy.std(CFs)
                    X[trueCF].append( muData )
                    Ycf[trueCF].append( mean-trueCF )
                    pull = 0.
                    if std>0.001 :
                        pull = (mean-trueCF)/std
                    else:
                        pull = 0
                        print('calibf.main ODD std=',std,'CFs',CFs)
                    Ycfn[trueCF].append( pull )

        #print('calibf.main X',X,'Ycf',Ycf,'Ycfn',Ycfn)
        fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1)
        for trueCF in trueCFlist:
            ax0.plot(X[trueCF],Ycf[trueCF],'o',label='trueCF {:.2f}'.format(trueCF))
            ax0.set_xlabel('mean NPE')
            ax0.set_ylabel('CF bias')

            ax1.plot(X[trueCF],Ycfn[trueCF],'o',label='trueCF {:.2f}'.format(trueCF))
            ax1.set_xlabel('mean NPE')
            ax1.set_ylabel('CF pull')
        ax0.legend(loc='best')
        ax0.grid()
        ax1.legend(loc='best')
        ax1.grid()
        plt.tight_layout()
        plt.show()

                
        return
if __name__ == '__main__' :
    P = calibf()
    P.main()
