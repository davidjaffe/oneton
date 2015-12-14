#!/usr/bin/env python
'''
calculate particle range in a material
units: MeV, cm, g
20140722
'''
import sys
import os
import numpy
import math
import graphUtils
from ROOT import TFile, TF1, gROOT

class scint_time_dist():
    def __init__(self):

        

        return
    def getAve(self,meanPE,decayTime=1.,Nexpts=1000):
        '''
        determine average decay time assuming meanPE are detected
        '''

        t,t2,N = 0., 0., 0.
        for i in range(Nexpts):
            n = numpy.random.poisson(meanPE)
            if n>0:
                N += 1.
                tmeas = min(numpy.random.exponential(decayTime,n))
                t += tmeas
                t2+= tmeas*tmeas
        ave = -1.
        sgm = -1.
        if N>1:
            ave = t/N
            sgm = (N/(N-1.)*(t2/N-ave*ave))
            sgm = math.sqrt(sgm/N)
        return ave,N,sgm
    def earliestTime(self,meanPE,f1=1.,t1=1.,t2=-1.):
        '''
        return earliest time from time distribution f1*g(t1) + (1.-f1)*g(t2)
        where g(tau) = exp(-t/tau)/tau
        if t2<=0, set f1=1.
        '''
        n = numpy.random.poisson(meanPE)
        if n>0:
            if t2<0 or f1==1. or numpy.random.random()<f1:
                times = numpy.random.exponential(t1,n)
            else:
                times = numpy.random.exponential(t2,n)
            t = min(times)
            return t
        
        return None


    def timeDist(self,x,par):
        t = float(x[0])
        tau1 = float(par[0])
        tau2 = float(par[1])
        f1   = float(par[2])
        norm = float(par[3])
        #print t,tau1,tau2,f1,'input'
        f1 = min(1.,max(0.,f1))
        f = f1/tau1*math.exp(-t/tau1)
        if tau2>0 and f1<1:
            f += (1.-f1)/tau2*math.exp(-t/tau2)
        #print t,tau1,tau2,f1,f
        f = norm*f # 100 bins?
        return f
    def plotTimeDist(self,t1=1.,t2=10.,f1=0.7,norm=1.,name="timeDist",tmax=-1.,nx=100):
        tma = 5.*max(t1,t2)
        if tmax>0: tma = tmax
        f = TF1(name,self.timeDist,0.,tma,4)
        #print 'plotTimeDist: t1,t2,f1',t1,t2,f1
        f.SetParameter(0,t1)
        f.SetParName(0,"tau1")
        f.SetParameter(1,t2)
        f.SetParName(1,"tau2")
        f.SetParameter(2,f1)
        f.SetParName(2,"f1")
        f.SetParameter(3,norm)
        f.SetParName(3,"Normalization")
        f.SetNpx(nx)
        #for i in range(4): print 'plotTimeDist:',name,'i,par',i,f.GetParName(i),f.GetParameter(i)
        return f
    def myfunction(self,x, par):
        xx = float(x[0])
        #print 'x',x,'xx',xx,'par',par
        f = abs(float(par[0])*math.sin(float(par[1])*xx))
#        f = 1.
        return f
    def myfunc(self):

        f1 = TF1("myfunc",self.myfunction,0,10,2)
        f1.SetParameters(2,1)
        f1.SetParNames("constant","coefficient")
        #f1.Draw()
        return f1

    
    def simple(self,Ntries=10):
        '''
        study bias in lifetime and mean PE estimator for single exponential
        run Ntries toy experiments for each combination of #of triggers(Nexpts) and meanPE
        20151030
        '''
        gU = graphUtils.graphUtils()

        makeTable = True
        Graphs = {}
        Hists  = {}
        icol = 0
        name = 'tauMeas_vs_muTrue'
        tmg1 = Graphs[name] = gU.makeTMultiGraph(name)
        name = 'tauMeas_vs_muMeas'
        tmg2 = Graphs[name] = gU.makeTMultiGraph(name)
        name = 'muMeas_vs_muTrue'
        tmg3 = Graphs[name] = gU.makeTMultiGraph(name)
        name = 'tauBias_vs_muTrue'
        tmg4 = Graphs[name] = gU.makeTMultiGraph(name)
        name = 'muBias_vs_muTrue'
        tmg5 = Graphs[name] = gU.makeTMultiGraph(name)
        decayTime = 1.
        #Ntries = 100
        if makeTable: print 'Ntrig & True mean & bias(mean) & bias(time) \\\\ '
        for Nexpts in [ 10000, 50000, 100000, 500000]:
            tauMeas,muMeas,muTrue = [],[],[]
            dtM, dmM, dmT = [],[],[]
            tauB, muB, dtauB, dmuB,goodPE = [],[],[],[],[]
            if makeTable: print ' \\hline '
            for meanPE in [0.001, 0.01, 0.02, 0.03, 0.05, 0.1, 0.3, 1.]:
                tBias,muBias,tBias2,muBias2 = [],[],[],[]
                for tries in range(Ntries):
                    ave,Nevt,sgm = std.getAve(meanPE,Nexpts=Nexpts,decayTime=decayTime)
                    if Nevt<Nexpts:
                        eff = float(Nevt)/float(Nexpts)
                        mu = -math.log(1.-eff)
                        dmu = 0.
                        if eff<1.: dmu = math.sqrt(eff/Nexpts/(1.-eff))
                        #print Nexpts, meanPE,'mu',mu,dmu,'tau',ave,sgm

                        tauMeas.append(ave)
                        dtM.append(sgm)
                        muMeas.append(mu)
                        dmM.append(dmu)
                        muTrue.append(meanPE)
                        dmT.append(0.)
                        tBias.append(ave-decayTime)
                        tBias2.append(math.pow(ave-decayTime,2))
                        muBias.append(mu-meanPE)
                        muBias2.append(math.pow(mu-meanPE,2))
                name = 'tauBias_meanPE_'+str(int(meanPE*1000))+'mille_'+str(Nexpts)+'trigs'
                Hists[name] = gU.makeTH1D(tBias,name.replace('_',' '),name,nx=Ntries)
                name = 'muBias_meanPE_'+str(int(meanPE*1000))+'mille_'+str(Nexpts)+'trigs'
                Hists[name] = gU.makeTH1D(muBias,name.replace('_',' '),name,nx=Ntries)
                NN = float(len(tBias))
                if NN>1.:
                    tB = sum(tBias)/NN
                    stB= sum(tBias2)/NN
                    stB= math.sqrt(NN/(NN-1.)*(stB-tB*tB))/NN
                    NN = float(len(muBias))
                    mB = sum(muBias)/NN
                    smB= sum(muBias2)/NN
                    smB= math.sqrt(NN/(NN-1.)*(smB-mB*mB))/NN
                    if makeTable:
                        print '{0:d} & {1:.4f} & {2:.4f} $\\pm$ {3:.4f} & {4:.4f} $\\pm$ {5:.4f} \\\\'.format(Nexpts,meanPE,mB,smB,tB,stB)
                    else:                
                        print 'Nexpts {0:d} meanPE {1:.4f} bias(meanPE) {2:.4f}+-{3:.4f} bias(tau) {4:.4f}+-{5:.4f}'.format(Nexpts,meanPE,mB,smB,tB,stB)
                    tauB.append(tB)
                    dtauB.append(stB)
                    muB.append(mB)
                    dmuB.append(smB)
                    goodPE.append(meanPE)


            icol += 1

            x,ex = goodPE,[0. for q in range(len(goodPE))]
            y,ey = tauB,dtauB
            name = 'tauBias_vs_muTrue_'+str(Nexpts)+'trigs'
            g = Graphs[name] = gU.makeTGraph(x,y,name.replace('_',' '),name,ex=ex,ey=ey)
            gU.color(g,icol,icol,setMarkerColor=True)
            tmg4.Add(g)

            y,ey = muB,dmuB
            name = 'muBias_vs_muTrue_'+str(Nexpts)+'trigs'
            g = Graphs[name] = gU.makeTGraph(x,y,name.replace('_',' '),name,ex=ex,ey=ey)
            gU.color(g,icol,icol,setMarkerColor=True)
            tmg5.Add(g)

            x,ex = muTrue,dmT
            y,ey = tauMeas,dtM
            name = 'tauMeas_vs_muTrue_'+str(Nexpts)+'trigs'
            title = name.replace('_',' ')
            g = Graphs[name] = gU.makeTGraph(x,y,title,name,ex=ex,ey=ey)
            gU.color(g,icol,icol,setMarkerColor=True)
            tmg1.Add(g)

            x,ex = muMeas,dmM
            y,ex = tauMeas,dtM
            name = 'tauMeas_vs_muMeas_'+str(Nexpts)+'trigs'
            title = name.replace('_',' ')
            g = Graphs[name] = gU.makeTGraph(x,y,title,name,ex=ex,ey=ey)
            gU.color(g,icol,icol,setMarkerColor=True)
            tmg2.Add(g)

            x,ex = muTrue,dmT
            y,ey = muMeas,dmM
            name = 'muMeas_vs_muTrue_'+str(Nexpts)+'trigs'
            title = name.replace('_',' ')
            g = Graphs[name] = gU.makeTGraph(x,y,title,name,ex=ex,ey=ey)
            gU.color(g,icol,icol,setMarkerColor=True)
            tmg3.Add(g)

        ofn = "std.root"
        outfile = TFile(ofn,"RECREATE")
        for g in Graphs: outfile.WriteTObject( Graphs[g] )
        for h in Hists:  outfile.WriteTObject( Hists[h] )
        outfile.Close()
        print "%Wrote",len(Graphs),"graphs and",len(Hists),"hists to",ofn

        for tmg in [tmg1,tmg2,tmg3,tmg4,tmg5]:
            gU.drawMultiGraph(tmg, SetLogx=True, abscissaIsTime = False, drawLines=False)
        gU.drawMultiGraph(tmg3, SetLogx=True, SetLogy=True, abscissaIsTime = False, drawLines=False)

        return
    def moreRealistic(self):
        '''
        more realistic comparison of true and measured distributions.
        Assume double exponential for scint time dist
        
        '''
        gU = graphUtils.graphUtils()
        Objects = []
        if 0:
            f1 = std.myfunc()
            Objects.append(f1)
            f2 = std.plotTimeDist()
            Objects.append(f2)

        t1,t2 = 1., 10.
        tmax = 10.*t2
        nbins = 100
        Ntrig = 100000
        icol = 0
        TMG = {}
        for f1 in [0.05, 0.25, 0.5, 0.75, 0.95]:
            truePE, bias, ebias = [],[],[]
            icol += 1
            for meanPE in [0.01, 0.03, 0.05, 0.10, 0.50, 1.0]:
                truePar = [t1,t2,f1]
                v = []
                for i in range(Ntrig):
                    t = std.earliestTime(meanPE,t1=t1,f1=f1,t2=t2)
                    if t is not None:
                        #print i,t
                        v.append(t)

                descrip = 't1_is_'+str(int(1000*t1))+'ps_t2_is_'+str(int(1000*t2))+'ps_f1_is_'+str(int(100*f1))+'percent_meanPE_is_'+str(int(1000*meanPE))+'mille_'+str(Ntrig)+'trigs'
                name = 'Meas_time_dist_'+descrip
                h = gU.makeTH1D(v,name.replace('_',' '),name,xmi=0.,xma=tmax,nx=nbins)
                Objects.append(h)


                name = 'True_time_dist_'+descrip
                f3 = std.plotTimeDist(name=name,t1=t1,t2=t2,f1=f1,norm=float(len(v)),tmax=tmax,nx=nbins)
                Objects.append(f3)

                name = 'Fitted_time_dist_'+descrip
                f4 = std.plotTimeDist(name=name,t1=t1,t2=t2,f1=f1,norm=float(len(v)),tmax=tmax,nx=nbins)
                gROOT.ProcessLine("gROOT->SetBatch()")  # no popup
                h.Fit(f4,"LQ")
                Objects.append(f4)
                #fitpar = []
                b,eb = [],[]
                for i in range(3):
                    #fitpar.append( f4.GetParameter(i) )
                    b.append(f4.GetParameter(i)-truePar[i])
                    eb.append(f4.GetParError(i))

                #print '%',f4.GetName(),fitpar
                truePE.append(meanPE)
                bias.append( b )
                ebias.append( eb )

                if 0:
                    KS = h.KolmogorovTest(f3.GetHistogram())
                    print 'KS',KS,descrip

            x,ex = truePE,[0. for g in range(len(truePE))]
            descrip = 't1_is_'+str(int(1000*t1))+'ps_t2_is_'+str(int(1000*t2))+'ps_f1_is_'+str(int(100*f1))+'percent_'+str(Ntrig)+'trigs'

            for i in range(len(truePar)):
                yname = f4.GetParName(i)
                y,ey = [],[]
                for b,eb in zip(bias,ebias):
                    y.append(b[i])
                    ey.append(eb[i])
                tmgname = yname + 'Bias_vs_muTrue'
                if tmgname not in TMG:
                    print 'tmgname',tmgname
                    TMG[tmgname] = gU.makeTMultiGraph(tmgname)
                name = tmgname+'_'+descrip
                g = gU.makeTGraph(x,y,name.replace('_',' '),name,ex=ex,ey=ey)
                gU.color(g,icol,icol,setMarkerColor=True)
                Objects.append(g)
                TMG[tmgname].Add(g)

        writeROOTFile = True
        if writeROOTFile:
            ofn = 'timedist.root'
            outfile = TFile(ofn,"RECREATE")
            for h in Objects:
                print "%",h.GetName()
                outfile.WriteTObject(h)
            for h in TMG:
                print "%Multigraph",TMG[h].GetName()
                outfile.WriteTObject(TMG[h])
            outfile.Close()
            print "%Wrote",len(Objects)+len(TMG),"objects to",ofn
        else:
            print '%NO ROOT FILE WRITTEN'

        # draw fitted distributions
        for h in Objects:
            if 'Meas' in h.GetName():
                gU.drawFit(h,SetLogy=True, figdir='Figures/')

        # mysteriously crashes when exitting after drawing multigraphs
        drawMG = False
        if drawMG: 
            for h in TMG:
                gU.drawMultiGraph(TMG[h], abscissaIsTime=False , figdir='Figures/')
                gU.drawMultiGraph(TMG[h], abscissaIsTime=False, SetLogx=True, figdir='Figures/')
        else:
            print "%NO MULTIGRAPHS DRAWN"

                
        return
if __name__ == '__main__' :
    std = scint_time_dist()
    #std.simple(Ntries=100) 
    std.moreRealistic()
