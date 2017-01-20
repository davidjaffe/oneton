#!/usr/bin/env python
'''
analyze waveforms presented as numpy arrays
20160114
'''
import numpy
import graphUtils
from ROOT import TFile
import sys
from collections import Counter
import math

class wfanal():
    def __init__(self):
        print 'wfanal: initialized'
        return
    def pulseAnal(self,v,detel,pedsub=True,debug=0,nsd = 5.):
        '''
        main routine
        return (flat) pedestal estimate, std.dev. of pedestal estimate,
         # of pulses, # of sub-pulses in each pulse,
        estimated area of each pulse and estimated time of each pulse
        '''

        # pulse definition parameters        
        top = 5
        if pedsub: top = 100
        
        ipre,ipost = 4,4
        minWidth = 3
        hFrac    = 0.3

        if debug>0:
            print '\nwfanal.pulseAnal counter',detel

        # get pedestal and std dev
        ped,pedsd = self.getPed(v,top=top)
        if debug>0:
            print 'wfanal.pulseAnal top,ped,pedsd {0:} {1:.2f} {2:.3f}'.format(top,ped,pedsd)

        # find pulses
        iPulse,subPperP = self.findPulse(v,ped,pedsd,nsd=nsd)
        if debug>0:
            #print 'wfanal.pulseAnal iPulse',iPulse,'subPperP',subPperP
            self.dumpPulse(v,iPulse,subPperP,ped,ipre,ipost,'after findPulse')

        # join nearby pulses, record number of sub-pulses in each pulse
        jPulse,jsubPperP = self.joinPulses(iPulse,subPperP,ipre=ipre,ipost=ipost)
        if debug>0:
            #print 'wfanal.pulseAnal jPulse',jPulse,'jsubPperP',jsubPperP
            self.dumpPulse(v,jPulse,jsubPperP,ped,ipre,ipost,'after joinPulses')

        # get rid of pulses that are too narrow
        iPulse,subPperP = self.winnowPulses(jPulse,jsubPperP,minWidth=minWidth)
        if debug>0:
            self.dumpPulse(v,iPulse,subPperP,ped,ipre,ipost,'after winnowPulses minWidth='+str(minWidth))

        # calculate area and time for each pulse
        pArea,pTime = self.calcATB(v,iPulse,ipre=ipre,ipost=ipost, hFrac=hFrac)
        if debug>0:
            if len(iPulse)>0:
                print 'wfanal.pulseAnal: area,times for hFrac',hFrac
                for i,a in enumerate(pArea):
                    t = pTime[i]
                    print 'Pulse# {0:} a {1:.2f} t {2:.2f}'.format(i,a,t),
                print ''
        return ped,pedsd,iPulse,subPperP,pArea,pTime
    def calcATB(self,v,iPulse,ipre=3,ipost=3,hFrac=0.3):
        '''
        return Area and Time for input pulses defined by iPulse, ipre,ipost in
        WFD array v

        First compute baseline using linear extrapolation between points
        Then create baseline subtracted pulse
        Compute the area of both the pulse and the 'total' pulse (with pulse range extended by pre,ipost)
        Determine the MINIMUM PH of pulse (pulse is negative going)
        Determine time of pulse as crossing point at fraction hFrac of minimum PH
        '''
        pArea, pTime, pBaseline = [],[],[]
        ptArea, pminPH = [],[] # total area, min pulse height
        if len(iPulse)>0:
            for pair in iPulse:
                i1,i2 = pair
                i0 = max(0,i1-ipre)
                i3 = min(i2+ipost,len(v))
                X,Y = [],[]
                for i in range(i0,i1,1):
                    X.append(float(i))
                    Y.append(v[i])
                for i in range(i2+1,i3,1):
                    X.append(float(i))
                    Y.append(v[i])
                X,Y = numpy.array(X),numpy.array(Y)
                #slope,yInt = self.fitLine(X,Y)
                slope,yInt = self.fitLine(X,Y)
                pBaseline.append([slope,yInt])
                #print 'wfANAL.calcATB baseline slope,yInt',slope,yInt,'points(x,y,f(x)):',['{0} {1:.2f} {2:.2f},'.format(a,b,slope*a+yInt) for a,b in zip(X,Y)]
                #print 'wfANAL.calcATB baseline Qslope,QyInt',Qslope,QyInt,'points(x,y,f(x)):',['{0} {1:.2f} {2:.2f},'.format(a,b,Qslope*a+QyInt) for a,b in zip(X,Y)]
                
                X,Y = [],[]
                area,tArea,minPH,imin = 0., 0., 1.e20, None
                for i in range(i0,i3,1):
                    x = float(i)
                    dx = 1.
                    y = slope*x + yInt
                    h = v[i]-y
                    X.append(x)
                    Y.append(h)
                    if h<minPH:
                        minPH = h
                        imin  = i
                    tArea += h*dx
                    if i1<=i and i<=i2: area += h*dx
                pArea.append(area)
                ptArea.append(tArea)
                pminPH.append(minPH)
                i = imin
                fthres = hFrac*minPH
                while i>i0 and Y[i-i0]<fthres:
                    #print 'i',i,'Y[i-i0]',Y[i-i0],'fthres',fthres
                    i = i - 1
                ##print 'wfanal.calcATB i,i0',i,i0
                U = numpy.array(X[i-i0:i+2-i0])
                V = numpy.array(Y[i-i0:i+2-i0])
                s,yI = self.fitLine(U,V)
                t = (fthres-yI)/s
                #print 'wfanal.calcATB t',t,'s,yI',s,yI,'fthres',fthres,'i0,i1,i2,i3',i0,i1,i2,i3,'i,imin',i,imin,'U,V',U,V
                pTime.append(t)
        return pArea,pTime
    def fitLineOld(self,X,Y):
        '''
        least squares solution for line with points X,Y
        '''
        N  = float(len(X))
        ##print 'wfanal.fitLine:X',X,'Y',Y
        m,b = 0.,0.
        if N==0.:
            print 'wfanal.fitLine: ERROR zero length input array X=',X,'Y=',Y
        else:
            sy = sum(Y)
            syy= sum(Y*Y)
            sx = sum(X)
            sxy= sum(X*Y)
            if sy==0.:
                print 'wfanal.fitLine: sy=',sy
                sxx = sum(X*X)
                b = (sy/sx - sxy/sxx) / (N/sx - sx/sxx)
                m = (sy - b*N)/sx
            else:
                m = (sy/N - syy/sy)/(sx/N -sxy/sy)
                b = (sy-m*sx)/N
        return m,b
    def fitLine(self,X,Y):
        '''
        least squares solution for line with points X,Y
        calculate chi-square like
        '''
        N  = float(len(X))
        ##print 'wfanal.fitLine:X',X,'Y',Y

        sy = sum(Y)
        syy= sum(Y*Y)
        sx = sum(X)
        sxy= sum(X*Y)
        sxx= sum(X*X)
        results = []
        for i in range(2):
            m,b,chi = 0.,0.,1.e20
            if i==0:
                den = (N*sxx - sx*sx)
                if sx!=0. and den!=0.:
                    b = (sy*sxx- sxy*sx) / den
                    m = (sy - b*N)/sx
                    Z = Y-(m*X+b)
                    chi = sum(Z*Z)
            else:
                den = (sx*sy-sxy*N)
                if N!=0. and den!=0.:
                    m = (sy*sy- syy*N) / den
                    b = (sy-m*sx)/N
                    Z = Y-(m*X+b)
                    chi = sum(Z*Z)
            results.append( [m,b,chi] )
        #print 'wfanal.fitLine: results',results
        if results[0][2]<results[1][2]:
            m,b = results[0][0:2]
        else:
            m,b = results[1][0:2]
        return m,b
    def getPed(self,v,top=100):
        '''
        estimate pedestal for input waveform v using mode
        mode computed from weighted average of top most common entries in v
        returns pedestal and standard deviation in pedestal
        top = 5 if v is integer array (without subtraction of online pedestals)
        top = 100 if v is float array
        '''
        G = Counter(v).most_common(top)
        g = numpy.array(G,'d')
        g0,g1 = g[:,0],g[:,1]
        sg1 = sum(g1)
        wa = sum(g0*g1)/sg1
        var= sum(g1*(g0-wa)*(g0-wa))/sg1*(float(top)-1.)/float(top)
        sd = math.sqrt(var)
        return wa,sd
    def findPulse(self,v,ped,sd,nsd=5.):
        '''
        return indices of start,end of each pulse in input v
        given pedestal ped and standard deviation estimate of pedestal distribution
        pulse starts(ends) when it is more(less) than nsd standard deviations from pedestal
        also initialize return list of subPperP = sub-Pulses per Pulse
        '''
        w = abs(v-ped)/sd/nsd
        wI = numpy.asarray(w,'i')   # render as integer array
        wNZ= wI.nonzero()           # find all non-zero entries
        K = wNZ[0]
        i2 = i1 = 0
        iPulses = []   # indices of start,end of each pulse
        for i,x in enumerate(K):
            if i>0:
                if x==K[i-1]+1:
                    i2 = i
                else:
                    if i2!=i1 : iPulses.append( [ K[i1],K[i2] ] )
                    i2 = i1 = i
        subPperP = [1 for x in range(len(iPulses))]
        return iPulses,subPperP
    def joinPulses(self,iPulse,subPperP,ipre=3,ipost=3):
        '''
        join pulses with a start(end) that is within
        ipre(ipost) bins of a previous(subsequent) pulse
        given indices of start,end of pulses in iPulse
        '''
        L = len(iPulse)
        #print 'wfanal.joinPulses at entry L',L,'iPulse',iPulse,'subPperP',subPperP
        if L<=1: return iPulse,subPperP
        jPulse = []
        jsubPperP   = []

        ipair = 0
        while ipair<L:
            pair = iPulse[ipair]
            nsub = subPperP[ipair]
            i1,i2 = pair
            I1,I2 = pair
            if ipair>0:
                i1pre,i2pre = iPulse[ipair-1]
                if i1-ipre<=i2pre: #join to previous pulse
                    I1 = min(I1,i1pre)
                    nsub += 1
            if ipair+1<L:
                i1post,i2post = iPulse[ipair+1]
                if i2+ipost>=i1post: #join to subsequent pulse
                    I2 = max(I2,i2post)
                    ipair += 1
                    nsub += 1
            if [I1,I2] not in jPulse:
                jPulse.append([I1,I2])
                jsubPperP.append(nsub)
            ipair += 1

        # iterate?
        LJ = len(jPulse)
        #print 'wfanal.joinPulses before iterate LJ',LJ,'jPulse',jPulse,'jsubPperP',jsubPperP
        if LJ>1 and LJ<L: jPulse,jsubPperP = self.joinPulses(jPulse,jsubPperP,ipre=ipre,ipost=ipost)
        return jPulse,jsubPperP
    def winnowPulses(self,iPulse,subPperP,minWidth=3):
        '''
        only keep pulses from list iPulse that have a width greater than minWidth
        '''
        jPulse, jsubPperP = [],[]
        for pair,nsub in zip(iPulse,subPperP):
            i1,i2 = pair
            if i2-i1+1>minWidth:
                jPulse.append(pair)
                jsubPperP.append(nsub)
        #print 'wfanal.winnowPulses iPulse,subPperP',iPulse,subPperP,'jPulse,jsubPperP',jPulse,jsubPperP
        return jPulse,jsubPperP
    def dumpPulse(self,v,iPulse,subPperP,ped,ipre,ipost,words=''):
        '''
        dump pulse info
        '''
        print 'wfanal.dumpPulse',words,len(iPulse),'pulses found. ipre,ipost',ipre,ipost
        for pair,nsub in zip(iPulse,subPperP):
            i1,i2 = pair
            print 'Pulse#',iPulse.index(pair),'NsubP',nsub,'range',i1,i2,'pulse',
            j1 = max(0,i1-ipre)
            j2 = min(i2+ipost,len(v))
            for j,x in enumerate(v):
                if j>=j1 and j<=j2:
                     if j==i1: print '|',
                     print '{0:.1f}'.format(x-ped),
                     if j==i2: print '|',
            print ''
        return
    
                    
