#!/usr/bin/env python
'''
utilities associated with particle trajectories

20140521
'''
import math
import sys
import random
import copy

class trajectory():
    def __init__(self):
        FWHM = 1. # cm at 2000 MeV measured in NSRL12C
        self.sigmaBeam = FWHM/2.35
        self.sigmaBeam = 2.0 # cm = sigmaH = sigmaV measured at 2000 MeV in NSRL13A (475MeV sigmaH = 2.5, sigmaV = 2.0 cm)
        self.bigSigma = 2.
        return

    def makeTrack(self,c,zc,option='Pencil'):
        '''
        make a track that passes through two hodoscopes with transverse
        dimensions c X c
        track is in Z direction between planes defined by zc
        '''
        dz = zc[1]-zc[0]
        if dz==0. :
            print 'makeTrack: ERROR zc',zc
            sys.exit(0)
        if option=='Pencil':
            x1 = y1 = x2 = y2 = 0.
        elif option=='Straight':
            x1 = x2 = random.uniform(-c/2.,c/2.)
            y1 = y2 = random.uniform(-c/2.,c/2.)
        elif option=='Spray':
            x1 = random.uniform(-c/2.,c/2.)
            y1 = random.uniform(-c/2.,c/2.)
            x2 = random.uniform(-c/2.,c/2.)
            y2 = random.uniform(-c/2.,c/2.)
        elif option=='Gaussian':
            x1 = random.gauss(0.,self.sigmaBeam)
            while abs(x1)>c/2.: x1 = random.gauss(0.,self.sigmaBeam)
            y1 = random.gauss(0.,sigmaBeam)
            while abs(y1)>c/2.: y1 = random.gauss(0.,self.sigmaBeam)
            x2 = x1
            y2 = y1
        elif option=='BigGauss':
            x1 = y1 = c
            while abs(x1)>c/2.: x1 = random.gauss(0., self.bigSigma)
            while abs(y1)>c/2.: y1 = random.gauss(0.,  self.bigSigma)
            x2 = x1
            y2 = y1
        elif option=='BigGaussLeft':
            x1 = y1 = c
            while abs(x1)>c/2.: x1 = random.gauss(c/2., self.bigSigma)
            while abs(y1)>c/2.: y1 = random.gauss(0.,  self.bigSigma)
            x2 = x1
            y2 = y1
        elif option=='BigGaussRight':
            x1 = y1 = c
            while abs(x1)>c/2.: x1 = random.gauss(-c/2.,self.bigSigma)
            while abs(y1)>c/2.: y1 = random.gauss(0.,self.bigSigma)
            x2 = x1
            y2 = y1
        elif option=='UpperRight':
            x1 = x2 = -c/2.
            y1 = y2 = c/2.
        elif option=='LowerRight':
            x1 = x2 = -c/2.
            y1 = y2 = -c/2.
        elif option=='UpperLeft':
            x1 = x2 = c/2.
            y1 = y2 = c/2.
        elif option=='LowerLeft':
            x1 = x2 = c/2.
            y1 = y2 = -c/2.
        else:
            print 'makeTrack: ERROR Invalid option',option
            sys.exit(1)
        track = [ [x1, y1, zc[0]], [x2, y2, zc[1]] ]
        return track

    def makePhoton(self,track, zIN, costheta, option='Cerenkov'):
        '''
        take an input track, generate an optical photon starting on the track and propagate it to the plane defined
        by z=zIN
        Assumes input track starts at z=0
        '''
        x1,y1,z1 = track[0][0],track[0][1],track[0][2]
        x2,y2,z2 = track[1][0],track[1][1],track[1][2]
        zg = random.uniform(z1,z2) # z of gamma emission point
        xg = x1 + (x2-x1)/(z2-z1)*(zg-z1)
        yg = y1 + (y2-y1)/(z2-z1)*(zg-z1)
        phi = random.uniform(0,2*math.pi)

        if option=='Cerenkov':
            sintheta = math.sqrt(1.-costheta*costheta)
            r = (zIN - zg)/costheta
            xk = xg + r*sintheta*math.cos(phi)
            yk = yg + r*sintheta*math.sin(phi)
            zk = zIN
        elif option=='Scint': # track either forward or backward based on random costheta
            ct = random.uniform(-1.,1.)
            st = math.sqrt(1.-ct*ct)
            zf = 0.
            if ct>0.: zf = zIN
            r = abs(zf-zg)/ct
            xk = xg + r*st*math.cos(phi)
            yk = yg + r*st*math.sin(phi)
            zk = zf
        else:
            print 'makePhoton: ERROR option',option
            sys.exit(3)
        return [ [xg,yg,zg], [xk,yk,zk] ]
    def getImpact(self,photonTrack, Cylinder):
        '''
        put the endpoint of the input photoTrack on the cylinder (radius, ztop, zbottom)
        '''
        X1,X2 = copy.copy( photonTrack )
        radcyl,ztop,zbottom = Cylinder
        x2,y2,z2 = X2
        r = math.sqrt(x2*x2 + y2*y2)
        if r<=radcyl: return photonTrack # hits bottom or top of cylinder
        x1,y1,z1 = X1
        # calculate distance to impact in x,y plane
        dxy = math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))
        ux = (x2-x1)/dxy
        uy = (y2-y1)/dxy
        b = x1*ux + y1*uy
        s = b*b - (x1*x1 + y1*y1 - radcyl*radcyl)
        if s<=0.:
            words = 'trajectory.getImpact: cannot calculate impact on cylinder wall'
            sys.exit(words)
        else:
            s = math.sqrt(s)
            dp = -b + s
            dm = -b - s
            if dp>=0. and dm<0.:
                d = dp
            elif dp<0. and dm>=0.:
                d = dm
            else:
                words = 'trajectory.getImpact: ERROR two solutions. dpos='+str(dp)+' dneg='+str(dm)
                sys.exit(words)
            xhit = x1 + ux*d
            yhit = y1 + uy*d
            zhit = z1 + (z2-z1)*d/dxy
        X2 = [xhit,yhit,zhit]
        # numskull check. initial and final directions same?
        col,dotprod = self.collinear( photonTrack,  [X1,X2]  )
        if not col:
            print 'trajectory.getImpact: ERROR inital',self.makeVec(photonTrack),'final',self.makeVec( [X1,X2] ),'photon track not collinear. dotprod',dotprod
        return [ X1, X2 ]
    def collinear(self, V1, V2):
        '''
        return 
        true if input vectors (defined by endpoints) are collinear
        and dotproduct
        '''
        col = self.trackAngle( V1, V2 )
        return col>(1.-1e-10), col
    def mag(self,v):
        m = 0.
        for i in range(len(v)):
            m += v[i]*v[i]
        return math.sqrt(m)
    def dot(self,v1,v2):
        m = 0.
        if len(v1)!=len(v2):
            print 'dot: ERROR len(v1),len(v2)',len(v1),len(v2),'v1',v1,'v2',v2
            sys.exit(5)
        for i in range(len(v1)):
            m += v1[i]*v2[i]
        m = m/self.mag(v1)/self.mag(v2)
        return m
    def trackAngle(self,track1, track2): # return cosine of angle between 2 tracks
        v1 = self.makeVec(track1)
        v2 = self.makeVec(track2)
        return self.dot(v1,v2)
    def makeVec(self,t):
        return [ t[1][0] - t[0][0],  t[1][1] - t[0][1],  t[1][2] - t[0][2] ]
    
