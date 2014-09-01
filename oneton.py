#!/usr/bin/env python
'''
simulation of the 1ton prototype
units mm, MeV, nm
20140521
'''
import math
import sys
import random
import trajectory
import photons
import datetime
import Range
import cerenkov
import numpy
import copy
import ParticleProperties

class oneton():
    def __init__(self):

        self.debug = 0

        # initialize external classes. Note self.Photons initialized below
        self.traj = trajectory.trajectory()
        self.Range = Range.Range()
        self.cerenkov = cerenkov.cerenkov()
        self.ParticleProperties = ParticleProperties.ParticleProperties()

        # defaults for definition of detector
        self.detectorName = None
        self.pawDir = None
        self.Detector = None
        self.PMTradius = None
        self.dEdx = None
        self.ScintYield = None
        self.generationLog = None
                

        # for defining cosmic ray hodoscope external to vessel
        c = 1.
        self.HodoWidth = c
        self.Hodoscope = [-c/2., c/2, -c/2, c/2] # x1,y1 , x2,y2 = edges of hodoscope


        self.ior = 1.33  # index of refraction


        #self.wlRange = [300., 600. ] # range of wavelengths for cerenkov light relevant for bialkalai pmt (nm)
        self.wlRange = [250., 700. ] # wider range of wavelengths for cerenkov light relevant for bialkalai pmt (nm)

        self.hc = 2.*math.pi*197.326972 # eV*nm or MeV*fm
        self.Photons = photons.photons(wl1=self.wlRange[0], wl2=self.wlRange[1] )

        self.downward = [0., 0., -1.] # unit vector pointing down
        self.upward   = [0., 0.,  1.] # unit vector pointing up
        self.twopi = 2.*math.pi
        self.options = {'Cerenkov':0, 'Scint':1}

        self.me = self.ParticleProperties.getMass('electron')
        self.betaSpec = {}
        self.betaBins = 10000

        self.savedEvents = {}
        self.savedEventKeys = {}


        return
    def defineDetector(self,detName='oneton'):
        '''
        module to define detector parameters
        sets a bunch of internal parameters
        '''
        
        if detName.lower()=='oneton':
            self.detectorName = 'OneTon'
            self.pawDir = '/Users/djaffe/work/paw/ONETON/'
            # construct detector
            # coordinate system has z=0 at top, positive z pointing up
            ID = 995.
            IH = 1250.
            thick = 25.4
            self.Detector = [ID/2., 0., -IH] # [radius,ztop, zbottom]
            self.PMTradius = 25.4
            self.dEdx = 0.2 # MeV/mm
            self.ScintYield = 100 # photons/MeV
        if detName.lower()=='dayabay':
            # for WbLS-loaded dayabay
            self.detectorName = 'DayaBay'
            self.pawDir = '/Users/djaffe/work/paw/DAYABAY/'
            # construct detector
            # coordinate system has z=0 at top, positive z pointing up
            ID = 3000. 
            IH = 3000. 
            thick = 10.
            self.Detector = [ID/2., 0., -IH] # [radius,ztop, zbottom]
            self.PMTradius = 200./2 # 8" PMT
            self.dEdx = 0.2 # MeV/mm
            self.ScintYield = 1000 # photons/MeV
        if self.detectorName is None:
            words = 'oneton.defineDetector detName ' + str(detName)
            sys.exit(words)
            
        # file to contain relation of output ntuple 
        self.generationLog = self.pawDir + 'AAAGeneration.log'

        return
    def betaSpectrum(self,T,Tmax):
        '''
        return value of beta spectrum for input kinetic energy and
        max KE. Neglect matrix element.
        from Elements of Nuclear Physics, Meyerhof, eqn4-141
        '''
        if T<=0. or T>=Tmax: return 0.
        p = math.sqrt(T*T + 2.*T*self.me)
        r = p * (T+self.me) * (Tmax-T)*(Tmax-T)
        return r
    def closest(self,a,l):
        f = lambda a,l:min(l,key=lambda x:abs(x-a))
        return f(a,l)
    def genBeta(self,Tmax):
        '''
        randomly generate electron KE from a beta decay given the max KE (Tmax).
        Check to see if a normalized cumulative integral spectrum has already been generated for the input Tmax,
        if not, generate it.
        '''
        dT = Tmax/float(self.betaBins)
        if Tmax in self.betaSpec:
            cumInt = self.betaSpec[Tmax]
        else:
            cumInt = []
            c = 0.
            for iT in range(self.betaBins):
                T = dT/2. + float(iT)*dT
                c += self.betaSpectrum(T,Tmax)
                cumInt.append( c )
            for i,b in enumerate(cumInt):
                cumInt[i] = b/c
            self.betaSpec[Tmax] = cumInt
        ## generate a beta energy
        a = random.uniform(0.,1.)
        f = self.closest(a,cumInt)
        d = a-f
        i = cumInt.index(f)
        if d>0: j = min(i+1,self.betaBins-1)
        if d<0: j = max(0,i-1)
        Ebeta = float(i)*dT + d/abs(cumInt[i]-cumInt[j])*dT
        return Ebeta
    def oneSimpleEvent(self, nCerenkov=1, nScint=1):
        zc = [self.Detector[1], self.Detector[2]]
        beamTrack = self.traj.makeTrack(self.HodoWidth, zc, option='Pencil')
        nC = nCerenkov
        if nC<0:
            nC = int( self.getMeanCerenkovPhotons(abs(zc[0]-zc[1]),self.wlRange[0],self.wlRange[1]) )
        nS = nScint
        if nS<0:
            nS = int( self.getMeanScintPhotons(length=abs(zc[0]-zc[1])) )
        
        ###print 'nC',nC,'nS',nS
        beta = 1.0 # velocity of particle
        costheta = 1/self.ior/beta
        photons = []
        for i in range(nC+nS):
            option = 'Cerenkov'
            if i>=nC: option = 'Scint'
            iopt = self.options[option]
            photonTrack = self.traj.makePhoton(beamTrack, zc[1], costheta, option=option)
            photonTrack = self.traj.getImpact(photonTrack, self.Detector)
            photons.append( [photonTrack, option, iopt] ) # [ [x1,y1,z1], [x2,y2,z2] ], 'words', 0/1

        return photons
    def normVector(self,vin):
        mag = self.traj.mag ( vin )
        vout = []
        for u in vin: vout.append( u/mag)
        if self.debug>2: 
            print 'oneton.normVect vin',vin,'vout',vout
        return vout
    def makeVector(self, tStart, tDir):
        '''
        determine endpoints of vector from tStart in direction tDir
        that intersects the detector vessel
        '''
        radius, ztop, zbottom = self.Detector
        x1,y1,z1 = tStart
        unitV = tDir 
        if z1<zbottom or ztop<z1 or radius<math.sqrt(x1*x1+y1*y1):
            words = 'oneton.makeVector ERROR Start position outside detector. x,y,z= ' + str(x1) + ', ' + str(y1) + ', ' + str(z1)
            sys.exit(words)

        if unitV[2]==0. : # deal with the perfectly horizontal trajectory
            z2 = z1
            ux,uy = unitV[0],unitV[1]
            B = 2*x1*ux + 2*y1*uy
            A = 1.
            C = x1*x1 + y1*y1 - radius*radius
            radical = math.sqrt(B*B - 4.*A*C)
            lp = (-B + radical)/2./A
            lm = (-B - radical)/2./A
            goodp = lp>0.
            x2p,y2p,x2m,y2m = None,None,None,None
            if goodp:
                x2p = ux*lp + x1
                y2p = uy*lp + y1
                goodp = abs(x2p*x2p + y2p*y2p - radius*radius)<1.e-6
            goodm = lm>0.
            if goodm:
                x2m = ux*lm + x1
                y2m = uy*lm + y1
                goodm = abs(x2m*x2m + y2m*y2m - radius*radius)<1.e-6
            if goodm and goodp:
                print 'oneton.makeVector: ERROR Two valid solutions? X1',x1,y1,z1,\
                      'X2p',x2p,y2p,z2,'X2m',x2m,y2m,z2,'unitV',unitV[0],unitV[1],unitV[2]
                sys.exit(3)
            elif not goodm and not goodp:
                print 'oneton.makeVector: ERROR No valid solutions? X1',x1,y1,z1,\
                      'unitV',unitV[0],unitV[1],unitV[2],\
                      'lp',lp,'x2p',x2p,'y2p',y2p,
                if x2p is not None: print 'drp',x2p*x2p + y2p*y2p - radius*radius,
                print 'lm',lm,'x2m',x2m,'y2m',y2m,
                if x2m is not None: print 'drm',x2m*x2m + y2m*y2m - radius*radius
                print ''
                sys.exit(4)
            else:
                if goodm: x2,y2 = x2m,y2m
                if goodp: x2,y2 = x2p,y2p
            vector = [ [x1,y1,z1], [x2,y2,z2] ]
            if self.debug>2: print 'oneton.makeVector perfectly horizontal trajectory. vector',vector
        else:  # upward or downward going trajectory, extrapolate to plane containing an endplate
            if (z1-zbottom)*unitV[2]<0:
                z2 = zbottom
            else:
                z2 = ztop
            dist = abs( (z1-z2)/unitV[2] )
            x2 = unitV[0]*dist + x1
            y2 = unitV[1]*dist + y1
            vector = [ [x1,y1,z1], [x2,y2,z2] ]
            col,dotprod = self.traj.collinear( vector, [[0,0,0],tDir])
            if not col:
                print 'oneton.makeVector: ERROR vector',vector,'tDir',tDir,'NOT COLLINEAR prior to getImpact. dotprod',dotprod
                print 'oneton.makeVector: dist',dist,'z2-z1',z2-z1,'unitV[2]',unitV[2]
            if self.debug>2:
                print 'oneton.makeVector: tStart',tStart,'tDir',tDir,'initial vector',vector
            # place end of vector on detector wall or end            
            vector = self.traj.getImpact( vector, self.Detector )
            if self.debug>2: print 'oneton.makeVector: up-/downward trajectory. final vector',vector
        return vector
    def makeOpticalPhotonTrack(self, tStart,tDir, D, gDir):
        '''
        given track start point tStart and direction, create photon trajectory
        starting at distance D along track and going in direction gDir to
        impact vessel
        '''
        X1 = []
        for u,t in zip(tDir,tStart):
            X1.append( u*D + t )
        vector = self.makeVector( X1, gDir)
        col,dotprod = self.traj.collinear( vector, [ [0,0,0], gDir])
        if not col:
            print 'oneton.makeOpticalPhotonTrack: tStart',tStart,'tDir',tDir,'D',D,'gDir',gDir,'X1',X1,'vector',vector,'col',col,'dotprod',dotprod
            print 'oneton.makeOpticalPhotonTrack: ERROR photon track',vector,'not collinear with unit vector',gDir,'dotprod',dotprod
        return vector
    def createUnitVector(self,cost=None,phi=None):
        '''
        create a unit vector using input spherical variables costheta and phi.
        if input variables are not present, randomly generate
        '''
        ct = cost 
        if cost is None: ct = random.uniform(-1.,1.)
        f = phi
        if phi is None: f = random.uniform(0,self.twopi)
        st = math.sqrt(1. - ct*ct)
        U = [math.cos(f)*st,math.sin(f)*st,ct]
        return U
    def makeRUV(self, tDir, cost, phi):
        '''
        return unit vector gDir that is at angle theta (acos(cost)) wrt tDir and rotated
        about tDir by phi
        '''
        ## define axis of theta rotation as unit vector orthogonal to tDir
        t = numpy.array(tDir)
        u = numpy.array([1.,0.,0.])
        r = numpy.cross(t,u)
        if sum(r)==0.:
            u = numpy.array([0.,1.,0.])
            r = numpy.cross(t,u)
        ## now do rotation about r
        v = self.ERvector(t,r,math.acos(cost))
        ## now do rotation about t of angle phi
        g = self.ERvector(v,t,phi)
        ### resolve ambiguity(?)
        #print 'makeRUV:t',t,'g',g,'numpy.dot(t,g)',numpy.dot(t,g),'cost',cost
        if abs(numpy.dot(t,g)-cost)>1e-10:
            print 'vector reversed cost',cost,'phi',phi
            g = -g
        # render as list instead of array
        gDir = [ g[0], g[1], g[2] ]
        return gDir
    def ERvector(self,x1,AXIS,phi):
        '''
        return unit vector x2 = rotation of x1 about AXIS by angle phi
        using Euler-Rodriques vector formulation
        '''
        a = math.cos(phi/2.)
        axis = numpy.array( AXIS )
        w = math.sin(phi/2.)*axis/math.sqrt(numpy.dot(axis,axis))
        wx = numpy.cross(w,x1)
        x2 = x1 + 2.*a*wx + 2*numpy.cross(w,wx)
        return x2
    def niceString(self,v):
        '''
        return  nice string for input list
        '''
        s = ''
        l = len(v)
        for i,f in enumerate(v):
            if type(f) is str:
                s += f
            elif type(f) is float:
                s += '{0:4.2f}'.format(f)
            else:
                s += str(f)
            if i<l-1: s+=', '
        return s
    def oneEvent(self, particle, material, KE, tStart, tDir=[0.,0.,-1.], processes=[]):
        '''
        generate an event taking into account energy loss, initial start point and
        initial direction.
        return a list of photon tracks and a list of wavelengths
        '''
        ### has event with these input parameters already been generated?
        inputPar = [particle, material, KE, tStart, tDir]
        found = False
        if len(self.savedEventKeys)<1:
            eKey = 1
            self.savedEventKeys[eKey] = inputPar
        else:
            for eKey in self.savedEventKeys:
                if inputPar==self.savedEventKeys[eKey]:
                    found = True
                    break
            if not found:
                eKey += 1
                self.savedEventKeys[eKey] = inputPar
            else:
                # deepcopy in next line is necessary to prevent track (and other dict contents) from alternation
                tDir, track, totR, finalKE, samples = copy.deepcopy( self.savedEvents[eKey] )
                if self.debug>1: print 'oneton.oneEvent: eKey',eKey,'retrieved track',track

        if not found:
            ### make sure direction vector is unit vector
            tDir = self.normVector( tDir )

            ### determine material thickness from start point given starting direction
            track = self.makeVector( tStart, tDir)
            thickness = self.traj.mag( self.traj.makeVec( track ) )

            ### determine total range along initial trajectory (no scattering)
            totR, finalKE, samples = self.Range.propagate(particle, material, KE, thickness, maxdEperSample=0.1)
            originalTrack = copy.deepcopy( track )
            self.savedEvents[eKey] = tDir, originalTrack, totR, finalKE, samples
            if self.debug>1: print 'oneton.oneEvent: eKey',eKey,'saving track',track
            if self.debug>0: 
                print '\n totR',totR,'finalKE',finalKE,'thickness',thickness
                if self.debug>1: self.Range.printSampleTable( samples )

        ### now create track with start and endpoints
        X1, X2 = track
        for i,u in enumerate(tDir):
            X2[i] = totR*u + X1[i]
        track = [ X1, X2]
        if self.debug>2: print 'oneton.oneEvent track',track

        if self.debug>0: 
            print 'oneton.oneEvent totR',totR,'finalKE',finalKE,'nSamples',len(samples)
            print 'oneton.oneEvent tStart',tStart,'tDir',tDir,'track',track

        ### generate cerenkov and scintillation photons on track using samples
        ### Cerenkov: for each sample generate the energy and cos(theta) of each photon
        ### Scint: generate the energy and a placeholder of -2. for cos(theta) for each photon
        ### then generate a trajectory for each photon. Label trajectory by type of photon.
        totalCerYield, totalSciYield = 0.,0.
        OpticalPhotons, OPcost = [],[] # energy, cos(theta)
        PhotonVectors = [] # [ [Xstart, Xend], option, iopt]
        for i,sample in enumerate(samples):
            KE, pathlength = sample
            for option in processes:
                if self.debug>1: print 'oneton.oneEvent process',option
                if option == 'Cerenkov':
                    y,OPElimits,OP,cost = self.cerenkov.getCerenkovYield(KE,pathlength,particle,material,makeOP=True)
                    totalCerYield += y
                    if self.debug>1: print 'oneton.oneEvent process',option,'yield',y
                    for ct in cost:
                        dist = random.uniform(0.,pathlength) # start OP at random point on this sample of track
                        ranphi = random.uniform(0.,self.twopi)
                        gDir= self.makeRUV(tDir,ct,ranphi)
                        vector = self.makeOpticalPhotonTrack(X1,tDir,dist,gDir)
                        PhotonVectors.append( [vector, option, self.options[option] ] )
                    OpticalPhotons.extend( OP )
                    OPcost.extend( cost )

                if option == 'Scint':
                    KEs = [KE]
                    if i+1<len(samples): KEs.append( samples[i+1][0] )
                    yScint,OPscint = self.getScintYield(KEs,material,makeOP=True)
                    totalSciYield += yScint
                    if self.debug>1: print 'oneton.oneEvent process',option,'yield',yScint
                    fake = []
                    for j in range(len(OPscint)):
                        dist = random.uniform(0.,pathlength)
                        gDir = self.makeRUV(tDir,random.uniform(-1.,1),random.uniform(0.,self.twopi))
                        vector = self.makeOpticalPhotonTrack(X1,tDir,dist,gDir)
                        PhotonVectors.append( [vector, option, self.options[option] ] )
                        fake.append( -2. )
                    OpticalPhotons.extend( OPscint )
                    OPcost.extend( fake )

            # move up the start point for the next sample
            for i,u in enumerate(tDir): 
                X1[i] = u*pathlength + X1[i]

                
        if self.debug>0: 
            print '\nGenerate photons: total Cerenkov Yield',totalCerYield,\
                  'total Scint Yield',totalSciYield,'NGen',len(OpticalPhotons)
            if self.debug>1 and len(OpticalPhotons)>0:
                self.cerenkov.printCerenkovPhotons( OpticalPhotons, OPcost, binned=True,binsize=0.2)
        
        return OpticalPhotons, PhotonVectors
    def getScintYield(self, KEs,material,makeOP=False):
        '''
        return mean number of scintillation photons given initial,final kinetic energy,
        optionally return a list of energies corresponding to generated photons
        '''
        if self.debug>2: 
            print 'oneton.getScintYield KEs',KEs,'material',material,'makeOP',makeOP
        y = self.getMeanScintPhotons(KEs=KEs)
        OP = []
        if makeOP and y>0.:
            n = numpy.random.poisson(y)
            for i in range(n):
                wl = self.Photons.getWavelength(process='Scint',medium=material)
                OP.append( self.hc/wl )
        return y,OP
    def getMeanCerenkovPhotons(self,length,wl1,wl2,beta=1.):
        E1 = self.hc/wl1
        E2 = self.hc/wl2
        costheta = 1/self.ior/beta
        s2t = 1. - costheta*costheta
        nG = 369.9 * s2t * abs(E1-E2) * length
        return nG
    def getMeanScintPhotons(self,length=None, KEs=None):
        '''
        determine mean number of scintillation photons created
        given either the length of the track or the initial and final kinetic energy (or just
        the final kinetic energy).
        (Future?): add particle type dependence, add energy dependence
        '''
        if (length is None and KEs is None) or (length is not None and KEs is not None):
            words = 'oneton.getMeanScintPhotons ERROR Inconsistent inputs length ' + str(length) + ' KEs ' + str(KEs)
            sys.exit(words)
        if length is not None:
            nG = length * self.dEdx * self.ScintYield
        if KEs is not None:
            if isinstance(KEs,list):
                if len(KEs)==1:
                    E1,E2 = KEs[0],0.
                else:
                    E1,E2 = KEs
            else:
                E1,E2 = KEs,0.
            nG = abs(E1-E2) * self.ScintYield
        return nG
    def placePMTs(self,outFile=None):
        radcyl, ztop, zbottom = self.Detector
        xyzPMT = []
        if self.detectorName.lower()=='oneton':
            ## bottom
            for i in range(10):
                x = 2.*float(i)*self.PMTradius
                if x+self.PMTradius<=radcyl:
                    xyzPMT.append( [x, 0., zbottom, self.PMTradius, False] ) # x,y,z,r, side or top/bottom pmt
                    if abs(x)>self.PMTradius:
                        xyzPMT.append( [-x, 0., zbottom, self.PMTradius, False] ) # x,y,z,r, side or top/bottom pmt

            ## side
            for i in range(30):
                z = zbottom + (2.*float(i)+1)*self.PMTradius
                if z+self.PMTradius<=ztop:
                    xyzPMT.append( [radcyl, 0., z, self.PMTradius, True] )
                    xyzPMT.append( [-radcyl, 0., z, self.PMTradius, True] )
            ## top
            for i in range(10):
                x = radcyl - (2.*float(i)+1)*self.PMTradius            
                if x-self.PMTradius>0:
                    xyzPMT.append( [x, 0., ztop, self.PMTradius, False] )
                    xyzPMT.append( [-x, 0., ztop, self.PMTradius, False] )
        if self.detectorName.lower()=='dayabay':
            position = 'Side'
            nrow = 8
            ncol = 24
            sep = (ztop-zbottom)/float(nrow)
            dphi = self.twopi/float(ncol)
            h = zbottom+sep/2.
            for irow in range(nrow):
                for iphi in range(ncol):
                    phi = float(iphi)*dphi
                    x = math.cos(phi)*radcyl
                    y = math.sin(phi)*radcyl
                    #print 'irow',irow,'phi',phi,'x',x,'y',y
                    xyzPMT.append( [x,y,h, self.PMTradius, True] )
                h += sep
                
        # produce a list of indices sorted by z position
        ZofPMT = [] # z,index (will be sorted)
        for i,pmt in enumerate(xyzPMT):
            ZofPMT.append( [ pmt[2],i ] )
        ZofPMT.sort()
        
        # make a table
        ncol = 4
        line = ''
        for icol in range(ncol):
            line += '{0:4} {1:>8} {2:>8} {3:>9}      '.format('PMT#','x','y','z')
        print line
        line = ''
        for i,pmt in enumerate(xyzPMT):
            position = 'Vert'
            if pmt[4]: position = 'Side'
            line += '{0:4} {1:>8.3f} {2:>8.3f} {3:>9.3f} {4}'.format(i,pmt[0],pmt[1],pmt[2],position)
            if (i+1)%ncol==0:
                print line
                line = ''
        if len(line)>0: print line
        # write to ntuple?
        if outFile is not None:
            iopt = -3
            for ipmt,e in enumerate(xyzPMT):
                X1 = e[:3] # position
                X2 = [e[3], 0., 0.] # radius
                option = 'Vert'
                if e[4] : option = 'Side'
                photon = [ [X1,X2], option, iopt ]
                self.makeNtuple(outFile, photon, ipmt, -1., -1., -1., -1)
        # check for mistakes, do PMTs overlap?
        self.checkPMToverlap(xyzPMT)
        return xyzPMT,ZofPMT
    def checkPMToverlap(self,xyzPMT):
        '''
        check if PMTs overlap. Error and exit, if so.
        xyzPMT = [ [x,y,z,r,T/F=side/top-bottom]

        Note tiny adjustment made to square of sum of radii to take into account
        numerical precision when PMTs are placed so that they are touching.
        '''
        for i,pmt1 in enumerate(xyzPMT):
            x1,y1,z1,r1,side1 = pmt1
            for j,pmt2 in enumerate(xyzPMT):
                if j>i:
                    x2,y2,z2,r2,side2 = pmt2
                    if side1==side2:
                        if (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) < (r1+r2)*(r1+r2)-1.e-10:
                            words = 'oneton.checkPMToverlap: ERROR pmt# ' + str(i) + ' overlaps pmt# '+ str(j) +' \n'
                            words+= 'x1='+str(x1)+' x2='+str(x2)+' x1-x2=' + str(x1-x2) +  ' \n'
                            words+= 'y1='+str(y1)+' y2='+str(y2)+' y1-y2=' + str(y1-y2) +  ' \n'
                            words+= 'z1='+str(z1)+' z2='+str(z2)+' z1-z2=' + str(z1-z2) +  ' \n'
                            words+= 'r1='+str(r1)+' r2='+str(r2)+' r1+r2=' + str(r1+r2) +  ' \n'
                            words+= 'sum of squares of differences='+str((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))+ ' square of sums of radii='+str((r1+r2)*(r1+r2))+' difference='+str((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)-(r1+r2)*(r1+r2))
                            
                            print words
                            sys.exit(words)
        return
    def hitPMT(self,photonT,xyzPMT,ZofPMT):
        '''
        determine hit PMT for input photon track photonT = [ [startPoint], [endPoint] ]
        given PMT positions in xyzPMT
        ZofPMT is list of [z,pmt index] sorted by z
        Check for >1 PMT hit supplanted by checkPMToverlap
        '''
        hits = []
        nhit = 0

        X1,X2 = photonT
        x,y,z = X2
        dist = 1.e20
        iPMT = -1 #  no PMT hit
        # first check if hit is below (above) first (last) PMT in z
        zPMT,jPMT = ZofPMT[0]
        xPMT,yPMT,zPMT,rPMT,side = xyzPMT[jPMT]
        if z<zPMT-rPMT: return iPMT
        zPMT,jPMT = ZofPMT[-1]
        xPMT,yPMT,zPMT,rPMT,side = xyzPMT[jPMT]
        if z>zPMT+rPMT: return iPMT
        # search from lowest to highest z, break when past the z of the hit
        for pair in ZofPMT:
            zPMT,jpmt = pair
            D = (z-zPMT)*(z-zPMT)
            if D>dist:
                return iPMT
            dist = D
            xPMT, yPMT, zPMT, rPMT, side = xyzPMT[jpmt]

            dz = 0.00001
            if side : dz = rPMT
            if abs(z-zPMT)<dz:
                if ( (x-xPMT)*(x-xPMT) + (y-yPMT)*(y-yPMT) + D )<rPMT*rPMT:
                    return jpmt
        return iPMT
    def makeNtuple(self,fopen, photon, ipmt, wl, qe, att, iEvt):
        '''
        for each photon   start     end       0/1   
        ntuple variables: x1,y1,z1, x2,y2,z2, iopt, ipmt, wl, qe, attenuation factor
        iopt = 0 = Cerenkov, 1 = Scintillation, -2 = drawDet, -3 = PMT positions
        for iopt = 3: x1,y1,z1 = PMT center, x2 = PMT radius
        
        replaced string concatentation with join of list
        '''
        track, option, iopt = photon
        X1, X2 = track
        L = []
        L.extend(X1)
        L.extend(X2)
        for u in [iopt, ipmt, wl, qe, att, iEvt, "\n"]: L.append(u)
        s = " ".join(str(a) for a in L)
        fopen.write(s)
        return
    def readNtuple(self,inFile=None,writeAnaNtuple=True):
        '''
        read ntuple
        
        '''
        if inFile is None:
            sys.exit('oneton.readNtuple ERROR No input file')
        f = open(inFile,'r')
        fana = None
        if writeAnaNtuple:
            unique = '_{0}'.format(datetime.datetime.now().strftime("%Y%m%d%H%M_%f"))
            anaName = inFile.replace('.nt',unique + '.ana')
            fana = open(anaName,'w')
            print 'opened',anaName
        PMTxyz = {}
        evtNum = None
        Event = []
        for line in f:
            w = line.split()
            iopt = int(w[6])
            if iopt==-3: # PMT positions
                x,y,z,r =  float(w[0]),float(w[1]),float(w[2]),float(w[3])
                ipmt = int(w[7])
                if ipmt in PMTxyz:
                    words = 'oneton.readNtuple: ERROR found duplicate PMT number ' + str(ipmt)
                    sys.exit(words)
                PMTxyz[ipmt]  = [x,y,z, r] 
            elif iopt==-2: # ignore, for detector drawing
                continue
            else: # should be data. iopt=0=cerenkov photons, iopt=1=scint photons
                ievt = int(w[-1])
                if ievt<0:
                    words = 'oneton.readNtuple: ERROR negative event in line: ' + line
                    sys.exit(words)
                # new event?
                if ievt!=evtNum: 
                    if evtNum is not None:
                        # analyze previous event
                        self.anaNtupleEvent(PMTxyz, Event, anaFile=fana)
                    Event = [ievt]
                    evtNum = ievt
                # add one photon to event
                xg1,yg1,zg1,xg2,yg2,zg2 = float(w[0]),float(w[1]),float(w[2]),float(w[3]),float(w[4]),float(w[5])
                ipmt = int(w[7])
                wl,qe,att = float(w[8]),float(w[9]),float(w[10])
                if qe<0. or att<0.:
                    words = 'oneton.readNtuple: ERROR negative QE or attenuation in line: ' + line
                    sys.exit(words)
                Event.append( [ [xg1,yg1,zg1], [xg2,yg2,zg2], iopt, ipmt, wl, qe, att] )
        f.close()
        # analyze final event
        if len(Event)>0:
            self.anaNtupleEvent(PMTxyz, Event, anaFile=fana)
            if fana is not None: fana.close()
        return
    def anaNtupleEvent(self, PMTxyz, Event, anaFile=None):
        '''
        analyze one event in ntuple
        given pmt positions in dict PMTxyz[ipmt] = [x,y,z,radius]
        and photons in list, one element = [ X1,X2, gtype, PMT#, wavelength, QE, attenuation ]
        gtype = 0 = cerenkov
        gtype = 1 = scint
        '''
        COGnum = {0:[0.,0.,0.], 1:[0.,0.,0.]}
        COGden = {0:0., 1:0.}
        nGamma = {0:0., 1:0.}

        report = False

        # am I an idiot?
        checkMe = False
        if checkMe:
            gtype = 0
            wt = 1.
            for ipmt in PMTxyz:
                PMTpos = PMTxyz[ipmt][:3]
                COGden[gtype] += wt
                for i in range(len(COGnum[gtype])): COGnum[gtype][i] += PMTpos[i]*wt
            wt = COGden[gtype]
            for i,e in enumerate(COGnum[gtype]):
                COGnum[gtype][i] = e/wt
            print 'oneton.anaNtupleEvent: checkMe gtype',gtype,'COG',COGnum[gtype],'wt',wt
            COGnum[gtype] = [0.,0.,0.]
            COGden[gtype] = 0.
            
        evtNum = None
        Xfirst= None
        for photon in Event:
            if type(photon) is int:
                evtNum = photon
            else:
                X1, X2, gtype, ipmt, wl, qe, att = photon
                if Xfirst is None: Xfirst = X1
                Xlast = X1

                if ipmt>-1:
                    nGamma[gtype] += 1
                    PMTpos = PMTxyz[ipmt][:3]
                    wt = qe*att
                    COGden[gtype] += wt
                    for i in range(len(COGnum[gtype])):
                        COGnum[gtype][i] += PMTpos[i]*wt
        if report: print 'oneton.anaNtupleEvent: Event #',evtNum
        for gtype in COGnum:
            wt = COGden[gtype]
            if wt>0:
                for i,e in enumerate(COGnum[gtype]): COGnum[gtype][i] = e/wt
            if report: print 'oneton.anaNtupleEvent: gtype',gtype,'COG',COGnum[gtype],'nGamma',nGamma[gtype],'wtGamma',wt
        recDir = [COGnum[1],COGnum[0]]
        truDir = [Xfirst,Xlast]
        cosAngle = self.traj.trackAngle(recDir, truDir)
        if report: print 'oneton.anaNtupleEvent: first,last point on track',Xfirst,Xlast,'cos(recDir,truDir)',cosAngle
        if anaFile is not None:
            anaList = []
            for gtype in COGnum:
                anaList.extend(COGnum[gtype])
            for X in truDir:
                anaList.extend(X)
            for gtype in COGnum:
                anaList.append(COGden[gtype])
            for gtype in nGamma:
                anaList.append(nGamma[gtype])
            anaList.append(cosAngle)
            anaList.append(evtNum)
            anaList.append(' \n')
            s = " ".join(str(a) for a in anaList)
            anaFile.write(s)
        return
    def drawDet(self,fopen):
        '''
        add entries to ntuple to enable drawing of detector.
        Need more than just extrema if only a portion of detector is to be drawn
        '''
        iopt = -2
        ipmt, wl, qe, att = -99,-1,-1,-1
        option = 'drawDet'
        X2 = [0,0,0]
        radius, ztop, zbottom = self.Detector
        nz = 10
        dz = (ztop-zbottom)/float(nz)
        nphi = 10
        dphi = self.twopi/float(nphi)
        ## side
        for iz in range(nz):
            z = float(iz)*dz + zbottom
            for iphi in range(nphi):
                phi = float(iphi)*dphi 
                x = radius*math.cos(phi)
                y = radius*math.sin(phi)
                X1 = [x,y,z]
                track = [X1,X2]
                photon = [track, option, iopt]
                self.makeNtuple(fopen,photon,ipmt,wl,qe,att,-1)
        ## top and bottom
        for z in [ztop, zbottom]:
            for irad in range(3):
                rad = radius/float(irad+1)
                for iphi in range(nphi):
                    phi = float(iphi)*dphi 
                    x = rad*math.cos(phi)
                    y = rad*math.sin(phi)
                    X1 = [x,y,z]
                    track = [X1,X2]
                    photon = [track, option, iopt]
                    self.makeNtuple(fopen,photon,ipmt,wl,qe,att,-1)
                
        return
    def ranPosDir(self):
        '''
        return random position vector and direction unit vector.
        '''
        ranphi = random.uniform(0.,self.twopi)
        ranrad = math.sqrt( random.uniform(0.,self.Detector[0]*self.Detector[0]) )
        ranz   = random.uniform(self.Detector[2],0.)
        tStart = [ranrad*math.cos(ranphi),ranrad*math.sin(ranphi),ranz]
        tDir = [random.uniform(-1,1),random.uniform(-1,1),random.uniform(-1,1)]
        tDir = self.normVector(tDir)
        return tStart,tDir
    def standardRun(self, nE, nC, nS, prefn='photons', Save='All', mode='Simple', ioption=-1):
        '''
        Simple mode:
        create a run of nE events with nC Cerenkov photons and nS Scint. photons per event
        output ntuple file prefix prefn
        if Save is All, then put all photons in ntuple
        otherwise only record the photons that hit a PMT

        Better mode:
        
        '''
        
        unique = '_{0}'.format(datetime.datetime.now().strftime("%Y%m%d%H%M_%f"))
        filename = self.pawDir  + prefn + unique + '.nt'
        fopen = open(filename,'w')
        print 'oneton.standardRun: Opened',filename
        self.drawDet(fopen)
        xyzPMT,ZofPMT = self.placePMTs(outFile=fopen)
        for iEvt in range(nE):
            printStuff = ( iEvt<10 or iEvt%(nE/10)==0 )
            if printStuff: print 'Event#',iEvt,'Save',Save,'mode',mode,'ioption',ioption,
            
            if mode.lower()=='simple':
                photons = self.oneSimpleEvent(nCerenkov=nC, nScint=nS)
            elif mode.lower()=='better':
                particle = 'e-'
                material = 'water'
                KE = 2.28
                tStart = [self.Detector[0]-10.,0.,self.Detector[2]/2]

                msg = 'ioption='+str(ioption)
                # good for directionality check
                if ioption==1:
                    particle = 'proton'
                    KE = 2000.
                    tStart = [0.,0.,self.Detector[2]/2.]
                    tDir = [-1., 0., 0.]
                    msg += ' dir check with p'
                elif ioption==2:
                    particle = 'e-'
                    KE = 250.
                    tStart = [-self.Detector[0]/3,-self.Detector[0]/2,self.Detector[2]/2]
                    tDir = [-1.,-1.,-1.]
                    msg += ' dir check with e-'
                # random position and direction
                elif ioption==3:
                    particle = 'e-'
                    KE = 5.
                    tStart, tDir = self.ranPosDir()
                    msg += ' e- random pos,dir'
                # cosmic
                elif ioption==4:
                    particle = 'mu+'
                    KE = 10000.
                    tStart = [0., 0., 0.]
                    tDir = self.downward
                    msg += ' cosmic mu+'
                # point at top
                elif ioption in [5,6,7,8]:
                    particle = 'e+'
                    KE = 5.
                    if ioption==5: tStart = [0.5*self.Detector[0],0.5*self.Detector[0], -100.]
                    if ioption==6: tStart = [0.5*self.Detector[0],-0.5*self.Detector[0], -100.]
                    if ioption==7: tStart = [-0.5*self.Detector[0],0.5*self.Detector[0], -100.]
                    if ioption==8: tStart = [-0.5*self.Detector[0],-0.5*self.Detector[0], -100.]
                    tDir = self.upward
                    msg = 'upward at top'
                # point at bottom
                elif ioption in [9,10,11,12]:
                    particle = 'e+'
                    KE = 5.
                    if ioption== 9: tStart = [0.5*self.Detector[0],0.5*self.Detector[0], -1100.]
                    if ioption==10: tStart = [0.5*self.Detector[0],-0.5*self.Detector[0], -1100.]
                    if ioption==11: tStart = [-0.5*self.Detector[0],0.5*self.Detector[0], -1100.]
                    if ioption==12: tStart = [-0.5*self.Detector[0],-0.5*self.Detector[0], -1100.]
                    tDir = self.downward
                    msg = 'downward at bottom'
                # point at side
                elif ioption in [13,14,15]:
                    particle = 'e+'
                    KE = 5.
                    if ioption==13:
                        tStart = [0.5*self.Detector[0],0.5*self.Detector[0], -555.]
                        tDir   = [1., 1., 0.1]
                    if ioption==14:
                        tStart = [-0.5*self.Detector[0],-0.5*self.Detector[0], -555.]
                        tDir   = [-1.,-1.,-0.1]
                    if ioption==15:
                        tStart = [-0.25*self.Detector[0],-0.5*self.Detector[0], -555.]
                        tDir   = [-.5,-1.,-0.]
                    msg = 'outward at side'
                # beta spectrum
                # Tmax of 2.2801 corresponds to 90Y39 decay (90Sr daughter)
                # Tmax of 3.541  corresponds to 106Rh decay (106Ru daughter)
                elif ioption in [16,17]:
                    particle = 'e-'
                    KE = 0.
                    isotope = 'NONE'
                    if ioption==16: 
                        Tmax = 2.2801
                        isotope = '90Sr'
                    if ioption==17:
                        Tmax = 3.541
                        isotope = '106Ru'
                    KE = self.genBeta(Tmax)
                    tStart = [0.,25.,self.Detector[2]+100.]
                    tDir = self.downward
                    msg = isotope + ' near bottom, pointed down'
                # fake elastic scatter electrons
                elif ioption in [18,19,20,21,22]:
                    particle = 'e-'
                    KE = 4.
                    r = 0.45*self.Detector[0]*float(ioption-20)/math.sqrt(2)
                    tStart = [r,r,-555.]
                    tDir = [-1., -1., -0.1]
                    msg = 'towards -ve x,y'
                
                processes = []
                if nC>0: processes.append('Cerenkov')
                if nS>0: processes.append('Scint')

                msg += ' '+particle+' in '+material+' KE '
                msg += '{0:.3f}'.format(KE)
                msg += ' Pos'+self.niceString(tStart)+' & Dir '+self.niceString(tDir)+' '+self.niceString(processes)
                if printStuff:
                    print msg
                    sys.stdout.flush()
                energys,photons = self.oneEvent(particle, material, KE, tStart, tDir=tDir,processes=processes)
                
            for j,photon in enumerate(photons):
                option = photon[1]
                
                ipmt = self.hitPMT( photon[0], xyzPMT, ZofPMT)

                if Save=='All' or ipmt>-1 :
                    if mode.lower()=='simple':
                        wl = self.Photons.getWavelength( beta=1., medium='water', process=photon[1])
                    elif mode.lower()=='better':
                        wl = self.hc/energys[j]
                    qe = self.Photons.getQE(wl)
                    length = self.traj.mag( self.traj.makeVec( photon[0] ) )
                    att = self.Photons.getAtten( wl, length, medium='water')
                    self.makeNtuple(fopen, photon, ipmt, wl, qe, att, iEvt)
        fopen.close()
        print 'exec oneton#valid',filename,' ! ','\'',msg,'\''
        
        # write a line to generation log
        with open(self.generationLog,'a') as gFile:
            gFile.write(filename + ' ioption='+str(ioption)+ ' ' + msg)
            
        # do some analysis
        self.readNtuple(filename)
        return
    def control(self,nC=1,nS=1,nE=1,ioption=1,Save='All',mode='better',det='oneton'):
        '''
        consolidate control into a module
        '''
        if len(sys.argv)>1: nC = int(sys.argv[1])
        if len(sys.argv)>2: nS = int(sys.argv[2])
        if len(sys.argv)>3: nE = int(sys.argv[3])
        if len(sys.argv)>4: Save = sys.argv[4]
        if len(sys.argv)>5: ioption = int(sys.argv[5])
        if len(sys.argv)>6: det = sys.argv[6]
        self.defineDetector(detName=det)
        self.standardRun(nE,nC,nS,Save=Save,mode=mode,ioption=ioption)
        return
if __name__ == '__main__' :
    if len(sys.argv)>1:
        if 'help' in sys.argv[1].lower():
            print 'oneton.py nCerenkov nScint nEvents Save(All,Hits) ioption oneton/dayabay'
            sys.exit()

    ton = oneton()
    ton.control()
    print ''
    #ton.standardRun(nE,nC,nS,Save=Save,mode=mode,ioption=ioption)

