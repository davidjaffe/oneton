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

class oneton():
    def __init__(self):

        self.debug = 0

        # initialize external classes. Note self.Photons initialized below
        self.traj = trajectory.trajectory()
        self.Range = Range.Range()
        self.cerenkov = cerenkov.cerenkov()

        # output directory for ntuple
        self.pawDir = '/Users/djaffe/work/paw/ONETON/'



        # construct detector
        # coordinate system has z=0 at top, positive z pointing up
        ID = 995.
        IH = 1250.
        thick = 25.4
        self.Detector = [ID/2., 0., -IH] # [radius,ztop, zbottom]

        c = 1.
        self.HodoWidth = c
        self.Hodoscope = [-c/2., c/2, -c/2, c/2] # x1,y1 , x2,y2 = edges of hodoscope

        self.PMTradius = 25.4

        self.ior = 1.33  # index of refraction

        self.dEdx = 0.2 # MeV/mm
        self.ScintYield = 100 # photons/MeV

        #self.wlRange = [300., 600. ] # range of wavelengths for cerenkov light relevant for bialkalai pmt (nm)
        self.wlRange = [250., 700. ] # wider range of wavelengths for cerenkov light relevant for bialkalai pmt (nm)

        self.hc = 2.*math.pi*197.326972 # eV*nm or MeV*fm
        self.Photons = photons.photons(wl1=self.wlRange[0], wl2=self.wlRange[1] )

        self.downward = [0., 0., -1.] # unit vector pointing down
        self.upward   = [0., 0.,  1.] # unit vector pointing up
        self.twopi = 2.*math.pi
        self.options = {'Cerenkov':0, 'Scint':1}

        self.savedEvents = {}
        self.savedEventKeys = {}
        return
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
            dist = (z1-z2)/unitV[2]
            x2 = unitV[0]*dist + x1
            y2 = unitV[1]*dist + y1
            vector = [ [x1,y1,z1], [x2,y2,z2] ]
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
        if abs(numpy.dot(t,g)-cost)>1e-10: g = -g
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
        print 'oneton.oneEvent:',particle,'in',material,'KE(MeV)',KE,'Initial pos',self.niceString(tStart),' & Direction',self.niceString(tDir),'processes',self.niceString(processes)
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
                tDir, track, totR, finalKE, samples = self.savedEvents[eKey]
                if self.debug>2: print 'oneton.oneEvent: retrieved track',track

        if not found:
            ### make sure direction vector is unit vector
            tDir = self.normVector( tDir )

            ### determine material thickness from start point given starting direction
            track = self.makeVector( tStart, tDir)
            thickness = self.traj.mag( self.traj.makeVec( track ) )

            ### determine total range along initial trajectory (no scattering)
            totR, finalKE, samples = self.Range.propagate(particle, material, KE, thickness, maxdEperSample=0.1)
            self.savedEvents[eKey] = tDir, copy.deepcopy( track ), totR, finalKE, samples
            if self.debug>2: print 'oneton.oneEvent: saving track',track
            if self.debug>0: 
                print '\n totR',totR,'finalKE',finalKE,'thickness',thickness
                if self.debug>1: self.Range.printSampleTable( samples )

        ### now create track with start and endpoints
        X1, X2 = track
        for i,u in enumerate(tDir):
            X2[i] = totR*u + X1[i]
        track = [ X1, X2]
        if self.debug>2: print 'oneton.oneEvent track',track

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
    def placePMTs(self):
        radcyl, ztop, zbottom = self.Detector
        xyzPMT = []
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
        return xyzPMT
    def hitPMT(self,photonT,xyzPMT):
        '''
        determine hit PMT for input photon track photonT = [ [startPoint], [endPoint] ]
        given PMT positions in xyzPMT
        error message if >1 hit PMT
        '''
        hits = []
        nhit = 0

        X1,X2 = photonT
        x,y,z = X2
        for pmt in xyzPMT:
            xPMT, yPMT, zPMT, rPMT, side = pmt
            hit = 0
            dz = 0.00001
            if side : dz = rPMT
            if abs(z-zPMT)<dz:
                if math.sqrt( (x-xPMT)*(x-xPMT) + (y-yPMT)*(y-yPMT) + (z-zPMT)*(z-zPMT) )<rPMT:
                       hit = 1
            hits.append(hit)
            nhit += hit
        if nhit>1:
            print 'oneton.hitPMT ERROR too many hit PMTs. NPMT',nhit
            print 'photonT',photonT
            print 'hits',hits
            for hit,pmt in zip(hits,xyzPMT):
                if hit>0: print 'hit pmt',pmt
        return hits
    def makeNtuple(self,fopen, photon, ipmt, wl, qe, att):
        '''
        for each photon   start     end       0/1   
        ntuple variables: x1,y1,z1, x2,y2,z2, iopt, ipmt, wl, qe, attenuation factor
        iopt = 0 = Cerenkov, 1 = Scintillation, -2 = drawDet
        '''
        track, option, iopt = photon
        X1, X2 = track
        s = ''
        for u in X1: s += str(u) + ' '
        for u in X2: s += str(u) + ' '
        s += str(iopt) + ' ' 
        s += str(ipmt) + ' ' 
        s += str(wl)   + ' '
        s += str(qe)   + ' '
        s += str(att) 
        s += ' \n'
        fopen.write(s)
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
        for iz in range(nz):
            z = float(iz)*dz + zbottom
            for iphi in range(nphi):
                phi = float(iphi)*dphi 
                x = radius*math.cos(phi)
                y = radius*math.sin(phi)
                X1 = [x,y,z]
                track = [X1,X2]
                photon = [track, option, iopt]
                self.makeNtuple(fopen,photon,ipmt,wl,qe,att)
        return
    def ranPosDir(self):
        '''
        return random position vector and direction unit vector
        '''
        ranphi = random.uniform(0.,self.twopi)
        ranrad = random.uniform(0.,self.Detector[0])
        ranz   = random.uniform(self.Detector[2],0.)
        tStart = [ranrad*math.cos(ranphi),ranrad*math.sin(ranphi),ranz]
        tDir = [random.uniform(-1,1),random.uniform(-1,1),random.uniform(-1,1)]
        tDir = self.normVector(tDir)
        return tStart,tDir
    def standardRun(self, nE, nC, nS, prefn='photons', Save='All', mode='Simple'):
        '''
        Simple mode:
        create a run of nE events with nC Cerenkov photons and nS Scint. photons per event
        output ntuple file prefix prefn
        if Save is All, then put all photons in ntuple
        otherwise only record the photons that hit a PMT

        Better mode:
        
        '''
        
        xyzPMT = self.placePMTs()
        unique = '_{0}'.format(datetime.datetime.now().strftime("%Y%m%d%H%M"))
        filename = self.pawDir  + prefn + unique + '.nt'
        fopen = open(filename,'w')
        print 'oneton.standardRun: Opened',filename
        self.drawDet(fopen)
        for iEvt in range(nE):
            print 'Event#',iEvt,'Save',Save,'mode',mode
            if mode.lower()=='simple':
                photons = self.oneSimpleEvent(nCerenkov=nC, nScint=nS)
            elif mode.lower()=='better':
                particle = 'e-'
                material = 'water'
                KE = 2.28
                #particle = 'proton'
                #KE = 475. # 2000.
                tStart = [0., 0., -1250+12.]
                tDir = [.0, 0.0, -1.]
                #tStart = [0., 0., -12.]
                tStart = [self.Detector[0]-10.,0.,self.Detector[2]/2]

                # good for directionality check
                if 0:
                    particle = 'proton'
                    KE = 2000.
                    tStart = [0.,0.,self.Detector[2]/2.]
                    tDir = [-1., 0., 0.]
                if 1:
                    particle = 'e-'
                    KE = 250.
                    tStart = [-self.Detector[0]/3,-self.Detector[0]/2,self.Detector[2]/2]
                    tDir = [-1.,-1.,-1.]

                # random position and direction
                if 0:
                    particle = 'e-'
                    KE = 150.
                    tStart, tDir = self.ranPosDir()

                # cosmic
                if 0:
                    particle = 'mu+'
                    KE = 10000.
                    tStart = [0., 0., 0.]
                    tDir = self.downward


                
                processes = []
                if nC>0: processes.append('Cerenkov')
                if nS>0: processes.append('Scint')
                energys,photons = self.oneEvent(particle, material, KE, tStart, tDir=tDir,processes=processes)
                
            for j,photon in enumerate(photons):
                option = photon[1]
                
                hits = self.hitPMT( photon[0], xyzPMT)
                ipmt = -1
                if 1 in hits :
                    ipmt = hits.index(1)
                if Save=='All' or ipmt>-1 :
                    if mode.lower()=='simple':
                        wl = self.Photons.getWavelength( beta=1., medium='water', process=photon[1])
                    elif mode.lower()=='better':
                        wl = self.hc/energys[j]
                    qe = self.Photons.getQE(wl)
                    length = self.traj.mag( self.traj.makeVec( photon[0] ) )
                    att = self.Photons.getAtten( wl, length, medium='water')
                    self.makeNtuple(fopen, photon, ipmt, wl, qe, att)
        fopen.close()
        return
if __name__ == '__main__' :
    ton = oneton()
    if 0:
        particle = 'e-'
        material = 'water'
        KE = 2.28
        tStart = [0., 0., -1230.]
        tDir = [0., 0., -1.]
        processes = ['Scint','Cerenkov']
        OPenergys,OPtracks = ton.oneEvent(particle, material, KE, tStart, tDir=tDir,processes=processes)
        for E,track in zip(OPenergys,OPtracks):
            print E,track
    
    if 1: 
        print ''
        nC = 1
        nS = 1
        nE = 1
        Save = 'All'
        mode = 'Better'
        if len(sys.argv)>1: nC = int(sys.argv[1])
        if len(sys.argv)>2: nS = int(sys.argv[2])
        if len(sys.argv)>3: nE = int(sys.argv[3])
        if len(sys.argv)>4: Save = sys.argv[4]
        ton.standardRun(nE,nC,nS,Save=Save,mode=mode)

