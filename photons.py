#!/usr/bin/env python
'''
generation of optical photons for scintillation and cerenkov processes
uses optical properties from daya bay via ProcessMaterialProperty
20140102
'''
import math
import sys
import random
import ProcessMaterialProperty

class photons():
    def __init__(self, wl1=200., wl2=800.):
        self.MatProp = {}
        self.CerenkovWLint = {}
        m = 'scintillator'
        for p in ['FASTCOMPONENT','SLOWCOMPONENT','RINDEX','ABSLENGTH']:
            k = m.lower() + '_' + p.lower()
            self.MatProp[k] = [m,p]

        for m in ['water', 'acrylic']:
            for p in ['RINDEX', 'ABSLENGTH']:
                k = m.lower() + '_' + p.lower()
                self.MatProp[k] = [m,p]

        # PMT QE
        m = 'bialkali'
        p = 'EFFICIENCY'
        k = m.lower() + '_' + p.lower()
        self.MatProp[k] = [m,p]
        

        self.pmp = ProcessMaterialProperty.ProcessMaterialProperty()
        for k in self.MatProp:
            m = self.MatProp[k][0]
            p = self.MatProp[k][1]
            self.pmp.iniRandom(material=m, property=p)

        self.wlRange = [wl1, wl2, 1.] # min,max, step size
            
        return
    def getQE(self, wl, material='bialkali'):
        QE = self.pmp.getMaterialProperty(wl, material=material, property='EFFICIENCY')
        return QE
    def getAtten(self, wl, pathL, medium='water'):
        '''
        compute attenuation of light in medium given wavelength wl and
        photon path length pathL in mm
        '''
        lam = self.pmp.getMaterialProperty(wl, material=medium, property='ABSLENGTH')
        A = math.exp(-pathL/lam)
        #print 'photons.getAtten wl',wl,'pathL',pathL,'attLen',lam,'attenuation',A
        return A
    def transFraction(self, wl, gVec, nVec, med1='water', med2='acrylic', process='Cerenkov'):
        '''
        compute the transmittance of light of wavelength wl going from med1 to med2
        with incident light vector gVec and normal vector nVec between media
        take into account polarization of Cerenkov light
        '''
        n1 = self.pmp.getMaterialProperty(wl,material=med1,property='RINDEX')
        n2 = self.pmp.getMaterialProperty(wl,material=med2,property='RINDEX')
        gn = nn = gg = 0
        for g,n in zip(gVec,nVec):
            gn += g*n
            nn += n*n
            gg += g*g
        cosi = gn/math.sqrt( gg * nn)
        sini = math.sqrt(1. - cosi*cosi)
        #print 'n1,n2,cosi,sini',n1,n2,cosi,sini
        f = n1/n2*sini
        if abs(f)>=1.: return 0. # total internal reflection
        cost = math.sqrt(1. - f*f)
        Rs = (n1*cosi - n2*cost) / (n1*cosi + n2*cost) # polarized perpendicular to plane of incidence
        Rs = Rs*Rs
        Rp = (n1*cost - n2*cosi) / (n1*cost + n2*cosi) # polarized parallel to plane of incidence
        Rp = Rp*Rp
        #print 'cost,Rs,Rp',cost,Rs,Rp
        if process.lower()=='cerenkov':
            T = 1. - Rp
        else:
            T = 1. - (Rp+Rs)/2.
        return T
        
    def getWavelength(self, beta=1., process='Cerenkov', medium='water'):
        if process.lower()=='cerenkov':
            wl = self.drawCerenkov(beta=beta, medium=medium)
        elif 'scint' in process.lower():
            wl = self.pmp.drawRandom(property='FASTCOMPONENT', material='scintillator')
        else:
            w = 'photons.getWavelength ERROR invalid option '+option
            sys.exit(w)
        return wl
    def drawCerenkov(self,beta=1., medium='water'):
        '''
        return random wavelength wl for particle with velocity beta in medium
        '''
        ik = self.initCerenkov(beta=beta, medium=medium)
        aboveThreshold = self.CerenkovWLint[ik][4]
        if not aboveThreshold: return -1 # not above cerenkov threshold
        
        wavelengths = self.CerenkovWLint[ik][2]
        cumint      = self.CerenkovWLint[ik][3]
        v = random.uniform( cumint[0], cumint[-1] )
        j = min(range(len(cumint)), key=lambda i: abs(cumint[i]-v))
        if v<cumint[j]:
            k = max(0,j-1)
        else:
            k = min(len(cumint)-1,j+1)
        return self.pmp.yint(cumint[j],wavelengths[j], cumint[k],wavelengths[k], v)
        
    def initCerenkov(self,beta=1., medium='water'):
        # check if initialization has already been done
        for ik in self.CerenkovWLint:
            l = self.CerenkovWLint[ik]
            if l[0]==beta and l[1]==medium: return ik # already done

        # 
        wl1,wl2,dwl = self.wlRange

        wl = wl1
        wavelengths = []
        cumint = []
        rsum = 0.
        aboveThreshold = False
        
        while wl<=wl2:
            n = self.pmp.getMaterialProperty(wl, material=medium, property='RINDEX')
            r = 1./wl/wl * (1. - 1./beta/beta/n/n)
            if r>0. : aboveThreshold = True
            rsum += max(0.,r)
            cumint.append(rsum)
            wavelengths.append(wl)
            wl += dwl

        ik = len(self.CerenkovWLint)+1
        self.CerenkovWLint[ik] = [beta,medium, wavelengths,cumint, aboveThreshold]
        print 'photons.initCerenkov key#',ik,'beta',beta,'medium',medium
        return ik

if __name__ == '__main__' :
    ph = photons()
    print ''
    for ycomp in range(0,100,10):
        gVec = [0,float(ycomp)/10.,1]
        nVec = [0,0,1]
        nn = gg = gn = 0
        for g,n in zip(gVec,nVec):
            gg += g*g
            nn += n*n
            gn += g*n
        cosi = gn/math.sqrt(gg*nn)
        a = math.acos(cosi)*180./math.pi
        angle = '{0:.3f}'.format(a)

        pathL = 6.

        
        for wl in range(200,900,300):
            for med1 in ['water']: #, 'scintillator']:
                med2 = 'acrylic'
                Tc = ph.transFraction(wl,gVec,nVec,process='Cerenkov',med1=med1,med2=med2)
                Ts = ph.transFraction(wl,gVec,nVec,process='Scint',med1=med1,med2=med2)
                A  = ph.getAtten(wl, pathL, medium=med1)
                qe = '{0:.4f}'.format(ph.getQE(wl,material='bialkali'))
                qeR7723 = '{0:.4f}'.format(ph.getQE(wl,material='QE_R7723'))
                print 'wl',wl,med1,'to',med2,'angle',angle,'Tc',Tc,'Ts',Ts,'pathL',pathL,'attenuation',A,'QE',qe,'QE_R7723',qeR7723

    print ''
    for wl in range(200,800,50):
        qe = '{0:.4f}'.format(ph.getQE(wl,material='bialkali'))
        qeR7723 = '{0:.4f}'.format(ph.getQE(wl,material='QE_R7723'))
        print 'wl',wl,'QE',qe,'QE_R7723',qeR7723
            
