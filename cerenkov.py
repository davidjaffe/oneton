#!/usr/bin/env python
'''
calculate Cerenkov yield
units: particle energy in MeV, pathlenght in cm, wavelength in nm, optical photon energy in nm
notation: OP = Optical Photon
20140723
'''
import sys
import math
import ProcessMaterialProperty
import ParticleProperties
import numpy


class cerenkov():
    def __init__(self):
        self.hc = 1239.841930 # nm per eV
        self.cerconst = 369.907 # cerenkov OP/eV/cm. Eqn 30.43 RPP2012
        self.pmp = ProcessMaterialProperty.ProcessMaterialProperty() # for IOR
        self.pp  = ParticleProperties.ParticleProperties() # for masses
        return
    def getWavelength(self,E):
        return self.hc/E
    def getIOR(self,OPenergy,material):
        '''
        get index of refraction given OP energy in eV and material
        '''
        wl = self.getWavelength(OPenergy)
        IOR = self.pmp.getMaterialProperty(wl,material=material,property='RINDEX')
        return IOR
    def getBeta(self,KE,particle):
        m = self.pp.getMass(particle)
        T = max(0.,KE)
        return math.sqrt(T*(T+2*m))/(m+T)
    def getCerenkovYield(self,KE,pathlength,particle,material,OPE1=1.7,OPE2=6.2,OPdE=0.01,makeOP=False):
        '''
        return total number of cerenkov photons, valid energy range
        given energy range (OPE1,OPE2), particle in material with KE and pathlength
        returned energy range assumes index of refraction is monotonically increasing
        or decreasing in input energy range (OPE1,OPE2).
        if makeOP, then generate cerenkov photons and return a list of energies of the generated
        optical photons with corresponding cos(theta_Cerenkov)
        '''
        beta = self.getBeta(KE,particle)
        prefactor = max(0.,pathlength) * self.cerconst
        OPE = OPE1
        OPEvec = []
        OpticalPhotons, OPcost = [], []
        rsum = 0.
        while OPE<=OPE2:
            n = self.getIOR(OPE,material)
            r = (1. - 1./beta/beta/n/n)
            if r>0.:
                rsum += r*OPdE
                OPEvec.append(OPE)
                if makeOP:
                    nOP = numpy.random.poisson(prefactor * r * OPdE)
                    for iOP in range(nOP):
                        e = numpy.random.uniform(OPE,OPE+OPdE)
                        OpticalPhotons.append( e )
                        cost = 1./beta/self.getIOR(e,material)
                        OPcost.append( cost )
            OPE += OPdE
        yld = prefactor * rsum
        OPElimits = None
        if len(OPEvec)>0: OPElimits = [min(OPEvec), max(OPEvec)]
        
        return yld, OPElimits, OpticalPhotons, OPcost
    def printCerenkovPhotons(self, OpticalPhotons, OPcost, binned=False,binsize=0.1):
        # sort lists based on wavelength
        # http://stackoverflow.com/questions/24152993/how-to-sort-one-item-to-match-another-but-applying-the-sort-elsewhere
        # cost = -2. signals a scintillation photon
        OpticalPhotons,OPcost = zip(*sorted(zip(OpticalPhotons,OPcost)))
        if binned:
            Elo = min(OpticalPhotons) - binsize/2.
            Ehi = max(OpticalPhotons)
            nOP, avect, aveE = 0, 0., 0. # cerenkov
            nSci, aveSciE = 0, 0. # scintillation
            print '\n{0:>6} {1:>6} {2:>6} or Scintillation'.format('nOP','<E> eV','<cost>')
            for E,ct in zip(OpticalPhotons,OPcost):
                if E>Elo+binsize or E>Ehi:
                    if nOP>0:
                        avect = avect/float(nOP)
                        aveE = aveE/float(nOP)
                        print '{1:6d} {0:6.4f} {2:6.4f}'.format( aveE,nOP,avect )
                    if nSci>0:
                        aveSciE = aveSciE/float(nSci)
                        print '{1:6d} {0:6.4f} {2:>6}'.format( aveSciE,nSci, 'Scint')
                    Elo += binsize
                    nOP,avect,aveE = 0, 0., 0.
                    nSci, aveSciE = 0, 0.
                if ct>-2.:
                    nOP += 1
                    avect += ct
                    aveE += E
                else:
                    nSci += 1
                    aveSciE += E
            if nOP>0: # flush cerenkov photons
                avect = avect/float(nOP)
                aveE = aveE/float(nOP)
                print '{1:6d} {0:6.4f} {2:6.4f}'.format( aveE,nOP,avect )
            if nSci>0: # flush scintillation photons
                aveSciE = aveSciE/float(nSci)
                print '{1:6d} {0:6.4f} {2:>6}'.format( aveSciE,nSci, 'Scint')
                
        else:
            print '\n{0:>6} {1:>6} {2:>6} or Scintillation'.format('i','Ei(eV)','cost_i')
            for i,pair in enumerate(zip(OpticalPhotons,OPcost)):
                print '{0:6d} {1:6.4f} {2:6.4f}'.format(i,pair[0],pair[1])
        return 
        
if __name__ == '__main__':
    c = cerenkov()
    KE = 1.
    pathlength = 1.0
    particle = 'e+'
    material = 'water'
    OPE2 = c.hc/300.
    OPE1 = c.hc/600.
    y,OPElimits,OpticalPhotons,OPcost = c.getCerenkovYield(KE,pathlength,particle,material,OPE1=OPE1,OPE2=OPE2,makeOP=True)
    print '\nyield',y,'OPE energy range(eV)',OPElimits,'Ngen(OP)',len(OpticalPhotons)
    c.printCerenkovPhotons(OpticalPhotons, OPcost)
    c.printCerenkovPhotons(OpticalPhotons, OPcost, binned=True,binsize=0.25)
    
    if 0: # check of calculations
        beta = c.getBeta(KE,particle)
        n1 = c.getIOR(OPE1,material)
        n2 = c.getIOR(OPE2,material)
        n3 = 1.34
        y1 = c.cerconst * pathlength * (OPE2-OPE1) * (1. - 1./beta/beta/n1/n1)
        y2 = c.cerconst * pathlength * (OPE2-OPE1) * (1. - 1./beta/beta/n2/n2)
        y3 = c.cerconst * pathlength * (OPE2-OPE1) * (1. - 1./beta/beta/n3/n3)
        print 'beta',beta,'y1',y1,'n1',n1,'y2',y2,'n2',n2,'y3',y3,'n3',n3,'cerconst',c.cerconst,'pathlength',pathlength

    if 0:
        import Range
        r = Range.Range()
        ke = 2.283
        t  = 10.
        totR, finalKE, samples = r.propagate('electron','water',ke,t)
        print '\n initial KE(MeV)',ke,'thickness(cm)',t,'total range(cm)',totR,'final KE(MeV)',finalKE
        r.printSampleTable(samples)
        totYield = 0.
        totpath = 0.
        totCpath = 0.
        for sample in samples:
            KE,pathlength = sample
            y,OPElimits = c.getCerenkovYield(KE,pathlength,'electron','water',OPE1=OPE1,OPE2=OPE2)
            totYield += y
            totpath += pathlength
            if OPElimits is not None: totCpath += pathlength
        print 'total cerenkov yield',totYield,'total pathlength',totpath,'total pathlength with Cerenkov',totCpath
