#!/usr/bin/env python
'''
calculate particle range in a material
units: MeV, cm, g
20140722
'''
import sys
import os

class Range():
    def __init__(self):
        # where is the range, dE/dx data?
        self.ESTARdir = 'ESTAR/'
        self.PSTARdir = 'PSTAR/'
        self.MUSTARdir= 'MUONSTOP/muonloss_'
        self.dirs = {self.ESTARdir: ['electron', 'positron', 'e-', 'e+'],\
                     self.PSTARdir: ['proton', 'antiproton', 'p', 'pbar'],\
                     self.MUSTARdir:['muon','antimuon','mu+','mu-']\
                     }
        # what is the meaning of each column in the data files
        self.columns = {self.ESTARdir: [ {'KE':0}, {'TotalStoppingPower': 3}, {'CSDA Range': 4} ],\
                        self.PSTARdir: [ {'KE':0}, {'TotalStoppingPower': 3}, {'CSDA Range': 4} ],\
                        self.MUSTARdir:[ {'KE':2}, {'TotalStoppingPower': 9}, {'CSDA Range': 10} ]\
                            }
        # map of material name to aliases, valid particle types and density (g/cm3)
        self.Materials = {'water' : [ ['liquidwater', 'water'] , ['electron','positron'], 1.00 ],\
                          'polyethylene' : [ ['polyethylene'], ['electron','positron'], 0.97 ], \
                          'toluene': [ ['toluene'], ['electron', 'positron'], 8.66900E-01 ] }
        
        self.tiny = 1.e-16

        self.Table = {}
        self.KEtable = {}
        self.dEdxtable={}
        self.CSDAtable={}
        self.TableKeys = {}
        self.KeyNumber = 0

        return
    def getParentMaterial(self,matAlias):
        '''
        return parent material corresponding to material name alias
        return null if invalid alias given
        '''
        for mat in self.Materials:
            aliases = self.Materials[mat][0]
            for a in aliases:
                if a.lower()==matAlias.lower(): return mat
        return None
    def genKey(self,particle,matAlias):
        '''
        generate or retrieve key associated with particle and material name alias
        exit if invalid material alias name given
        '''
        mat = self.getParentMaterial(matAlias)
        if mat is not None:
            pair = [particle,mat]
            for key in self.TableKeys:
                if pair==self.TableKeys[key]:
                    return key
            self.KeyNumber += 1
            self.TableKeys[self.KeyNumber] = pair
            return self.KeyNumber
        else:
            words = 'range.genKey: Invalid material alias name ' + matAlias
            sys.exit(words)
        return None # should never get here
    def getDir(self,particle):
        '''
        return directory path for input particle
        '''
        for d in self.dirs:
            for p in self.dirs[d]:
                if particle.lower()==p.lower(): return d
        return None
    def getFilename(self,particle,matAlias):
        '''
        return the filename associated with the particle and material alias name
        '''
        directory = self.getDir(particle)
        mat = self.getParentMaterial(matAlias)
        for name in self.Materials[mat][0]:
            fn = directory + name + '.txt'
            if os.path.isfile(fn): return fn
        return None
    def getTable(self,particle,matAlias):
        '''
        return the key for the dict corresponding to the table of dE/dx and range.
        generate the key and fill dict with ESTAR, PSTAR, ASTAR, etc. table if necessary
        '''
        key = self.genKey(particle,matAlias)
        if key not in self.Table:
            self.Table[key] = []
            self.KEtable[key] = []
            self.dEdxtable[key] = []
            self.CSDAtable[key] = []
            mat = self.getParentMaterial(matAlias)
            density = self.Materials[mat][2]
            fn = self.getFilename(particle,matAlias)
            cols = self.columns[self.getDir(particle)]
            iKE, iStoppingPower, iRange = -1,-1,-1
            for d in cols:
                if 'KE' in d: iKE = d['KE']
                if 'TotalStoppingPower' in d: iStoppingPower = d['TotalStoppingPower']
                if 'CSDA Range' in d: iRange = d['CSDA Range']
            f = open(fn,'r')
            ignoredLines = ''
            for line in f:
                if line[0] is not '*':
                    a = line.split(' ')
                    #print line[:-1]
                    #print a
                    KE = float(a[iKE])
                    # convert stopping power to MeV/cm and range to cm
                    StoppingPower = float(a[iStoppingPower]) * density
                    try:
                        Range = float(a[iRange])/density
                    except:
                        ignoredLines += 'range.getTable Ignored '+line
                    else:
                        self.Table[key].append( [KE, StoppingPower, Range] )
                        self.KEtable[key].append( KE )
                        self.dEdxtable[key].append( StoppingPower )
                        self.CSDAtable[key].append( Range )
            f.close()
            print 'range.getTable for particle',particle,'& material',matAlias,'with key',key,'from filename',fn
            if len(ignoredLines)>0: print ignoredLines
        else:
            pass
            #print 'range.getTable key',key,'exists for particle',particle,'& material',matAlias
        return key
    def setStepSize(self,KE):
        '''
        set the step size based on particle kinetic energy
        '''
        step = 0.1 # 1mm steps for table
        if KE < 100. : step = 0.01
        if KE < 10.  : step = 0.001
        if KE < 1.   : step = 0.0001
        return step
    def xint(self,x1,y1,x2,y2,yt):
        dx = x2-x1
        if dx==0. : return x1
        m = (y2-y1)/dx
        if m==0. : return (y1+y2)/2.
        b = y2 - m*x2
        return (yt-b)/m
    def getdEdx(self,key,KE):
        '''
        return dE/dx and maximum step size given table key and kinetic energy

        Special treatment at both ends of table,
        otherwise linearly interpolate using neighboring entries
        '''
        dEdx, maxstep = None, None
        if KE<self.KEtable[key][0]:
            dEdx = self.dEdxtable[key][0]
            dist = KE/dEdx
            maxstep = min(dist,self.CSDAtable[key][0])
        elif KE>=self.KEtable[key][-1]:
            dEdx = self.dEdxtable[key][-1]
            maxstep = self.CSDAtable[key][-1]
        else:
            for i,Ei in enumerate(self.KEtable[key]):
                j = i+1
                Ej = self.KEtable[key][j]
                if Ei<=KE and KE<Ej:
                    dEdxi = self.dEdxtable[key][i]
                    dEdxj = self.dEdxtable[key][j]
                    dEdx = self.xint(dEdxi,Ei,dEdxj,Ej,KE)
                    maxstep = min(self.CSDAtable[key][i], self.CSDAtable[key][j] )
                    return dEdx,maxstep
        return dEdx,maxstep
        
    def propagate(self,particle,matAlias,initialKE,thickness,nomStep=-1,nomSample=0.1,maxdEperSample=0.2):
        '''
        propagate particle thru material given thickness and initial kinetic energy
        return total range in material, final kinetic energy and a list of [T_i, dx_i]
        where T_i is the kinetic energy at the beginning of the ith sample and
        dx_i is the ith sample size
        The nominal sample size is nomSample unless the energy loss per sample exceeds maxdEperSample.
        (The defaults of nomSample and maxdEperSample are set such that maxdEperSample takes
        over in the region where dE/dx ~ 1/beta^2; that is, particle energies below minimum ionizing.)
        The nominal step size is nomStep unless nomStep<0, then the step size is
        determined dynamically based on the kinetic energy.
        '''
        key = self.getTable(particle,matAlias)
        KE = initialKE
        thick = thickness
        samples = []
        while KE>0. and thick>0.:
            step = nomStep
            if step<0: step = min( self.setStepSize(KE), nomSample, thick)
            sam = 0.
            dEperSample = 0.
            startKE = KE
            # make tiny adjustment to deal with numerical precision
            while KE>0. and sam<nomSample-self.tiny and sam<thick and dEperSample<maxdEperSample:
                dedx,maxstep = self.getdEdx(key,KE)
                step = min(maxstep,step)
                eloss = dedx*step
                dEperSample += eloss
                KE -= eloss
                sam += step
            samples.append( [startKE, sam] )
            thick -= sam
        totalRange = thickness-thick
        finalKE = KE
        return totalRange, finalKE, samples
    def printSampleTable(self,samples):
        '''
        print samples in a table with nice format
        '''
        print '{0:>5} {1:>6} {2:>6} {3:>6}'.format('i','KE_i','dx_i','dx_tot')
        dx_tot = 0.
        for i,pair in enumerate(samples):
            dx_tot += pair[1]
            print '{0:5d} {1:6.3f} {2:6.3f} {3:6.2f}'.format(i,pair[0],pair[1],dx_tot)
        return
if __name__ == '__main__' :
    r = Range()
    if 0:
        k = r.getTable('muon','liquidwater')
        k = r.getTable('electron','liquidwater')
        k = r.getTable('electron','water')
        k = r.getTable('electron','toluene')
        k = r.getTable('p','toluene')
    etest = 0
    if etest:
        KEs = [2.283]
        thick = [10.]
        for ke,t in zip(KEs,thick):
            totR, finalKE, samples = r.propagate('electron','water',ke,t)
            print '\n initial KE(MeV)',ke,'thickness(cm)',t,'total range(cm)',totR,'final KE(MeV)',finalKE
            r.printSampleTable(samples)
    ptest = 1
    if ptest:
        level = 1
        KEs = [198., 199., 200., 201., 204., 207., 210., 213., 216., 219.]
        thick = [4.0] #[1.0, 1.5, 2.0, 2.5]
        poly = 20.831
        for mat in ['water', 'polyethylene' ]:
            if level>0: print '\n---- holder is',thick[0],'cm',mat,' ----'
            for ke in KEs:
                if level>1: print '\n INITIAL BEAM ENERGY',ke,'MeV -----------------'
                initialKE = ke
                for t in thick:
                    totR, finalKE, samples = r.propagate('proton',mat,ke,t)
                    Eloss1 = ke-finalKE
                    if level>1: print ' inital KE(MeV)',ke,mat,'thickness(cm)',t,'total range(cm)',totR,'final KE(mEV)',finalKE,'energy lost(MeV)',Eloss1
                    newKE = finalKE
                    totR, finalKE, samples = r.propagate('proton','polyethylene',newKE,poly)
                    if level>1:print ' inital KE(MeV)',newKE,'poly thickness(cm)',poly,'total range(cm)',totR,'final KE(mEV)',finalKE,'energy lost(MeV)',newKE-finalKE
                    newKE = finalKE
                    totR, finalKE, samples = r.propagate('proton',mat,newKE,t)
                    Eloss2 = newKE-finalKE
                    if level>1: print ' inital KE(MeV)',newKE,mat,'thickness(cm)',t,'total range(cm)',totR,'final KE(mEV)',finalKE,'energy lost(MeV)',Eloss2,'ratio(up/down)=',Eloss2/Eloss1

                    if level>0: print ' Initial KE(MeV) {0:.2f}, final KE {1:.2f} ratio(Rear/Front) {2:.2f}'.format(initialKE,finalKE,Eloss2/Eloss1)
            
    mutest = 0
    if mutest:
        KEs = [2000.]
        thick = [1.]
        for ke,t in zip(KEs,thick):
            totR, finalKE, samples = r.propagate('muon','water',ke,t)
            print '\n initial KE(MeV)',ke,'thickness(cm)',t,'total range(cm)',totR,'final KE(MeV)',finalKE
            r.printSampleTable(samples)
            

