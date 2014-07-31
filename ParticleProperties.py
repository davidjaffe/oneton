#!/usr/bin/env python
'''
particle properties
masses in MeV
lifetimes in ns (requires more code)
20140723
'''
import sys
import os

class ParticleProperties():
    def __init__(self):
        self.mass = {'K0': 497.614,\
                     'K+': 493.677,\
                     'lambda' : 1115.683,\
                     'proton' : 938.272046,\
                     'neutron': 939.565379,\
                     'pion'   : 139.57018,\
                     'pi0'    : 134.9766,\
                     'electron': 0.510998928,\
                     'muon'   : 105.6583715}
        self.aliases = {'K0'    : ['K0','K0L','K0S'],\
                        'K+'    : ['K+','K-'],\
                        'lambda': ['lambda'],\
                        'proton': ['proton','antiproton','p','pbar'],\
                        'neutron':['neutron','n'],\
                        'pion'  : ['pion','pi+','pi-'],\
                        'pi0'   : ['pi0','pizero'],\
                        'electron': ['electron','e-','positron','e+'],\
                        'muon'  : ['muon','mu+','mu-'] }
        return
    def getNameFromAlias(self,alias):
        for name in self.aliases:
            if alias.lower() in self.aliases[name] : return name
            for a in self.aliases[name]:
                if alias.lower()==a.lower() : return name
        return None
    def getMass(self,particle):
        name = self.getNameFromAlias(particle)
        if name is None:
            words = 'ParticleProperties.getMass: ERROR Invalid particle ' + str(particle)
            sys.exit(words)
        return self.mass[name]

if __name__ == '__main__':
    pp = ParticleProperties()
    print 'e- mass',pp.getMass('e-')
    print 'dog mass',pp.getMass('dog')
