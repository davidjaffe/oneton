#!/usr/bin/env python
'''
process material property files used by DYB/DetDesc
and get property as function of wavelength in nm
20140102
'''

import math
import sys
import random

class ProcessMaterialProperty():
    defaultDir ='/Users/djaffe/work/GIT/ONETON/DYBmaterials/'# 'DYBmaterials/' # '/Users/djaffe/Documents/Reactor/work/trunk/NuWa-trunk/dybgaudi/Detector/XmlDetDesc/DDDB/materials/'
    otherDir   = '/Users/djaffe/work/GIT/ONETON/InputData/'#'./' # /Users/djaffe/PythonScripts/'
    def __init__(self,dbg=0):
        self.debug =  dbg
        self.Property = {}    # map of maps
        self.sortedPropertyKeys = {}

        self.denseProperty = {} # map of maps created for drawing randomly
        self.sortedDensePropertyKeys = {}

        self.Integral = {} 
        self.sortedIntegralKeys = {}
        self.sortedIntegralValues = {}
        return
    def getUnit(self,line):
        s = line.strip().split('=')
        unit = s[1].strip('"')
        return unit
    def setDebug(self,value):
        self.debug = value
        return
    def otherRead(self,property='EFFICIENCY', material='water'):
        fn = self.otherDir + material + '.txt'
        f = open(fn,'r')
        if self.debug>0: print 'opened',fn
        result = {}
        ndata = 0
        wlmin = 1e6
        wlmax =-1e6
        for line in f:
            if line.find('#')<0: # avoids header
                a = line.replace('\n','').split(',')
                wl,qe = float(a[0]), float(a[1])/100.  # convert QE from percent to fraction
                wlmin = min(wl,wlmin)
                wlmax = max(wl,wlmax)
                result[wl] = qe
                ndata += 1
        print 'ProcessMaterialProperty.otherRead',ndata,'lines read. range(',wlmin,',',wlmax,')nm for',property,'& material',material
        return result
    def read(self,property='ABSLENGTH',material='water'):
        '''
        first try to get info from Daya Bay xml files, then
        try alternative
        '''
        fn = self.defaultDir + material.lower() + '.xml'
        try: 
            f = open(fn,'r')
        except IOError:
            return self.otherRead(property=property, material=material)
        
        if self.debug>0: print 'opened',fn
        foundType = False
        ready = False
        done  = False
        ndata = 0
        wlmin = 1e6
        wlmax =-1e6
        eVtonm = 2.*math.pi*197.3269718 # = 2*pi*hbar*c (hbar*c from 2014 PDG)
        xunit = None
        yunit = None
        result = {}
        for line in f:
            if self.debug>1: 
                print line[:80]
                print 'xunit',xunit,'yunit',yunit,'ready',ready,'foundType',foundType
            if line.find(property.upper())>-1 : foundType = True
            if line.find('xunit')>-1 : xunit = self.getUnit(line)
            if line.find('yunit')>-1 : yunit = self.getUnit(line)
            if ready:
                dot = line.find('.')>-1
                if line.find('<!--')>-1 : done = True # comment line, assume this is end of the table
                if not dot and ndata>0:  done = True # end of table found
                if dot and not done:
                    s = line.split()
                    if xunit.lower()=='ev':
                        energy = float(s[0])
                        wl = eVtonm/energy
                        wlmin = min(wlmin,wl)
                        wlmax = max(wlmax,wl)
                    else:
                        wl = float(s[0]) # no coversion
                    if yunit.lower()=='cm':
                        prop = 10.*float(s[1]) # convert cm to mm
                    else:
                        prop = float(s[1]) 
                    result[wl] = prop
                    ndata += 1
            elif done:
                ready = False
            else:
                ready = foundType and line.find('yaxis')>-1
        print 'ProcessMaterialProperty.read',ndata,'lines read. range(',wlmin,',',wlmax,')nm for',property,'& material',material
        return result
    def yint(self,x1,y1,x2,y2,xt):
        dx = x2-x1
        if dx==0. : return y1
        m = (y2-y1)/dx
        b = y2 - m*x2
        return m*xt+b
    def drawRandom(self,material='acrylic', property='abslength'):
        '''
        return random wl for material and property
        '''
        key = material.lower() + '_' + property.lower()
        b = self.sortedIntegralValues[key]
        v = random.uniform(b[0],b[-1])

        j = min(range(len(b)), key=lambda i: abs(b[i]-v))
        a = self.sortedIntegralKeys[key]
        if v<b[j]:
            k = max(0,j-1)
        else:
            k = min(len(b)-1,j+1)
        wlj = a[j]
        wlk = a[k]
        vj  = b[j]
        vk  = b[k]
        return self.yint(vj,wlj, vk,wlk, v)
        
    def getMaterialProperty(self,wl,material='Acrylic',property='ABSLENGTH'):
        '''
        return property at wl for material.
        first check if data has already been read from file, if not, get it.
        then find data points that neighbor input wl and linearly interpolate to find property value.
        if wl is outside range of data, then use closest property value
        '''
        key = material.lower() + '_' + property.lower()
        
        if key in self.Property:
            result = self.Property[key]
            a = self.sortedPropertyKeys[key]
        else:
            result = self.read(material=material,property=property)
            a = result.keys()
            a.sort()
            self.Property[key] = result
            self.sortedPropertyKeys[key] = a

        # original method used next line that traverses entire list
        #j = min(range(len(a)), key=lambda i: abs(a[i]-wl))
        # faster method that assumes monotonically increasing values
        j = None
        dist = 1.e20
        for i in range(len(a)):
            Q = abs(a[i]-wl)
            if Q>dist:
                break
            if Q<dist:
                j = i
                dist = Q
            
        if wl<a[j]:
            k = max(0,j-1)
        else:
            k = min(len(a)-1,j+1)
        wlj = a[j]
        wlk = a[k]
        propj = result[wlj]
        propk = result[wlk]
        return self.yint(wlj,propj, wlk,propk, wl)
    def makeDense(self,wl1=200., wl2=800., dwl=1., property='abslength',material='acrylic'):
        key = material.lower() + '_' + property.lower()
        if self.debug>0: print 'makeDense input key',key
        if key in self.denseProperty:
            if self.debug>0: print 'makeDense already exists for key',key
            return
        
        wl = wl1
        final = {}
        while wl<=wl2:
            final[wl] = self.getMaterialProperty(wl, material=material, property=property)
            wl += dwl
        self.denseProperty[key] = final
        a = final.keys()
        a.sort()
        self.sortedDensePropertyKeys[key] = a
        if self.debug>0: print 'completed makeDense for key',key

        self.makeIntegral(material=material, property=property)
        return
    def makeIntegral(self,material='acrylic', property='abslength'):
        KEY = material.lower() + '_' + property.lower()
        if KEY in self.Integral: return
        a = self.sortedDensePropertyKeys[KEY]
        result = self.denseProperty[KEY]
        final = {}
        rsum = 0.
        for wl in a:
            rsum += result[wl]
            final[wl] = rsum
        self.Integral[KEY] = final
        self.sortedIntegralKeys[KEY] = a

        b = final.values()
        b.sort()
        self.sortedIntegralValues[KEY] = b
        return
    
    def write(self,wl1=200., wl2=800., dwl=1., property='abslength',material='acrylic'):
        key =  material.lower() + '_' + property.lower() 
        fn = key + '.vec'
        
        f = open(fn,'w')

        self.makeDense( wl1=200., wl2=800., dwl=1., property=property,material=material)
        for wl in self.sortedDensePropertyKeys[key]:
            p = self.denseProperty[key][wl]
            f.write( str(wl) + ' ' + str(p) + '\n')
        f.close()
        print 'wrote',len(self.denseProperty[key]),'lines to',fn
        return True
    def iniRandom(self, property='abslength',material='acrylic'):
        wl1 = 200.
        wl2 = 800.
        dwl = 1.
        self.makeDense( wl1=wl1, wl2=wl2, dwl=dwl, property=property, material=material)
        return
    def churn(self):
        wl1 = 200.
        wl2 = 800.
        dwl = 1.
        for m in ['acrylic','water']:
            for p in ['ABSLENGTH', 'RINDEX']:
                self.write(wl1=wl1,wl2=wl2,dwl=dwl,property=p,material=m)
        m = 'scintillator'
        for p in ['FASTCOMPONENT','SLOWCOMPONENT']:
            self.write(wl1=wl1,wl2=wl2,dwl=dwl,property=p,material=m)

        return True
       
    
if __name__ == '__main__':
    print 'I am main'
    pmp = ProcessMaterialProperty(dbg=0)
    if 0: 
        wl = 246
        prop = pmp.getMaterialProperty(wl)
        print wl, prop
        pmp.write()
    if 0: 
        pmp.churn()

    m = 'scintillator'
    for p in ['FASTCOMPONENT', 'SLOWCOMPONENT']:
        pmp.iniRandom(property=p,material=m)
        print 'drawRandom results for',p,m
        for i in range(10):
            wl = pmp.drawRandom(property=p,material=m)
            print wl
    
