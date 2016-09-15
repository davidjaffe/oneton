#!/usr/bin/env python
'''
process uv-vis data
20160823
'''
import os
import graphUtils
import ROOT
import math
import sys

import numpy
from scipy.interpolate import interp1d

class uvvis_proc():
    def __init__(self):
        self.gU = graphUtils.graphUtils()
        self.Segelstein = None
        self.graphs = {}
        self.TMGs   = {}
        self.thick = 10.
        self.ln10 = math.log(self.thick)
        return
    def readFile(self,fn=None,Name=None,useFileName=True):
        '''
        return lists of wavelength and absorbance, also return name
        if useFileName = True, then name = basename without suffix
        if useFileName = False and Name is not None, name is taken from file header
        otherwise name = Name
        '''
        f = open(fn,'r')
        header = None
        wavelength,absorbance = [],[]
        for line in f:
            if header is None:
                header = line.strip('\n').strip('\r').strip('"').replace('-','').replace(' ','_').replace('__','_')
            elif line[0]!='"':
                s = line.split()
                wavelength.append( float(s[0]) )
                absorbance.append( float(s[1]) )
        f.close()
        if useFileName:
            header,ext = os.path.splitext(os.path.basename(fn))
        else:
            if Name is not None: header = Name
        return header,wavelength,absorbance
    def readAll(self,dirpath=None):
        '''
        read all files directory given by dirpath, create dict[name] = { WL,A } where
        WL,A = lists of wavelengths, absorbance
        name = basename of file
        '''
        dictUVmeas = {}
        fileList = []
        if dirpath is not None:
            onlyfiles = [os.path.join(dirpath,f) for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
            for x in onlyfiles:
                prefix,suffix = os.path.splitext(x)
                if suffix=='.txt': fileList.append(x)
        for fn in fileList:
            name,x,y = self.readFile(fn)
            name = name.replace('1T_','OneT_') ### root cannot handle graph names that start with integer
            dictUVmeas[name] = [x,y]
        return dictUVmeas
    def readSegelstein(self,fn='/Users/djaffe/work/UVvisData/segelstein81.dat'):
        '''
        return wavelength(nm) and absorption(1/cm) measured by Segelstein from
        http://omlc.org/spectra/water/data/segelstein81.dat on 20160823
        Also initialize global list variable 
        '''
        if self.Segelstein is None:
            f = open(fn,'r')
            i = 0
            wavelength,absorption = [],[]
            for line in f:
                i+=1
                if i>5:
                    s = line.split()
                    wavelength.append(float(s[0]))
                    absorption.append(float(s[1]))
            f.close()
            self.Segelstein = [wavelength,absorption]
        else:
            wavelength, absorption = self.Segelstein
        return wavelength,absorption
    def drawSegelstein(self,wavelength,absorption,wlmin=100.,wlmax=1000.):
        '''
        plot data from Segelstein
        '''
        name = 'Segelstein_water_transparency'
        title = name.replace('_',' ')
        if wlmin<wlmax:
            title += ' wavelengths=('+str(wlmin)+','+str(wlmax)+') nm'
            wl,a = [],[]
            for x,y in zip(wavelength,absorption):
                if wlmin<=x and x<=wlmax:
                    wl.append(x)
                    a.append(y)
        else:
            wl,a = wavelength,absorption
        g = self.gU.makeTGraph(wl,a,title,name)
        if name in self.graphs:
            sys.exit('uvvis_proc.drawSegelstein ERROR duplicate graph name ' + name)
        else:
            self.graphs[name] = g
        g.GetXaxis().SetTitle('Wavelength (nm)')
        g.GetYaxis().SetTitle('Absorption (1/cm)')
        self.gU.color(g,1,1,setMarkerType=False)
        self.gU.drawGraph(g,figDir='UVVIS_Results/')
        self.gU.drawGraph(g,figDir='UVVIS_Results/',SetLogy=True)
        self.gU.drawGraph(g,figDir='UVVIS_Results/',SetLogy=True,SetLogx=True)
        return
                                  
    def graphSome(self,dictUVmeas,match=None,descrip=None,figdir='UVVIS_Results/',yAxisLabel='Relative absorbance',SetLogy=False,wlmin=200.,wlmax=800.,draw=True):
        '''
        return tmultigraph of data in dictUVmeas with key that contains string match
        '''

        tmg = None
        tmgName = match
        if descrip is not None: tmgName += '_' + descrip
        if tmgName in self.TMGs:
            sys.exit('uvvis_proc.graphSome ERROR Duplicate multigraph name ' + tmgName)
        i = 0
        for key in dictUVmeas:
            if match.upper() in key.upper():
                if tmg is None: tmg = self.gU.makeTMultiGraph(tmgName)
                x,y = dictUVmeas[key]
                if wlmin<wlmax:
                    wavelength,absor = dictUVmeas[key]
                    x,y = [],[]
                    for wl,a in zip(wavelength,absor):
                        if wlmin<=wl and wl<=wlmax:
                            x.append(wl)
                            y.append(a)
                title = name = key
                if name in self.graphs:
                    g = self.graphs[name]
                else:
                    self.graphs[name] = g = self.gU.makeTGraph(x,y,title,name)
                self.gU.color(g,i,i,setMarkerType=False)
                i+=1
                tmg.Add(g)
        #print 'uvvis_proc.graphSome yAxisLabel',yAxisLabel,'tmg',tmg       
        if draw:
            for SetLogy in [False,True]: self.gU.drawMultiGraph(tmg,debugMG=False,figdir=figdir,xAxisLabel='Wavelength(nm)',yAxisLabel=yAxisLabel,SetLogy=SetLogy,abscissaIsTime=False)
        return tmg
    def calcAbsorption(self,wavelength,absorbance,refWL,refAbs,thick=None):
        '''
        Return absorption as fcn of wavelength given absorbance measurement relative to pure water,
        assuming pure water is described by Segelstein data
        thick = thickness of absorbance cell in cm
        '''
        if thick is None : thick = self.thick
        WL,A = self.readSegelstein()
        f    = interp1d( numpy.array(WL),    numpy.array(A) )
        fref = interp1d( numpy.array(refWL), numpy.array(refAbs) )
        Absorption = []
        for wl,a in zip(wavelength,absorbance):
            alpha = a/thick
            refAlpha = fref(wl)/thick
            if refAlpha>=alpha:
                b = f(wl) - (refAlpha - alpha)*self.ln10
            else:
                b = f(wl) 
            Absorption.append(b)
        return wavelength,Absorption
    def groupMeasurements(self,dictUVmeas,method='byDate',avoid=[]):
        '''
        return groups of measurements in dictUVmeas organized by method
        grouping is done using key in dictUVmeas
        at present, only method is byDate
        '''
        debug = False
        groups = {}
        if method=='byDate':
            ## find unique date strings in measurement names
            goodDates = ['16', '2016', '2017']
             
            avoid.extend( [s.lower() for s in avoid] )
            dateStrings = []
            for name in dictUVmeas:
                if name not in avoid:
                    for s in name.split('_'):
                        OK = False
                        for gD in goodDates:
                            if s[0:len(gD)]==gD: OK = True
                        if OK:
                            if s not in dateStrings: dateStrings.append(s)
                            break
            dateStrings.sort()
            if debug : print 'uvvis_proc.groupMeasurements dateStrings',dateStrings
                
            ## now create groups with names that contain the same unique date string
            keyList = dictUVmeas.keys()
            for dS in dateStrings:
                if debug : print 'uvvis_proc.groupMeasurements dS',dS
                used = []
                for key in keyList:
                    if dS in key:
                        if debug : print 'uvvis_proc.groupMeasurements dS',dS,'in key',key
                        if dS not in groups: groups[dS] = []
                        groups[dS].append(key)
                        if debug : print 'uvvis_proc.groupMeasurements groups[dS]',groups[dS]
                        used.append(key)
                if key in used:
                    if debug : print 'uvvis_proc.groupMeasurements remove key',key
                    keyList.remove(key)
            if debug : print 'uvvis_proc.groupMeasurements groups.keys()',groups.keys()
        else:
            sys.exit('uvvis_proc.groupMeasurements ERROR Invalid method=',str(method))
        return groups
    def getRefInGroup(self,groups):
        '''
        return a dict with key = group name, value = list of reference measurements
        '''
        refInGroups = {}
        for gName in groups:
            refs = []
            for firstWord in ['OmegaWater', 'Baseline']:
                for msmt in groups[gName]:
                    if msmt.split('_')[0]==firstWord: refs.append(msmt)
                if len(refs)>0: break
            refInGroups[gName] = refs
        return refInGroups
    
    def plotRaw(self,dictUVmeas,groups,draw=True):
        '''
        make multigraphs from each group of measurements
        '''
        matches = groups.keys()
        for match in matches:
            tmg = self.graphSome(dictUVmeas,match=match,descrip='Raw',draw=draw)
            self.TMGs[tmg.GetName()] = tmg
        return
    def plotAbsorption(self,dictUVmeas,groups,refInGroups,draw=True):
        '''
        plot absorption coefficient (units = 1/cm) from relative absorbance of sample with respect to pure water
        '''
        debug = False
        for gName in groups:
            if debug : print 'uvvis_proc.plotAbsorption gName',gName
            if gName in refInGroups:
                mNames = groups[gName]    # all measurements
                refs = refInGroups[gName] # pure water measurements
                if debug : print 'uvvis_proc.plotAbsorption mNames',mNames
                if debug : print 'uvvis_proc.plotAbsorption refs',refs

                for ref in refs: mNames.remove(ref) # remove pure water msmts from list of all msmts
                if debug : print 'uvvis_proc.plotAbsorption mNames after removing refs',mNames
                for ref in refs:
                    if debug : print 'uvvis_proc.plotAbsorption ref',ref                    
                    refWL,refA = dictUVmeas[ref]
                    absorpt = {}
                    for mName in mNames:
                        mWL,mA = dictUVmeas[mName]
                        wl,a = self.calcAbsorption(mWL,mA,refWL,refA)
                        absorpt[mName] = [wl,a]
                        if debug : print 'uvvis_proc.plotAbsorption mName',mName
                    if debug : print 'uvvis_proc.plotAbsorption absorpt.keys()',absorpt.keys()
                    xdescrip = ''
                    if len(refs)>1: xdescrip = ref[ref.find(gName):].replace(gName,'')
                    tmg = self.graphSome(absorpt,match=gName,descrip='absorpt'+xdescrip,draw=draw,yAxisLabel='Absorption coeff (1/cm)',SetLogy=True)
                    self.TMGs[tmg.GetName()] = tmg
        return
        
if __name__ == '__main__' :
    UVP = uvvis_proc()
    if 0:
        print '\nSingle file test'
        fn = '../../UVvisData/1Ton_Filled_151217/Baseline_160509.txt'
        h,x,y = UVP.readFile(fn)
        print 'header',h,'len(x)',len(x),'len(y)',len(y),'x[3]',x[3],'y[3]',y[3]
        
    dn = '../../UVvisData/1Ton_Filled_151217/'
    d = UVP.readAll(dn)
    print 'Parsed',len(d),'data files in',dn

    groups = UVP.groupMeasurements(d)
    kg = groups.keys()
    kg.sort()
    print 'group names',kg
    rg = UVP.getRefInGroup(groups)
    for name in rg:
        extra = ''
        if len(rg[name])==0: extra = ' '.join([x for x in groups[name]])
        print 'name',name,'rg[name]',rg[name],extra

    
    if 1:
        print '\nProcess raw absorbance measurements'
        UVP.plotRaw(d,groups,draw=True)
                
    x,y = UVP.readSegelstein()
    UVP.drawSegelstein(x,y)

    if 1:
        print '\nTransform to absorption measurements'
        UVP.plotAbsorption(d,groups,rg,draw=True)

        
    rfn = 'UVVIS_Results/plots.root'
    rf = ROOT.TFile(rfn,"RECREATE")
    for g in UVP.TMGs: rf.WriteTObject(UVP.TMGs[g])
    for g in UVP.graphs: rf.WriteTObject(UVP.graphs[g])
    rf.Close()
    print 'Wrote',len(UVP.TMGs),'multigraphs and',len(UVP.graphs),'graphs to',rfn
