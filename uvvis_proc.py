#!/usr/bin/env python
'''
process uv-vis data
20160823
'''
import os
import graphUtils
import ROOT
import math

import numpy
from scipy.interpolate import interp1d

class uvvis_proc():
    def __init__(self):
        self.gU = graphUtils.graphUtils()
        self.Segelstein = None
        self.graphs = []
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
        self.graphs.append(g)
        g.GetXaxis().SetTitle('Wavelength (nm)')
        g.GetYaxis().SetTitle('Absorption (1/cm)')
        self.gU.color(g,1,1,setMarkerType=False)
        self.gU.drawGraph(g,figDir='UVVIS_Results/')
        self.gU.drawGraph(g,figDir='UVVIS_Results/',SetLogy=True)
        self.gU.drawGraph(g,figDir='UVVIS_Results/',SetLogy=True,SetLogx=True)
        return
                                  
    def graphSome(self,dictUVmeas,match=None,figdir='UVVIS_Results/',yAxisLabel='Relative absorbance',SetLogy=False,wlmin=200.,wlmax=800.):
        '''
        return tmultigraph of data in dictUVmeas with key that contains string match
        '''

        tmg = None
        i = 0
        for key in dictUVmeas:
            if match.upper() in key.upper():
                if tmg is None: tmg = self.gU.makeTMultiGraph(match)
                x,y = dictUVmeas[key]
                if wlmin<wlmax:
                    wavelength,absor = dictUVmeas[key]
                    x,y = [],[]
                    for wl,a in zip(wavelength,absor):
                        if wlmin<=wl and wl<=wlmax:
                            x.append(wl)
                            y.append(a)
                title = name = key
                g = self.gU.makeTGraph(x,y,title,name)
                self.graphs.append(g)
                self.gU.color(g,i,i,setMarkerType=False)
                i+=1
                tmg.Add(g)
        print 'uvvis_proc.graphSome yAxisLabel',yAxisLabel,'tmg',tmg       
        self.gU.drawMultiGraph(tmg,figdir=figdir,xAxisLabel='Wavelength(nm)',yAxisLabel=yAxisLabel,SetLogy=SetLogy,abscissaIsTime=False)
        return tmg
    def calcAbsorption(self,wavelength,absorbance,thick=10.):
        '''
        Return absorption as fcn of wavelength given absorbance measurement relative to pure water,
        assuming pure water is described by Segelstein data
        thick = thickness of absorbance cell in cm
        '''
        WL,A = self.readSegelstein()
        x = numpy.array(WL)
        y = numpy.array(A)
        f = interp1d(x,y)
        Absorption = []
        for wl,a in zip(wavelength,absorbance):
            if a>0: 
                b = f(wl) - math.log(a)/thick
            else:
                b = f(wl)
            if b<=0:
                print 'uvvis_proc.calcAbsorption f(wl)',f(wl),'a',a,'b',b
            Absorption.append(b)
        return wavelength,Absorption
            
        
if __name__ == '__main__' :
    UVP = uvvis_proc()
    if 0:
        fn = '../../UVvisData/1Ton_Filled_151217/Baseline_160509.txt'
        h,x,y = UVP.readFile(fn)
        print 'header',h,'len(x)',len(x),'len(y)',len(y),'x[3]',x[3],'y[3]',y[3]
    dn = '../../UVvisData/1Ton_Filled_151217/'
    d = UVP.readAll(dn)
    print 'Parsed',len(d),'data files in',dn
    TMG = []
    if 1:
        for match in ['_IN_','Baseline','_OUT_']:
            TMG.append( UVP.graphSome(d,match=match) )
    x,y = UVP.readSegelstein()
    UVP.drawSegelstein(x,y)

    absorpt = {}
    for key in d:
        wl,a = UVP.calcAbsorption(d[key][0],d[key][1])
        absorpt['A'+key] = [wl,a]
    for match in ['_IN_','Baseline','_OUT_']:
        TMG.append( UVP.graphSome(absorpt,match=match,yAxisLabel='Absorption coeff (1/cm)',SetLogy=False) )

    
    
        
    rfn = 'UVVIS_Results/plots.root'
    rf = ROOT.TFile(rfn,"RECREATE")
    for g in TMG: rf.WriteTObject(g)
    for g in UVP.graphs: rf.WriteTObject(g)
    rf.Close()
    print 'Wrote',len(TMG)+len(UVP.graphs),'multigraphs to',rfn
