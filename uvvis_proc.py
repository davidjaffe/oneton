#!/usr/bin/env python
'''
process uv-vis data
20160823
'''
import os
import graphUtils
import ROOT

class uvvis_proc():
    def __init__(self):
        self.gU = graphUtils.graphUtils()
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
            dictUVmeas[name] = [x,y]
        return dictUVmeas
    def readSegelstein(self,fn='/Users/djaffe/work/UVvisData/segelstein81.dat'):
        '''
        return wavelength(nm) and absorption(1/cm) measured by Segelstein from
        http://omlc.org/spectra/water/data/segelstein81.dat on 20160823
        '''
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
        g.GetXaxis().SetTitle('Wavelength (nm)')
        g.GetYaxis().SetTitle('Absorption (1/cm)')
        self.gU.color(g,1,1,setMarkerType=False)
        self.gU.drawGraph(g,figDir='UVVIS_Results/')
        self.gU.drawGraph(g,figDir='UVVIS_Results/',SetLogy=True)
        self.gU.drawGraph(g,figDir='UVVIS_Results/',SetLogy=True,SetLogx=True)
        return
                                  
    def graphSome(self,dictUVmeas,match=None,figdir='UVVIS_Results/'):
        '''
        return tmultigraph of data in dictUVmeas with key that contains string match
        '''
        graphs = []
        tmg = None
        i = 0
        for key in dictUVmeas:
            if match.upper() in key.upper():
                if tmg is None: tmg = self.gU.makeTMultiGraph(match)
                x,y = dictUVmeas[key]
                title = name = key
                g = self.gU.makeTGraph(x,y,title,name)
                graphs.append(g)
                self.gU.color(g,i,i,setMarkerType=False)
                i+=1
                tmg.Add(g)
                
        self.gU.drawMultiGraph(tmg,figdir=figdir,xAxisLabel='Wavelength(nm)',yAxisLabel='Relative absorbance',abscissaIsTime=False)
        return tmg
        
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
    for match in ['_IN_','Baseline','_OUT_']:
        TMG.append( UVP.graphSome(d,match=match) )
    x,y = UVP.readSegelstein()
    UVP.drawSegelstein(x,y)
    
        
    rfn = 'UVVIS_Results/plots.root'
    rf = ROOT.TFile(rfn,"RECREATE")
    for g in TMG: rf.WriteTObject(g)
    rf.Close()
    print 'Wrote',len(TMG),'multigraphs to',rfn
