#!/usr/bin/env python
'''
process uv-vis data
20160823
'''
import os
import graphUtils

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
    def plotSome(self,dictUVmeas,match=None):
        '''
        graph data in dictUVmeas with key that contains string match
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
                
        self.gU.drawMultiGraph(tmg,xAxisLabel='Wavelength(nm)',yAxisLabel='Relative absorbance',abscissaIsTime=False)
        return
        
if __name__ == '__main__' :
    UVP = uvvis_proc()
    fn = '../../UVvisData/1Ton_Filled_151217/Baseline_160509.txt'
    h,x,y = UVP.readFile(fn)
    print 'header',h,'len(x)',len(x),'len(y)',len(y),'x[3]',x[3],'y[3]',y[3]
    d = UVP.readAll('../../UVvisData/1Ton_Filled_151217/')
    print 'len(d)',len(d)
    UVP.plotSome(d,match='_IN_')
    UVP.plotSome(d,match='Baseline')
    UVP.plotSome(d,match='_OUT_')
