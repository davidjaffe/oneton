#!/usr/bin/env python
'''
process sphere and box fluorescence data from FelixGX
20160901
'''
import os,sys
import graphUtils
import ROOT
import math

import numpy
#from scipy.interpolate import interp1d

class felix():
    def __init__(self):
        self.gU = graphUtils.graphUtils()
        self.graphs = {}
        self.graphSuffixes = ['raw','corr','exInt','cByr']
        self.TMG = {}
        self.icol = 1

        self.figdir = 'Felix_Results/'

        self.fluors = ['bisMSB','POPOP','PPO']
        self.solvents = ['ETOH','CX','cyclo','LAB']

        return
    def clean(self,line):
        return line.strip('\n').strip('\r')
    def loopy(self):
        '''
        given top directory path, create list of all 'good' files, then process data in the files
        at present, only emission scan data is deemed 'good'
        '''
        printAllFileNames = False

        
        #filenames = [ '/Users/djaffe/work/GIT/QY/Henry/Sphere/bisMSB_LAB/EmissionScan_bisMSBinLAB_4.47mgL_ex310_2sec_160829.txt' ]
        filenames = []
        #dirpath = '/Users/djaffe/work/GIT/QY/Henry/Sphere/bisMSB_LAB/'
        dirpath = '/Users/djaffe/work/GIT/QY/Henry/'
        onlyfiles = self.get_filepaths(dirpath)
        for fn in onlyfiles:
            if self.goodFile(fn) : filenames.append(fn)
        #print onlyfiles
        results = {}
        for fn in filenames:
            if printAllFileNames: print fn
            key,value = self.getEmissionScan(fn=fn,debug=False)
            if value is not None:
                if key in results:
                    sys.exit('felix.loopy ERROR key ' + str(key) + ' is already in results')
                results[key] = value
        #self.plotEmissionScan(results)
        self.analyzeLikeResults(results)

        
        return
    def analyzeLikeResults(self,results):
        '''
        plot groups of similar results
        i.e., same solvent with same excitation wavelength
        '''
        debugThis = False
        
        like = self.collect(results)
        print 'like.keys()',like.keys()
        offset = 'EmissionScan'
        for key in ['Sphere_bisMSB','Sphere_PPO']:
            for base in ['LAB_ex','ETOH_ex']:
                for msmtName in like[key]:
                    if base in msmtName:
                        exWL1 = results[msmtName][4]
                        baseNames = [msmtName]
                        for name in like[key]:
                            if name!=msmtName:
                                exWL2 = results[name][4]
                                if exWL1==exWL2:
                                    baseNames.append(name)
                        if len(baseNames)>1:
                            if debugThis : print '\nfelix.analyzeLikeResults',msmtName
                            sDict,fDict = {},{} # dictionaries of like msmts for solvent,fluor. key=name solvent or fluor
                            for name in baseNames:
                                sf,isSolvent,isFluor = self.solventOrFluor(name,offset=offset)
                                if debugThis : print name,sf,isSolvent,isFluor
                                if isSolvent and isFluor:
                                    sys.exit('felix.analyzeLikeResults: ERROR Both solvent AND fluor for '+name)
                                if sf is not None:
                                    if isSolvent:
                                        if sf not in sDict: sDict[sf] = []
                                        sDict[sf].append(name)
                                    if isFluor:
                                        if sf not in fDict: fDict[sf] = []
                                        fDict[sf].append(name)
                                    if not isSolvent and not isFluor:
                                        sys.exit('felix.analyzeLikeResults: ERROR Neither solvent OR fluor for sf '+sf+' name '+name)
                                else:
                                    print 'felix.analyzeLikeResults: FAIL',name,'NOT SOLVENT OR FLUOR'
                            if debugThis:
                                for k in sDict: print k,' '.join([w for w in sDict[k]])
                                for k in fDict: print k,' '.join([w for w in fDict[k]])
                            for s in sDict:
                                for solvent in sDict[s]:
                                    for f in fDict:
                                        for fluor in fDict[f]:
                                            QYresults = self.getQY(solvent,fluor,results)
            

        return
    def getQY(self,solvent,scint,results,plot=True):
        '''
        interface to QY calculation
        inputs: solvent name=key, scint name=key, dict of results
                boolean plot
        output: QY, uncertainties
        '''
        debug = True
        
        ## stupidity check
        if solvent not in results or scint not in results:
            print 'felix.getQY: ERROR solvent',solvent,'or scint',scint,'is not a valid key for results'
            return None

        sWL, sRaw, sCorr, sexInt, sexWL, sdatetime = results[solvent]
        fWL, fRaw, fCorr, fexInt, fexWL, fdatetime = results[scint]

        ## additional stupidity checks
        if sexWL!=fexWL:
            print 'felix.getQY: ERROR mis-match in excitation WL. solvent',sexWL,'scint',fexWL
            return None
        goodWL = len(sWL)==len(fWL)
        if goodWL: goodWL = (sWL==fWL).all()
        if not goodWL:
            print 'felix.getQY: ERROR mis-match in length or content of WL for solvent and scint. Lengths',len(sWL),len(fWL)
            return None
        
        rawQY,L  = self.analQY1( sWL, sexWL, sRaw, fRaw  ,debug=debug)
        meanWLraw,sRawBL,fRawBL,Araw,Eraw = L
        if plot:
            self.icol = 1
            tmg = self.gU.makeTMultiGraph(solvent)
            tmg.Add( self.makeGraph(sWL, sRaw,   solvent, 'raw') )
            tmg.Add( self.makeGraph(sWL, sRawBL, solvent, 'rawBL') )
            tmg.Add( self.makeGraph(sWL, fRaw,   scint,   'raw') )
            tmg.Add( self.makeGraph(sWL, fRawBL, scint,   'rawBL') )
            tmg.Add( self.makeGraph(sWL, Araw,   scint,   'Abs') )
            tmg.Add( self.makeGraph(sWL, Eraw,   scint,   'Emiss') )
            self.TMG[solvent] = tmg
            self.gU.drawMultiGraph(tmg, figdir=self.figdir, abscissaIsTime=False, xAxisLabel = 'Wavelength (nm)')

        
        corrQY,L = self.analQY1( sWL, sexWL, sCorr,fCorr ,debug=debug)
        meanWLcorr,sCorrBL,fCorrBL,Acorr,Ecorr = L
        QYresults = [rawQY, corrQY]
        if debug:
            print '\nfelix.getQY:',solvent,scint
            print '      rawQY {0:.3f} corrQY {1:.3f}'.format(rawQY,corrQY)
        return QYresults
    def makeGraph(self,X,Y,Name,Suffix):
        '''
        return TGraph with name generated from Name and Suffix and line, point color from global icol
        add TGraph to global dict
        '''
        icol = self.icol
        title = Name.replace('_',' ') + ' ' + Suffix
        name  = Name.replace('.','x') + '_' + Suffix
        g = self.gU.makeTGraph(list(X),list(Y),title,name)
        self.graphs[name] = g
        g.GetXaxis().SetTitle('Wavelength (nm)')
        self.gU.color(g,icol,icol,setMarkerType=False)
        self.icol += 1
        return g
    def solventOrFluor(self,inputName,offset=None):
        '''
        return name of solvent or fluor in solution based on input name, booleans for solvent,flour
        First check for fluor in solution, then solvent
        '''
        name = inputName
        solvent,fluor = False,False
        if offset is not None:
            if offset in inputName: name = inputName[inputName.find(offset):]
        solvent,fluor = False,True
        for f in self.fluors:
            if f in name: return f,solvent,fluor
            if f.lower() in name.lower(): return f,solvent,fluor
        solvent,fluor = True,False
        for s in self.solvents:
            if s in name: return s,solvent,fluor
            if s.lower() in name.lower(): return s,solvent,fluor
        return None,False,False
    def collect(self,results,debug=False):
        '''
        try to collect similar measurements
        where similar means the first 2 words in name of measurement are identical
        '''
        like = {}
        used = []
        keylist = []
        for key in results.keys():
            if results[key] is not None: keylist.append(key)
        for i,k1 in enumerate(keylist):
            s1 = k1.split('_')
            for k2 in keylist[i+1:]:
                if k2 not in used:
                    s2 = k2.split('_')
                    if s1[0]==s2[0] and s1[1]==s2[1]:
                        lkey = s1[0] + '_' + s1[1]
                        if lkey in like:
                            like[lkey].append(k2)
                        else:
                            like[lkey] = [k1,k2]
                            used.append(k1)
                        used.append(k2)
        print 'felix.collect Found',len(like),'like measurements -------------- '
        if debug:
            for key in like:
                print key
                for v in like[key]: print v
                print ''
        return like

    def get_filepaths(self,directory):
        '''
        20160906 taken from http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
        This function will generate the file names in a directory 
        tree by walking the tree either top-down or bottom-up. For each 
        directory in the tree rooted at directory top (including top itself), 
        it yields a 3-tuple (dirpath, dirnames, filenames).
        '''
        file_paths = []  # List which will store all of the full filepaths.

        # Walk the tree.
        for root, directories, files in os.walk(directory):
            for filename in files:
                # Join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                file_paths.append(filepath)  # Add it to the list.

        return file_paths  # Self-explanatory.

    def goodFile(self,fn):
        '''
        return True if filename fn is likely to be an interesting emission scan
        '''
        avoidThese = ['BoxContam','1sec','DONOTUSE','NelsonLabNotebook']
        for words in avoidThese:
            if words in fn: return False
        if '.txt' in fn: return True
        return False
    def getMsmtName(self,fn):
        '''
        return measurement name based on full path name of file 
        '''
        sd = os.path.split(fn)
        name = sd[1].strip('.txt')
        for subdir in reversed(sd[0].split('/')):
            if subdir=='Henry':
                break
            name = subdir + '_' + name
        return name
    def getEmissionScan(self,fn=None,debug=False):
        '''
        return name and emission scan for given input file
        returns None if no Emission scan available in file
        '''
        if fn is None:
            sys.exit('felix.getEmissionScan ERROR No input file')
        else:
            msmtName = self.getMsmtName(fn) #os.path.basename(fn).strip('.txt')
            msmt = self.readEmissionScan(fn,debug=debug)
            return msmtName,msmt
    def analQY1(self,WL,exWL,solvent,scint,debug=False):
        '''
        attempt to compute quantum yield given sphere measurement spectra of solvent only and scintillator
        WL are the wavelengths in nm where measurements were made
        exWL is the excitation wavelength
        '''
        peakR = 5.0
        peakSideBand = 2.
        peakB = peakR + peakSideBand
        emissionWidth = 200.

        minEm = exWL + peakB
        maxEm = minEm + emissionWidth


        # analysis check. How close is measured excitation wavelength to specified wavelength
        meanWL = []
        for sample in [solvent,scint]:
            q,n = self.integrateBetween(WL,sample,exWL-peakR,exWL+peakR,getMeanX=True)
            meanWL.append( q )
        if debug: print 'felix.analQY1: mean excitation WL',meanWL,'excitation WL',exWL

        # absorption peak analysis: get baseline for each sample, subtract baselines, then absorption = solvent - scint
        alo=exWL-peakB
        ahi=alo+peakSideBand
        blo=exWL+peakR
        bhi=blo+peakSideBand
        Xsb,Ysb = self.getBaseline(WL,solvent,alo=alo,ahi=ahi,blo=blo,bhi=bhi,Method=1)
        Xfb,Yfb = self.getBaseline(WL,scint  ,alo=alo,ahi=ahi,blo=blo,bhi=bhi,Method=1)
        sA = numpy.subtract(solvent,Ysb)
        fA = numpy.subtract(scint  ,Yfb)
        A = numpy.subtract(sA,fA)
        absorption,n = self.integrateBetween(WL,A,exWL-peakR,exWL+peakR)
        
        # emission peak analysis: use solvent spectrum as baseline for scint
        E = numpy.subtract(scint,solvent) # for baseline subtraction of emission spectrum
        emission,n = self.integrateBetween(WL,E,minEm,maxEm)
        
        QY = -1.
        if absorption>0.: QY = emission / absorption

        if debug: print 'felix.analQY1 emission',emission,'absorption',absorption,'QY','{0:.3f}'.format(QY)
        
        return QY,[meanWL,Ysb,Yfb,A,E]
    def getBaseline(self,X,Y,alo=None,ahi=None,blo=None,bhi=None,Method=0):
        '''
        return numpy arrays U,V giving estimated baseline of Y from lower sideband (alo,ahi) and/or upper sideband (blo,bhi)
        based on Method
        Method = 0 : simple average of sidebands
        Method = 1 : linear interpolation of sidebands iff two sidebands given, otherwise revert to Method=0
        Method = other => use Method=0
        '''
        a,b = None,None
        if alo is not None and ahi is not None:
            a,na = self.integrateBetween(X,Y,alo,ahi)
            if Method==1: amean,na = self.integrateBetween(X,Y,alo,ahi,getMeanX=True)
        if blo is not None and bhi is not None:
            b,nb = self.integrateBetween(X,Y,blo,bhi)
            if Method==1: bmean,nb = self.integrateBetween(X,Y,blo,bhi,getMeanX=True)
        if a is None and b is None:
            sys.exit('felix.getBaseline: ERROR lower and upper sidebands improperly defined')

        method = Method
        if a is None:
            a,na = b,nb
            method = 0
        if b is None:
            b,nb = a,na
            method = 0

        y1,y2 = a/float(na),b/float(nb)
        if method==1:
            x1,x2 = amean,bmean
            m = (y2-y1)/(x2-x1)
            b = y1 - m*x1
        else:
            m = 0.
            b = (y1+y2)/2.
        U = numpy.array(X)
        V = numpy.add(numpy.multiply(X,m),b)
        return U,V
    
    def integrateBetween(self,x,y,xlo,xhi,getMeanX=False):
        '''
        if getMeanX is False: return integral of y for xlo<=x and x<=xhi and number of bins 
        if getMeanX is True : return mean x weighted by y for xlo<=x and x<=xhi and number of bins
        '''
        s,n = 0.,0
        if getMeanX: sx = 0.
        for u,v in zip(x,y):
            if u>xhi: break
            if xlo<=u:
                s += v
                n += 1
                if getMeanX: sx += u*v
        if getMeanX:
            if s!=0.: s = sx/s
                
        return s,n
    
    def plotEmissionScan(self,results):
        i = 0
        for msmtName in results:
            WL, raw, corr, exInt, exWL, mdatetime = results[msmtName]
            cByr = self.ratio(raw,corr)
            x = list(WL)
            #print x
            for Y,suf in zip([raw,corr,exInt,cByr],self.graphSuffixes):
                y = list(Y)
                #print y
                name = str(msmtName + '_' + suf).replace('.','x')
                title = name.replace('_',' ')
                g = self.gU.makeTGraph(x,y,title,name)
                self.graphs[name] = g
                g.GetXaxis().SetTitle('Wavelength (nm)')
                i += 1
                self.gU.color(g,i,i,setMarkerType=False)
                

        return
    def ratio(self,raw,corr):
        '''
        make new array that is ratio of corr to raw
        '''
        return numpy.divide(corr,raw)
    def readEmissionScan(self,fn=None,debug=False):
        '''
        return emission scan results as list of arrays: WL,raw,corr,exInt; float: exWL; list of strings: mdatetime
        WL, raw   =         (wavelength, raw data)
        WL, corr  =         (wavelength, corrected data)
        WL, exInt =        (wavelength, excitation intensity)
        exWL      =         excitation wavelength 
        mdatetime =         measurement date & time
        '''
        warningMsgs = False
        sessionTitle, Group = self.readFile(fn,debug=debug)
        mdatetime = sessionTitle.split()[-2:]
        gtitle = 'Detector1'
        exWL, emWL, A = Group[gtitle]
        if exWL is None:
            if warningMsgs: print 'felix.readEmissionScan',fn,'Not an emission scan'
            return None
        if A.shape[1]!=4:
            if warningMsgs: print 'felix.readEmissionScan FAIL',fn,'A.shape[1]',A.shape[1],'!=4 as expected for Emission Scan'
            return None
        WL   = A[:,0]
        raw  = A[:,1]
        corr = A[:,3]
        gtitle = 'RCQCSignal'
        nothing, notmuch, A = Group[gtitle]
        exInt = A[:,1]
        return [WL, raw, corr, exInt, exWL, mdatetime]
    def readFile(self,fn=None,debug=False):
        '''
        read and unpack a single file
        returns
        Session information (mainly date and time)
        `Group` information (measurement data, including excitation or emission wavelength)
        Measurement data as (2*Np,Ncol) arrays
        Excitation or emission wavelength as float
        Group `title'
        Results should be further unpacking in calling routine
        
        File format:
        <Session>
        Title (incudes date, time)
        <Group>
        Detector1
        # of pairs of columns = Np
        # of entries per column (has Np entries) = Nc
        Description of Np columns
        'X Y' repeated Np times
        Np pairs of columns
         :  :
        </Group>
        <Group>
         : :
        </Group>
        </Session>
        '''
        
        f = open(fn,'r')
        lines = f.readlines()
        f.close()
        if debug: print 'felix.readFile opened',f
                
        sessionTitle = None
        Group = {}
        i = 0
        while i<len(lines):
            line = lines[i]
            if '<Session>' in line:
                i += 1
                sessionTitle = self.clean(lines[i])
                if debug: print 'sessionTitle',sessionTitle
            elif '<Group>' in line:
                i += 1
                gtitle = self.clean(lines[i])
                if debug: print 'gtitle',gtitle
                i += 1
                Npairs = int(self.clean(lines[i]))
                if debug: print 'Npairs',Npairs
                i += 1
                s = self.clean(lines[i]).split()
                Ncols = []
                for x in s: Ncols.append(int(x))
                if debug: print 'Ncols',Ncols
                i += 1
                excitationWavelength, emissionWavelength = None,None
                header = self.clean(lines[i]).split()
                if header[0]=='D1':
                    s = header[1].split(':')
                    if '-' not in s[1]:
                        emissionWavelength = float(s[1])
                    if '-' not in s[0]:
                        excitationWavelength = float(s[0])
                    if debug: print 'excitationWavelength',excitationWavelength,'emissionWavelength',emissionWavelength
                i += 1
                columnLabels = self.clean(lines[i]).split()
                A = None
                for j in range(max(Ncols)):
                    i += 1
                    cols = []
                    for x in self.clean(lines[i]).split(): cols.append( float(x) )
                    if A is None:
                        A = numpy.array(cols)
                    else:
                        A = numpy.vstack( (A, numpy.array(cols)) )
                Group[gtitle] = [excitationWavelength, emissionWavelength, A]
                if debug: 
                    D = {}
                    for d in range(2*Npairs):
                        D[d] = A[:,d]
                    for d in range(2*Npairs):
                        print 'd',d,'D[d][0]',D[d][0],'D[d][-1]',D[d][-1]
            i += 1
        return sessionTitle, Group
if __name__ == '__main__' :
    F = felix()
    if 0:
        print F.getMsmtName('/Users/djaffe/work/GIT/QY/Henry/SyncScans/lampChange/empyIS_72.txt')
        sys.exit(1)
    F.loopy()
    rfn = 'Felix_Results/plots.root'
    rf = ROOT.TFile(rfn,"RECREATE")

    for g in F.graphs: rf.WriteTObject(F.graphs[g])
    for g in F.TMG: rf.WriteTObject(F.TMG[g])
    rf.Close()
    print 'Wrote',len(F.graphs)+len(F.TMG),'objects to',rfn

        
