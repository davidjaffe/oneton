#!/usr/bin/env python
'''
basic modules to make, print tgraphs
20150921 

'''
import time
import datetime
import sys
import os
import math
import ROOT 
from ROOT import TH1D, TFile, gROOT, TCanvas, TLegend, TGraph, TDatime, TMultiGraph, gStyle, TGraphErrors, TLine, TH2D
from array import array
import re # regular expression
import pipath


class graphUtils():
    def __init__(self,defMarkerSize=1):
        '''
        initialize colors and markers in order preferred by Danielle Berish 20170226
        set default marker size for graphs
        '''
        self.defaultMarkerSize = defMarkerSize
        # use in color()
        self.goodColors = [ROOT.kBlack, ROOT.kBlue, ROOT.kMagenta, ROOT.kGreen, ROOT.kOrange+7, ROOT.kCyan, ROOT.kViolet-1]
        clist = [1,2,3,4, 6,7,8,9, 11,12, 18] # no yellow(5) or white (0,10)
        clist.extend( [x for x in range(28,50)] )
        for c in clist:
            if c not in self.goodColors: self.goodColors.append( c )

        self.goodMarkers = [29, 20, 21, 22, 23, 33, 34]
        for m in range(20,31):
            if m not in self.goodMarkers: self.goodMarkers.append( m )

                
                
        self.pip = pipath.pipath()
        return
    
    def t2dt(self,t):
        '''
        convert struct_time to datetime object
        '''
        return datetime.datetime.fromtimestamp(time.mktime(t))
    def dt2t(self,dt):
        '''
        convert datetime object to struct_time (time object)
        '''
        fmt = "%Y %M %d %H %m %S"
        return time.strptime(dt.strftime(fmt),fmt)
    def addSeconds(self,t,seconds=0):
        '''
        add seconds to struct_time t by converting to datetime object,
        using timedelta and converting back
        '''
        dt = self.t2dt(t)
        dt += datetime.timedelta(seconds=seconds)
        return self.dt2t(dt)
    def convertTime(self,day,fmt,text):
        c = text.count(":")
        if c==1: return time.strptime(day+text,fmt+"%H:%M")
        if c==2: return time.strptime(day+text,fmt+"%H:%M:%S")
        sys.exit("graphUtils.convertTime ERROR Unknown input " + str(text))
        return
    def fixTimeDisplay(self,g,showDate=False,maybeShowDate=True):
        '''
        set time axis to display nicely
        '''
        if g:
            g.GetXaxis().SetTimeDisplay(1)
            g.GetXaxis().SetTimeFormat("%H:%M")
            if showDate:
                g.GetXaxis().SetTimeFormat("#splitline{%H:%M}{%y/%m/%d}")
            else:
                if maybeShowDate:
                    x1 = g.GetXaxis().GetXmin()
                    x2 = g.GetXaxis().GetXmax()
                    if x2-x1>24.*60.*60.:
                        g.GetXaxis().SetTimeFormat("#splitline{%H:%M}{%y/%m/%d}")
                        #print 'graphUtils.fixTimeDisplay: >1 day, so use splitline in SetTimeFormat'
            g.GetXaxis().SetNdivisions(-409)
            g.GetXaxis().SetLabelSize(0.025) #0.5*lx)
            g.GetXaxis().SetTimeOffset(0,"local") # what does this do?
#            g.GetXaxis().SetTimeOffset(0,"gmt") # using gmt option gives times that are only off by 1 hour on tgraph
        else:
            print 'graphUtils.fixTimeDisplay: WARNING Null pointer passed to fixTimeDisplay?????'
        return
    def normHist(self,Hin,divisor,makeNewHist=True,nameSuffix='N',titleSuffix=' per second'):
        '''
        return histogram scaled by 1/divisor
        new histogram is returned, unless makeNewHist=False, then Hin is replaced
        '''
        cn = Hin.ClassName()
        if 'TH1' in cn or 'TH2' in cn:
            if makeNewHist:
                title = Hin.GetTitle()
                name  = Hin.GetName()
                newname = (name + nameSuffix).replace(' ','_')
                hnew = Hin.Clone(newname)
                hnew.SetTitle(title + titleSuffix)
            else:
                hnew = Hin
            hnew.Scale(1./divisor)
        return hnew
    def makeTH2D(self,u,v,title,name,nx=40,xmi=1.,xma=-1.,ny=40,ymi=1.,yma=-1.):
        ''' book, fill 2d hist '''
        if xmi>xma :
            xmi,xma = min(u),max(u)
            dx = (xma-xmi)/float(nx)
            xmi -= dx
            xma += dx
        if ymi>yma :
            ymi,yma = min(v),max(v)
            dy = (yma-ymi)/float(ny)
            ymi -= dy
            yma -= dy
        h = TH2D(name,title,nx,xmi,xma,ny,ymi,yma)
        for x,y in zip(u,v): h.Fill(x,y)
        return h
    def makeTH1D(self,v,title,name,nx=100,xmi=1,xma=-1):
        if xmi>xma:
            xmi = min(v)
            xma = max(v)
            dx = (xma-xmi)/float(nx)
            xmi -= dx/2.
            xma += dx/2.
        h = TH1D(name,title,nx,xmi,xma)
        for y in v: h.Fill(y)
        return h
    def makeTH1Dwtd(self,x,y,title,Name='',NX=None,XMI=None,XMA=None):
        '''
        fill 1d hist with weights y
        given equal size, monotonically increasing bin centers x
        '''
        name = Name
        if Name=='': name = title.replace(' ','_').replace('.','_')
        nx = len(x)
        dx = x[1]-x[0]
        xmi = min(x)-dx/2.
        xma = max(x)+dx/2.
        if NX is not None: nx = NX
        if XMI is not None:xmi =XMI
        if XMA is not None:xma =XMA
        h = TH1D(name,title,nx,xmi,xma)
        for a,b in zip(x,y): h.Fill(a,b)
        ymi,yma = min(y),max(y)
        dy = (yma-ymi)/20.
        ymi,yma = ymi-dy/2.,yma+dy/2.
        h.SetMaximum(yma)
        h.SetMinimum(ymi)
        return h
    def printHistStats(self,h):
        '''
        print some stats for input hist
        '''
        N,mean,stddev,underflow,overflow = self.getHistStats(h)
        print h.GetTitle(),'mean',mean,'stddev',stddev,'Nentries',N,'uflow',underflow,'oflow',overflow
        return
    def getHistStats(self,h):
        '''
        return histogram stats
        '''
        axis = 1 # 1d hist only
        mean = h.GetMean(axis)
        stddev = h.GetStdDev(axis)
        N = h.GetEntries()
        underflow = h.GetBinContent(0)
        if axis==1: nbins = h.GetNbinsX()
        if axis==2: nbins = h.GetNbinsY()
        if axis==3: nbins = h.GetNbinsZ()
        overflow = h.GetBinContent(nbins+1)
        return N,mean,stddev,underflow,overflow
    def drawGraph(self,g,figDir="",SetLogx=False,SetLogy=False,option='APL', verbose=False,abscissaIsTime=False,xLimits=None,yLimits=None,noPopUp=False):
        '''
        output graph to file
        '''
        title = g.GetTitle()
        name  = g.GetName()
        if SetLogx: name += '_logx'
        if SetLogy: name += '_logy'
        if len(figDir) > 0 and figDir[-1] != os.path.sep:
             pdf = self.pip.fix(figDir + '/' + name + '.pdf')
        else:
             pdf   = figDir + name + '.pdf'
        ps = pdf.replace('.pdf','.ps')
        if verbose: print 'drawing Graph:',pdf

    
        xsize,ysize = 1100,850 # landscape style
#        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)

        if abscissaIsTime: self.fixTimeDisplay(g)
        if xLimits is not None:
            g.GetXaxis().SetRangeUser(xLimits[0],xLimits[1])
            g.GetXaxis().SetLimits(xLimits[0],xLimits[1])
        if yLimits is not None:
            g.GetYaxis().SetRangeUser(yLimits[0],yLimits[1])
            g.GetYaxis().SetLimits(yLimits[0],yLimits[1])
        
        g.Draw(option)

        if SetLogy: canvas.SetLogy(1)
        if SetLogx: canvas.SetLogx(1)

        if abscissaIsTime: self.fixTimeDisplay(g)

        if abscissaIsTime: self.fixTimeDisplay(g)
        if xLimits is not None:
            g.GetXaxis().SetRangeUser(xLimits[0],xLimits[1])
            g.GetXaxis().SetLimits(xLimits[0],xLimits[1])
        if yLimits is not None:
            g.GetYaxis().SetRangeUser(yLimits[0],yLimits[1])
            g.GetYaxis().SetLimits(yLimits[0],yLimits[1])

        title = os.path.basename(pdf).replace('.pdf','')
        self.finishDraw(canvas,ps,pdf,ctitle=title)
        return
    def drawFit(self,h,figdir='',SetLogy=False,SetLogx=False,extraName=None):
        '''
        draw histogram with fit parameters
        '''
        name = h.GetName()
        if extraName is not None: name += '_' + extraName
        title = h.GetTitle()
        if SetLogx: name += '_logx'
        if SetLogy: name += '_logy'
        if len(figdir)>0 and figdir[-1]!= os.path.sep:
            pdf = self.pip.fix( figdir + '/' + name + '.pdf')
            ps  = self.pip.fix( figdir + '/' + name + '.ps')
        else:
            pdf = figdir +  name + '.pdf' 
            ps  = figdir +  name + '.ps' 
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)
        
        gStyle.SetOptFit(1111)
        h.Draw()
        if SetLogy: canvas.SetLogy(1)
        if SetLogx: canvas.SetLogx(1)
        
        canvas.Draw()
        canvas.SetGrid(1)
        canvas.SetTicks(1)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        
        canvas.Print(ps,'Landscape')
        os.system('ps2pdf ' + ps + ' ' + pdf)
        if os.path.exists(pdf): os.remove(ps)

        return
    def makeCanvas(self,name='c1',StatSize=None):
        '''
        return standard canvas
        StatSize controls size of text box
        '''
        c1 = ROOT.TCanvas(name)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(1111)
        ROOT.gStyle.SetTitleX(0.8)
        if StatSize is not None:
            ROOT.gStyle.SetStatW(StatSize) # size of stats box (and text?)
            ROOT.gStyle.SetStatH(StatSize)
        c1.SetGrid(1)
        c1.SetTicks(1)
        return c1

        
    def finishDraw(self,canvas,ps,pdf,setGrid=True,setTicks=True,ctitle=None):
        '''
        standard nonsense to finish drawing
        ctitle can be considered 'global' title
        '''
        canvas.Draw()
        canvas.SetGrid(setGrid)
        canvas.SetTicks(setTicks)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        ct = None
        s = None
        if ctitle is not None:
            ct = ROOT.TText(0.5,0.975+0.013,ctitle) # 20170120 move title up a tiny bit
            ct.SetNDC(True) # this is needed to ensure norm'ed coordinates
            ct.SetTextAlign(20) # horizontally centered
            s = ct.GetTextSize()/4.
            ct.SetTextSize(s) # 20170120 make title text smaller to avoid interference with individual hist titles
            ct.Draw()
        if ps is not None:
            canvas.Print(ps,'Landscape')
            os.system('ps2pdf ' + ps + ' ' + pdf)
            if os.path.exists(pdf): os.remove(ps)

        #print 'graphUtils.finishDraw pdf',pdf,'ctitle',ctitle,'text size',s
        canvas.IsA().Destructor(canvas) # avoids seg fault?
        return
    def drawMultiHists(self,histlist,fname='',figdir='',statOpt=1111111,setLogy=False,setLogx=False,dopt='',abscissaIsTime=False,biggerLabels=True,fitOpt=None,Grid=False,forceNX=None,changeColors=False,addLegend=False):
        '''
        draw multiple histograms on single pdf output file
        20161228 histlist can be a list of lists. hists in innermost list are overlaid.
        '''
        nHist = len(histlist)
        if nHist<=0:
            print 'graphUtils.drawMultiHists: ERROR zero length histogram list'
            return
        if nHist==1:
            nX = nY = 1
        else:
            nX = 2
            if forceNX is not None: nX = forceNX
            nY = int(float(nHist)/float(nX) + 0.5)
            nY = int(math.ceil (float(nHist)/float(nX)))
            if nHist<4 and forceNX is None: nX,nY = 1,nHist

        #print 'nHist,nX,nY=',nHist,nX,nY
        # create output directory if it does not exist
        if len(figdir)>0:
            if os.path.isdir(figdir):
               pass
            else:
                try:
                    os.mkdir(figdir)
                except IOError,e:
                    print 'graphUtils.drawMultiHists',e
                else:
                    print 'graphUtils.drawMultiHists created',figdir
        # set output file name and canvas title
        base = figdir        
        ctitle = None
        if fname!='':
            ctitle = fname
        else:
            for h in histlist:
                if type(h) is list:
                    for hh in h:
                        name = hh.GetName()
                        ctitle += name
                else:
                    name = h.GetName()
                    ctitle += name
                if h!=histlist[-1]:
                    ctitle += '_'

        if setLogx:
            ctitle += '_logx'
        if setLogy:
            ctitle += '_logy'
        if len(base)>0 and base[-1] != os.path.sep:
            pdf = self.pip.fix(base + '/' + ctitle + '.pdf')
            ps = self.pip.fix(base + '/' + ctitle + '.ps')
        else:
            pdf = base + ctitle + '.pdf'
            ps = base + ctitle + '.ps'
        
        # open canvas, draw on it
        title = ''
        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)

        gStyle.SetOptStat(statOpt)
        if fitOpt is not None: gStyle.SetOptFit(fitOpt)
        spaceBtwnPads = 0.01 / 1000.
        canvas.Divide(nX,nY,spaceBtwnPads,spaceBtwnPads)
        for i,h in enumerate(histlist):
            canvas.cd(i+1).SetLogy(setLogy)
            canvas.cd(i+1).SetLogx(setLogx)
            canvas.cd(i+1).SetGridx(Grid)
            canvas.cd(i+1).SetGridy(Grid)


            if type(h) is list: # overlay
                prefix = ''
                
                if addLegend: lg = self.createLegend(h)
                for icol,hh in enumerate(h):
                    if abscissaIsTime : self.fixTimeDisplay(hh)
                    if changeColors: hh.SetLineColor(self.goodColors[icol])
                    hh.Draw(dopt+prefix)
                    if 'func' in dopt.lower() : hh.Draw('funcsame')
                    prefix = 'same'
                    self.biggerLabels(hh, biggerLabels)
                    if abscissaIsTime:self.fixTimeDisplay(hh)
                if addLegend: lg.Draw()
            else:
                if abscissaIsTime : self.fixTimeDisplay(h)
                h.Draw(dopt)
                if 'func' in dopt.lower() : h.Draw('funcsame')
                self.biggerLabels(h, biggerLabels)
                if abscissaIsTime : self.fixTimeDisplay(h)
            #print i+1,h.GetName()

        self.finishDraw(canvas,ps,pdf,ctitle=ctitle)
        return
    def drawMultiObjects(self,objlist,fname='',figdir='',statOpt=1111111,setLogy=False,setLogx=False,
                         dopt='',gopt='AP',abscissaIsTime=False,biggerLabels=True,fitOpt=None,Grid=False,
                         forceNX=None,changeColors=False,addLegend=False,debug=False,noPopUp=True):
        '''
        draw multiple objects (histograms, multiple hists, graphs, multigraphs) on single pdf output file
        20161228 histlist can be a list of lists. hists in innermost list are overlaid.
        '''

        
        nObj = len(objlist)
        if nObj<=0:
            print 'graphUtils.drawMultiObjects: ERROR zero length input list'
            return
        if nObj==1:
            nX = nY = 1
        else:
            nX = 2
            if forceNX is not None: nX = forceNX
            nY = int(float(nObj)/float(nX) + 0.5)
            nY = int(math.ceil (float(nObj)/float(nX)))
            if nObj<4 and forceNX is None: nX,nY = 1,nObj

        if debug : print 'graphUtils.drawMultiObjects nObj,nX,nY=',nObj,nX,nY
        # create output directory if it does not exist
        if len(figdir)>0:
            if os.path.isdir(figdir):
               pass
            else:
                try:
                    os.mkdir(figdir)
                except IOError,e:
                    print 'graphUtils.drawMultiObjects:',e
                else:
                    print 'graphUtils.drawMultiObjects: created',figdir
        # set output file name and canvas title
        base = figdir        
        ctitle = None
        if fname!='':
            ctitle = fname
        else:
            for h in objlist:
                if type(h) is list:
                    for hh in h:
                        name = hh.GetName()
                        ctitle += name
                else:
                    name = h.GetName()
                    ctitle += name
                if h!=objlist[-1]:
                    ctitle += '_'

        if setLogx:
            ctitle += '_logx'
        if setLogy:
            ctitle += '_logy'
        if len(base)>0 and base[-1] != os.path.sep:
            pdf = self.pip.fix(base + '/' + ctitle + '.pdf')
            ps = self.pip.fix(base + '/' + ctitle + '.ps')
        else:
            pdf = base + ctitle + '.pdf'
            ps = base + ctitle + '.ps'
        
        # open canvas, draw on it
        title = ''
        xsize,ysize = 1100,850 # landscape style
        #noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)

        gStyle.SetOptStat(statOpt)
        if fitOpt is not None: gStyle.SetOptFit(fitOpt)
        spaceBtwnPads = 0.01 / 1000.
        canvas.Divide(nX,nY,spaceBtwnPads,spaceBtwnPads)
        for i,h in enumerate(objlist):
            if type(setLogy) is list:
                canvas.cd(i+1).SetLogy(setLogy[i])
            else:
                canvas.cd(i+1).SetLogy(setLogy)
            canvas.cd(i+1).SetLogx(setLogx)
            canvas.cd(i+1).SetGridx(Grid)
            canvas.cd(i+1).SetGridy(Grid)

            if 'func' in dopt.lower() :
                ROOT.gStyle.SetOptFit(1111)
            else:
                ROOT.gStyle.SetOptFit(0)

            if type(h) is list: # overlay
                prefix = ''
                gopt1 = gopt
                
                if addLegend: lg = self.createLegend(h)
                for icol,hh in enumerate(h):
                    if abscissaIsTime : self.fixTimeDisplay(hh)
                    if changeColors: hh.SetLineColor(self.goodColors[icol])
                    if debug : print 'graphUtils.drawMultiObjects i+1,hh.ClassName()',i+1,hh.ClassName()
                    self.Draw(hh,dopt=dopt+prefix,gopt=gopt1,showFit='func' in dopt.lower(),debug=debug)
                    gopt1 = gopt1.replace('A','').replace('a','')
                    prefix = 'same'
                    self.biggerLabels(hh, biggerLabels)
                    if abscissaIsTime:self.fixTimeDisplay(hh)
                if addLegend: lg.Draw()
            else:
                if abscissaIsTime : self.fixTimeDisplay(h)
                if debug : print 'graphUtils.drawMultiObjects i+1,h.ClassName()',i+1,h.ClassName()
                self.Draw(h,dopt=dopt,gopt=gopt,showFit='func' in dopt.lower(),debug=debug)
                self.biggerLabels(h, biggerLabels)
                if abscissaIsTime : self.fixTimeDisplay(h)


        self.finishDraw(canvas,ps,pdf,ctitle=ctitle)
        return
    def objWasFit(self,obj):
        '''
        return True if input object was fitted
        '''
        wasFit = False
        for a in obj.GetListOfFunctions():
            if 'TF1' in a.ClassName():
                print 'graphUtils.objWasFit obj.GetName()',obj.GetName(),'a,a.ClassName(),a.GetName()',a,a.ClassName(),a.GetName()
                wasFit = True
        return wasFit
    def Draw(self,obj,dopt='',gopt='APL',showFit=True,debug=False):
        '''
        Draw input object using options dopt or gopt
        if showFit, add fit results
        '''
        cName = obj.ClassName()
        if debug : print 'graphUtils.Draw obj,cName,dopt,gopt,showFit',obj,cName,dopt,gopt,showFit
        wasFit = showFit and self.objWasFit(obj)
        if wasFit: ROOT.gStyle.SetOptFit(1111)
        if 'TH1' in cName or 'TH2' in cName:
            obj.Draw(dopt)
            if wasFit: obj.Draw("func same")
        elif 'TGraph' in cName or 'TMultiGraph'==cName:
            obj.Draw(gopt)
        else:
            print 'graphUtils.Draw WARNING Incapable of drawing',cName
        return
    def biggerLabels(self,h,Factor=2.0):
        '''
        increase axis label size

        default factor of 2.0 was determined empirically
        Check on type of factor for backward compat
        '''
        factor = Factor
        if type(Factor) is bool:
            if not Factor : return 
            factor = 2.0 # empirically determined
        sx = h.GetXaxis().GetLabelSize()
        h.GetXaxis().SetLabelSize(factor*sx)
        sy = h.GetYaxis().GetLabelSize()
        h.GetYaxis().SetLabelSize(factor*sy)

        return
        
    def drawMultiGraph(self,TMG,figdir='',SetLogy=False, SetLogx=False, abscissaIsTime = True, drawLines=True, xAxisLabel=None,yAxisLabel=None,NLegendColumns=None, maxTitleLength=9999999, debugMG = False):
        '''
        draw TMultiGraph with legend and output as pdf
        Default is that abscissa is calendar time.
        Returns canvas
        
        '''
        if debugMG: print 'graphUtils.drawMultiGraph:TMG',TMG,'TMG.GetListOfGraphs()',TMG.GetListOfGraphs(),\
            'SetLogy',SetLogy,'SetLogx',SetLogx,'drawLines',drawLines,'xAxisLabel',xAxisLabel,'yAxisLabel',yAxisLabel,'NLegendColumns',NLegendColumns
       
        if not TMG.GetListOfGraphs(): return  # empty
        title = TMG.GetTitle()
        name  = TMG.GetName()
        if SetLogx: name += '_logx'
        if SetLogy: name += '_logy'
        if debugMG: print 'graphUtils.drawMultiGraph',title,name,'TMG.GetListOfGraphs()',TMG.GetListOfGraphs(),'TMG.GetListOfGraphs().GetSize()',TMG.GetListOfGraphs().GetSize()
        nGraphs = TMG.GetListOfGraphs().GetSize()

        if len(title)>maxTitleLength and 'splitline' not in title:
            for i in range(maxTitleLength,len(title)):
                if title[i]==' ': break
            t1 = title[:i]
            t2 = title[i:]
            title = '#splitline{'+t1+'}{'+t2+'}'
            TMG.SetTitle(title)

        if len(figdir)>0 and figdir[-1] != os.path.sep:
            pdf = self.pip.fix(figdir + '/' + name + '.pdf')
            ps = self.pip.fix(figdir + '/' + name + '.ps')
        else:
            pdf = figdir + name  + '.pdf'
            ps  = figdir + name  + '.ps'
        #print 'ps=',ps
        #print 'pip.fix',self.pip.fix(figdir + '/' + name + '.ps')

        xsize,ysize = 1100,850 # landscape style
        noPopUp = True
        if noPopUp : gROOT.ProcessLine("gROOT->SetBatch()")
        canvas = TCanvas(pdf,title,xsize,ysize)
        canvas.SetLogy(SetLogy)
        canvas.SetLogx(SetLogx)

        # move title to left in order to put legend above plot
        gStyle.SetTitleX(0.3+0.025)
        x1 = 0.5
        x2 = x1 + .5
        y1 = 0.9
        y2 = y1 + .1
        lg = TLegend(x1,y1,x2,y2)
        NGraph = 0
        for g in TMG.GetListOfGraphs():
            NGraph += 1
            t = g.GetTitle()
            lg.AddEntry(g, t, "LP" )
            if abscissaIsTime : self.fixTimeDisplay(g)
        if NGraph>6: lg.SetNColumns(2)
        if NLegendColumns is not None: lg.SetNColumns(NLegendColumns)

        dOption = "AP"
        if drawLines: dOption += "L"

        # complicated monkey business because of idiotic way that logY is set
        if SetLogy:
            ymi,yma = 1.e20,1.e-20
            for g in TMG.GetListOfGraphs():
                x,y = self.getPoints(g)
                ymi = min(ymi,min(y))
                yma = max(yma,max(y))
            if ymi<=0: ymi = 0.1
            ymi = ymi/2.
            yma = 2.*yma
            TMG.SetMinimum(ymi)
            TMG.SetMaximum(yma)
            for g in TMG.GetListOfGraphs():
                g.SetMinimum(ymi)
                g.SetMaximum(yma)
                if "A" in dOption:
                    g.SetTitle( TMG.GetTitle() )
                    if xAxisLabel is not None: g.GetXaxis().SetTitle(xAxisLabel)
                    if yAxisLabel is not None: g.GetYaxis().SetTitle(yAxisLabel)
                g.Draw(dOption)
                if debugMG: print 'graphUtils.drawMultiGraph: Draw',g.GetName()
                dOption = dOption.replace("A","")
        else:
            TMG.Draw(dOption)
            if debugMG: print 'graphUtils.drawMultiGraph: Draw',TMG.GetName()
            if xAxisLabel is not None: TMG.GetXaxis().SetTitle(xAxisLabel)
            if yAxisLabel is not None: TMG.GetYaxis().SetTitle(yAxisLabel)
            if abscissaIsTime : self.fixTimeDisplay(TMG)


        self.labelTMultiGraph(TMG,xAxisLabel=xAxisLabel,yAxisLabel=yAxisLabel,debug=debugMG)
        lg.Draw()
        canvas.Draw()
        canvas.SetGrid(1)
        canvas.SetTicks(1)
        canvas.cd()
        canvas.Modified()
        canvas.Update()
        if 0:
            canvas.Print(pdf,'pdf')
        else:
            self.canvasPrint(canvas,ps)
#            canvas.Print(ps,'Landscape')
#            os.system('ps2pdf ' + ps + ' ' + pdf)
#            if os.path.exists(pdf): os.remove(ps)

        if debugMG: print 'graphUtils.drawMultiGraph',title,'complete'
        return canvas
    def canvasPrint(self,canvas,ps):
        canvas.Print(ps,'Landscape')
        pdf = ps.replace('.ps','.pdf')
        os.system('ps2pdf ' + ps + ' ' + pdf)
        if os.path.exists(pdf): os.remove(ps)
        return 
    def createLegend(self,objlist,x1=0.5,y1=0.8,popt="L"):
        
        '''
        return TLegend object from input list of objects
        x1,y1 = 0.5,0.9 places legend in upper left corner of pad
        '''
        x2 = x1 + .5
        y2 = y1 + .1
        lg = TLegend(x1,y1,x2,y2)
        for obj in objlist:
            t = obj.GetTitle()
            lg.AddEntry(obj, t, popt )
        return lg
    def makeTMultiGraph(self,name,tit=None,debug=False):
        title = tit
        if tit is None:title = name.replace('_',' ')
        tmg = TMultiGraph()
        tmg.SetName(name)
        tmg.SetTitle(title)
        if debug:
            print 'graphUtils.makeTMultiGraph:name',name,'title',title,'object',tmg
        return tmg
    def labelTMultiGraph(self,tmg,xAxisLabel=None,yAxisLabel=None,debug=False):
        name = tmg.GetName()
        if not tmg.GetXaxis() or tmg.GetYaxis() :
            if debug: 'graphUtils.labelTMultiGraph name',name,'not drawn so cannot label axes'
            return
        if debug : print 'graphUtils.labelTMultiGraph name',name,'tmg',tmg,'tmg.GetXaxis()',tmg.GetXaxis(),'tmg.GetYaxis()',tmg.GetYaxis()
        if 'vs' in name:
            s = name.split('_')
            xt = s[2]
            yt = s[0]
            xt = xt.replace('by','/')
            xt = xt.replace('BY','/')
            yt = yt.replace('by','/')
            yt = yt.replace('BY','/')
            if debug:
                print 'graphUtils.labelTMultiGraph: xt,yt',xt,yt,'tmg',tmg
                print 'tmg.GetXaxis()',tmg.GetXaxis(),'tmg.GetYaxis()',tmg.GetYaxis()
            if tmg.GetXaxis(): tmg.GetXaxis().SetTitle(xt)
            if tmg.GetYaxis(): tmg.GetYaxis().SetTitle(yt)
        if xAxisLabel is not None:
            tmg.GetXaxis().SetTitle(xAxisLabel)
            if debug: print 'graphUtils.labelTMultiGraph:xAxisLabel',xAxisLabel
        if yAxisLabel is not None:
            tmg.GetYaxis().SetTitle(yAxisLabel)
            if debug: print 'graphUtils.labelTMultiGraph:yAxisLabel',yAxisLabel
        return
    def makeTGraph(self,u,v,title,name,ex=None,ey=None,markerSize=None):
        '''
        make TGraph with axes u,v. ex,ey are uncertainties in u,v if not None
        '''
        if ex is None:
            g = TGraph(len(u),array('d',u), array('d',v))
        else:
            dy = ey
            if ey is None: dy = [0. for x in range(len(ex))]
            g = TGraphErrors(len(u),array('d',u),array('d',v),array('d',ex),array('d',dy))
        g.SetTitle(title)
        g.SetName(name)
        m = markerSize
        if m is None: m = self.defaultMarkerSize
        g.SetMarkerSize(m)
        return g
    def rebinByWeek(self,t,y,dt,dy,byDay=False):
        '''
        given timestamp, ordinate, timestamp uncertainty, ordinate uncertainty
        return same averaged by week or by day if byDay==True.
        uncertainty in timestamp is spread,
        ordinate uncertainty is weighted average uncertainty
        
        datetime.date.isocalendar() is an instance-method returning a tuple
        containing year, weeknumber and weekday in respective order for the given date instance.

        timetuple is a time.struct_time object with attributes
        tm_yday with range [1,366]
        tm_year i.e. 2016
        '''
        ts,week,WeekNs = [],[],[]
        for x in t:
            dt = datetime.datetime.fromtimestamp(x)
            ts.append( dt ) # convert to datetime stamp
            if byDay:
                w = dt.timetuple().tm_year*1000 + dt.timetuple().tm_yday
            else:
                w  = dt.isocalendar()[1] 
            week.append( w )
            if w not in WeekNs: WeekNs.append( w )
                
        T,Y,DT,DY = [],[],[],[]
        WeekNs.sort()
        for wknum in WeekNs:
            u,v,dv = [],[],[]
            for i,wk in enumerate(week):
                if wk==wknum:
                    u.append( t[i] )
                    v.append( y[i] )
                    dv.append( dy[i] )
            T.append( sum(u)/float(len(u)) )
            DT.append( 0.5*(max(u)-min(u)) )
            wt = []
            for x in dv: wt.append( 1./x/x )
            Y.append( sum([ a*b for a,b in zip(v,wt)])/sum(wt) )
            DY.append( 1./math.sqrt(sum(wt)))
        return T,Y,DT,DY
    def color(self,obj,n,M,setMarkerColor=False,setMarkerType=True):
        '''
        set line color and marker type for obj based on indices n and M
        if M=n then use M to set marker type, otherwise determine marker type from n
        unless setMarkerType is False
        '''
        debug = False
        LC = len(self.goodColors)
        LM = len(self.goodMarkers)
        c = n%LC
        obj.SetLineColor( self.goodColors[c] )
        if debug: print 'color: obj',obj,'n',n,'obj.IsA().GetName()',obj.IsA().GetName()
        if setMarkerType:
            oName = obj.IsA().GetName()
            if oName=='TGraph' or oName=='TGraphErrors':
                if M==n:
                    m = M%LM
                else:
                    m = int(float(n)/float(LC))%LM
                obj.SetMarkerStyle( self.goodMarkers[m] )
                if setMarkerColor: obj.SetMarkerColor( self.goodColors[c] )
                if debug: print 'color:',obj.GetName(),'m',m,'self.goodMarkers[m]',self.goodMarkers[m]
        return
    def getPoints(self,g,getErrors=False):
        '''
        return abscissa,ordinate values of input graph g
        also return errors if getErrors is True
        '''
        x,y = [],[]
        if getErrors: dx,dy = [],[]
        for i in range(g.GetN()):
            a,b = ROOT.Double(0),ROOT.Double(0)
            OK = g.GetPoint(i,a,b)
            if OK!=-1:
                x.append(a)
                y.append(b)
                if getErrors:
                    dx.append(g.GetErrorX(i))
                    dy.append(g.GetErrorY(i))
        if getErrors: return x,y,dx,dy
        return x,y
    def getTDatime(self,dt,fmt='%Y/%m/%d %H:%M:%S'):
        '''
        convert date/time text to TDatime object
        '''
        datetimeObj = self.getdatetime(dt,fmt=fmt)
        return TDatime( datetimeObj.strftime('%Y-%m-%d %H:%M:%S') ).Convert()
    def getdatetime(self,dt,fmt='%Y/%m/%d %H:%M:%S'):
        ''' convert timestamp dt to text '''
        return datetime.datetime.strptime(dt,fmt)
    def reportHist(self,h):
        '''
        write out some properties of hist h
        '''
        name = h.GetName()
        title = h.GetTitle()
        xa = h.GetXaxis()
        nx = xa.GetNbins()
        xmi= xa.GetXmin()
        xma= xa.GetXmax()
        xex= xa.CanExtend()
        nd = h.GetDimension()
        words = 'graphUtils.reportHist',name,title,'nx,xmi,xma',nx,xmi,xma
        if xex: words += 'can extend x-axis.'
        if nd>1:
            ya = h.GetYaxis()
            ny = ya.GetNbins()
            ymi= ya.GetXmin()
            yma= ya.GetXmax()
            yex= ya.CanExtend()
            words += 'ny,ymi,yma=',ny,ymi,yma
            if yex: words += 'can extend y-axis.'
        print words
        return
