#!/usr/bin/env python
'''
plot position of barrier as measured
azimuthal positions denoted by bolt hole number from 1 to 16
measured depth of top of barrier from top of flange and
radial displacement of top of barrier wrt acrylic vessel
vessel internal height
vessel inner diameter 
measurement units were inches
20151211
'''

import graphUtils
from ROOT import TFile,TMultiGraph

gU = graphUtils.graphUtils()
Graphs = {}

hDesign = 1224.6 # mm, design inner height
IDDesign= [994.70, 996.72] # mm, measured + design inner diameter

depth = [1./2., 11./16., 3./4., 3./4., 5./8., 11./16., 21./32., 17./32., 1./2., 5./8., 21./32., 5./8., 19./32., 1./2., 7./16., 7./16.]
raddis= [3./16., 1./6., 1./2., 5./8., 3./8., 17./32., 15./16., 1.+(15./32.), 15./16., 3./16., 1./8.,3./16., 1./8.,15./32., 11./16., 1./8.]
bolt = [a+1 for a in range(16)]

# inner diameter. key = bolt, value = id. duplicate values for bolt and bolt+8
ID = {1: 39.+5./16., 7: 39.+7./16., 5: 39.+5./16., 3: 39.+5./16.}
newID = {}
for b in ID:
    newID[b+8] = ID[b]
    newID[b] = ID[b]
ID = newID

# height. key=bolt, value =height}
# use average of multiple measurements at bolt#7
height = {1: 49.25, 3:49.+5./16., 5:49.+5./16., 7: 49.+1./3.*(7./16.+3./8.+5./16.), 9: 49.+5./16., 11: 49.+5./16., 13: 49.+5./16., 15:49.+5./16.}

# interpolate
newheight = {}
newID = {}
for b in bolt:
    if b in height:
        newheight[b] = height[b]
    else:
        blo = b - 1
        bhi = b + 1
        if bhi>max(sorted(height)): bhi = min(sorted(height))
        newheight[b] = 0.5*(height[blo]+height[bhi])
    if b in ID:
        newID[b] = ID[b]
    else:
        blo = b - 1
        bhi = b + 1
        if bhi>max(sorted(height)): bhi = min(sorted(height))
        newID[b] = 0.5*(ID[blo]+ID[bhi])
        
height = newheight
ID = newID

# barrier height
barH = {}
for b,d in zip(bolt,depth):
    barH[b] = (height[b]-d)


# barrier inner diameter
barID = {}
for b,dr in zip(bolt,raddis):
    barID[b] = ID[b] - dr

#### graphs
mgh = gU.makeTMultiGraph('heights')
mgd = gU.makeTMultiGraph('diameters')
mgm = gU.makeTMultiGraph('barrier_msmts')

# design values
X = [1,16]
Y = [hDesign, hDesign]
title = 'vessel design inner height mm'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgh.Add(g)

X = [1,16]
for i,ind in enumerate(IDDesign):
    Y = [ind,ind]
    title = 'vessel design IDlo mm'
    if i==1: title = 'vessel design IDhi mm'
    name = title.replace(' ','_')
    Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
    mgd.Add(g)

# barrier height vs bolt
X = bolt
Y = []
for b in sorted(barH):
    Y.append(barH[b]*25.4)
title = 'barrier height mm vs bolt'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgh.Add(g)

Y = []
for b in sorted(barID):
    Y.append(barID[b]*25.4)
title = 'barrier ID mm vs bolt'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgd.Add(g)    

X = bolt
Y = [a*25.4 for a in depth]
title = 'depth mm vs bolt'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgm.Add(g)

Y = [a*25.4 for a in raddis]
title = 'radial displacement mm vs bolt'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgm.Add(g)

X,Y = [],[]
for b in sorted(height):
    X.append(b)
    Y.append(height[b]*25.4)
title = 'vessel height mm vs bolt'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgh.Add(g)

X,Y = [],[]
for b in sorted(ID):
    X.append(b)
    Y.append(ID[b]*25.4)
title = 'vessel ID mm vs bolt'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X,Y,title,name)
mgd.Add(g)

# set line color and point types
icol = 0
ipt  = len(Graphs)
for g in Graphs:
    icol += 1
    gU.color(Graphs[g],icol,ipt,setMarkerColor=True)
    ipt -= 1

# for writing out the multigraphs
MGs = [mgd,mgm,mgh]
for mg in MGs:
    gU.drawMultiGraph(mg,abscissaIsTime=False,figdir='Figures/')


        
ofn = 'barrier.root'
outf= TFile.Open(ofn,'RECREATE')
for g in Graphs: outf.WriteTObject(Graphs[g])
for g in MGs: outf.WriteTObject(g)
outf.Close()
print 'Wrote',len(Graphs)+len(MGs),'graphs to',ofn
