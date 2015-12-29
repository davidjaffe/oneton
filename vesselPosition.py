#!/usr/bin/env python
'''
convert measured positions of vessel to positions in x,y
20151218
'''

import csv
import graphUtils
from ROOT import TFile
import math

gU = graphUtils.graphUtils()

## read measurement data from spread sheet
fn = 'Measured_positions_book3p35_20151218.csv'
#print 'list dialects:',csv.list_dialects()
columns = {}
colhdrs = []
with open(fn,'rU') as csvfile:
    rdr = csv.reader(csvfile,delimiter=',')
    #print rdr
    first = True
    for l in rdr:
        #print len(l),l
        if first:
            first = False
            for x in l:
                if x!='':
                    colhdrs.append(x)
                    columns[x] = []
        else:
            for i,x in enumerate(l):
                if i<len(colhdrs):
                    if x=='':
                        fx = 0.
                    else:
                        fx = float(x)
                        if i==0: fx = int(x)
#                    print 'i',i,'x',x,'fx',fx
                    columns[colhdrs[i]].append(fx)
#for h in colhdrs:print h,columns[h]


## additional measurements to better constrain radial position of vessel bottom
## U1,V1 and U7,V7 are positions wrt platform U(Y) is in X(Y) direction
U1 = [15.,18.,21.,24.,27.,30.,33.,36.] # inches
U1 = [x*2.54 for x in U1] # convert to cm
V1 = [20.6,16.7,14.2,12.8,12.8,14.0,16.1,19.9] # cm, estimated uncertainty is +-0.3cm

U7 = [15.,18.,21.,24.,27.,30.,33.,36.] # inches
U7 = [x*2.54 for x in U7] # convert to cm
V7 = [20.0,16.1,13.6,12.3,12.3,13.2,15.4,19.0] # cm

# only 6 measurements in spreadsheet should be used.
lastPosition = max(columns['position'])-1
lastIndex    = columns['position'].index(lastPosition)
#print 'lastPostion',lastPosition,'lastIndex',lastIndex
lastIndex += 1
#print 'new lastIndex',lastIndex,columns['position'][:lastIndex]

## calculate translation of U,V measurements to X,Y
## requires extrapolation of platform x-positions (called r,l in spreadsheet) to table x-positions
offset = columns['kcm'][0]
halfbar = 1.5/2.*2.54
print 'offset',offset,'halfbar',halfbar
sprime = [offset/2.-x-halfbar for x in columns['scm'][:lastIndex]]
mprime = [x+halfbar-offset/2. for x in columns['mcm'][:lastIndex]]
rprime = columns['rcm'][:lastIndex]
lprime = columns['lcm'][:lastIndex]
print 'sprime',sprime
print 'rprime',rprime
print 'mprime',mprime
print 'lprime',lprime

Slope,Intercept = {},{}
for Side in ['right','left']:
    R,S = rprime,sprime
    if Side=='left': R,S = lprime,mprime

    Srr,Sr,Ss,Srs = 0.,0.,0.,0.
    for r,s in zip(R,S):
        Sr += r
        Srr+= r*r
        Ss += s
        Srs+= r*s
    N = float(len(rprime))
    slope = (Sr*Sr - N*Srr)/(Sr*Ss - N*Srs)
    intercept = (Sr - slope*Srs)/N
    x0 = slope*offset/2.+intercept
    x7 = -slope*offset/2.+intercept
    print Side,'side: slope',slope,'intercept',intercept,'x0',x0,'x7',x7
    Slope[Side],Intercept[Side] = slope,intercept

    
# use right side to correct U1,V1, U7,V7
Xq,Yq = [],[]
for u1,v1 in zip(U1,V1):
    delta = Slope['right']*(-offset/2.) + Intercept['right']
    x = u1 - delta - offset/2.
    y = offset/2. - v1
    Xq.append(x)
    Yq.append(y)
for u7,v7 in zip(U7,V7):
    delta = Slope['right']*(offset/2.) + Intercept['right']
    x = u7 - delta - offset/2.
    y = -offset/2. + v7
    Xq.append(x)
    Yq.append(y)

#    put information into more convenient form
X,Y = [],[]
for w in ['Lcm','Rcm']:
    c = 'X'+w
    #print 'c',c,columns[c][:lastIndex]
    X.extend(columns[c][:lastIndex])
    c = 'Y'+w
    Y.extend(columns[c][:lastIndex])

# append additional measurements        
#X.extend(Xq)
#Y.extend(Yq)

# estimated measurement uncertainties
sY = [1./16.*2.54 for y in Y]
sX = []
for i,x in enumerate(X):
    if i%6==0 or i%6==5:
        sX.append(1./4.*2.54)
    else:
        sX.append(1./16.*2.54)


#print 'X',X,'Y',Y
#print 'len(X),len(Y)',len(X),len(Y)
# guess at center and radius of circle        
xave = sum(X)/float(len(X))
yave = sum(Y)/float(len(Y))
rave = 0.
for x,y in zip(X,Y):
    rave += (x-xave)*(x-xave) + (y-yave)*(y-yave)
rave = rave/float(len(X))
rave = math.sqrt(rave)
print 'Initial: xave,yave',xave,yave,'rave',rave

### find best R,x0,y0 using chi2 raster scan
chimin = 1.e20
bestRxy = []

dr = 1.
nr = 11
dx = 0.1
nx = 11
dy = 0.1
ny = 11
nloop = 3
debug = False
if debug: nloop = 1
for loops in range(nloop):
    if len(bestRxy)>0:
        dr = dr/float(nr-1)
        dx = dx/float(nx-1)
        dy = dy/float(ny-1)
        rave,xave,yave = bestRxy
    for ix in range(nx):
        x0 = xave + float(ix-nx/2)*dx
        for iy in range(ny):
            y0 = yave + float(iy-ny/2)*dy
            for i in range(nr):
                Rguess = rave + float(i-nr/2)*dr
                chi2 = 0.
                for j,x in enumerate(X):
                    y = Y[j]
                    sy=sY[j]
                    sx=sX[j]
                    r = math.sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) )
                    den = sx*sx + sy*sy
                    chi2 += (r-Rguess)*(r-Rguess)/den
                if debug: print 'i',i,'ix',ix,'iy',iy,'R,x0,y0',Rguess,x0,y0,'chi2',chi2
                if chi2<chimin:
                    chimin = chi2
                    bestRxy = [Rguess,x0,y0]
    print 'loops',loops,'best chi2',chimin,'R,x0,y0',bestRxy

## list contribution to chi2
Rguess,x0,y0 = bestRxy
chi2 = 0.
Phi,Residual = [],[]
for j,x in enumerate(X):
    y = Y[j]
    sy=sY[j]
    sx=sX[j]
    r = math.sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) )
    phi = math.atan2(y-y0,x-x0) + math.pi
    den = sx*sx + sy*sy
    residual = (r-Rguess)*(r-Rguess)/den
    chi2 += residual
    Phi.append(phi)
    Residual.append(residual)
    print 'x,y,phi,residual',x,y,phi,residual
print 'chi2',chi2

rf = TFile('vP.root','RECREATE')
title = 'residual vs phi'
name = title.replace(' ','_')
g = gU.makeTGraph(Phi,Residual,title,name)
gU.color(g,0,0,setMarkerColor=True)
gU.drawGraph(g,option='AP')
rf.WriteTObject(g)

tmg = gU.makeTMultiGraph('Fit_circle_and_measured_points') 

phi0 = -math.pi
nphi = 1000
dphi = 2.*math.pi/float(nphi)
Xfit,Yfit = [],[]
for i in range(nphi):
    phi = phi0 + float(i)*dphi
    x = x0 + Rguess*math.cos(phi)
    y = y0 + Rguess*math.sin(phi)
    Xfit.append(x)
    Yfit.append(y)
name = 'fitted_circle'
title = name.replace('_',' ')
g = gU.makeTGraph(Xfit,Yfit,title,name)
gU.color(g,2,2,setMarkerColor=False,setMarkerType=False)
tmg.Add(g)
rf.WriteTObject(g)

name = 'measured_points'
title = name.replace('_',' ')
g = gU.makeTGraph(X,Y,title,name)
gU.color(g,0,0,setMarkerColor=True)
g.SetLineColor(10) # white lines
tmg.Add(g)
gU.drawMultiGraph(tmg,abscissaIsTime=False,xAxisLabel='X(cm) in table coord system',yAxisLabel='Y(cm) in table coord system')
rf.WriteTObject(g)
rf.WriteTObject(tmg)

rf.Close()
