#!/usr/bin/env python
'''
convert measured positions of vessel to positions in x,y
20151218
20160201 more work
20160204 add bottom PMT positions
'''

import csv
import graphUtils
from ROOT import TFile
import math
import sys

gU = graphUtils.graphUtils()

### basic stuff
offsetINCH = 51.  # width,length of platform in inches
barINCH    = 1.5  # width of 80-20 bar in inches
offset = offsetINCH * 2.54 # width and length of platform
bar    = barINCH*2.54   # width of 80-20 bar
halfbar = bar/2.
fpt   = 0.25*2.54 # 80-20 fastener plate thickness in cm
Rmushield = 6.0/2.# mu-metal shield radius from R7723 ass'y drawing
Rpmt = 5.2/2. # nominal radius of PMT window from R7723 ass'y drawing (diameter=52 +-1.5mm)
dQpmt = 25. # spacing set between bottom PMTs

Rdesign = 104.58/2. # design OD/2

print 'width of platform',offset,'bar',bar,'design Radius',Rdesign
print 'fastener plate thickness',fpt,'radius of R7723 PMT mu-metal shield',Rmushield,'radius of R7723 PMT window',Rpmt,'spacing btwn bottom PMTs',dQpmt

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
print 'Data read from',fn
for h in colhdrs:print h,columns[h]
print ' '
# validity check of conversion from inches to cm
for h in colhdrs:
    if len(h)==1:
        hcm = h+'cm'
        inches = columns[h]
        cm     = columns[hcm]
        ok = True
        for i,c in zip(inches,cm):
            if i==c and i==0.:
                pass 
            elif abs(c/i-2.54)>1.e-4:
                ok = False
                j = inches.index(i)
                print 'ERROR',h,hcm,'inches',i,'cm',c,'ratio-2.54',c/i-2.54
        if ok: print h,hcm,'valid conversion of inches to cm'


# only 6 measurements in spreadsheet should be used.
lastPosition = max(columns['position'])-1
lastIndex    = columns['position'].index(lastPosition)
#print 'lastPostion',lastPosition,'lastIndex',lastIndex
lastIndex += 1
#print 'new lastIndex',lastIndex,columns['position'][:lastIndex]

# 20160201 add new columns q, k for right,left horiz positions of vertical 80/20 bars of base 
#                          n, t for right,left vertical positions of vertical 80/20 bars of base 

colhdrs.append('q')
columns['q'] = [offsetINCH - barINCH - r for r in columns['r'][:lastIndex]]
#colhdrs.append['k']  # overwriting, so don't append
columns['k'] = [l + barINCH for l in columns['l'][:lastIndex]]
colhdrs.append('t')
columns['t'] = [offsetINCH - s + barINCH/2. for s in columns['s'][:lastIndex]]
colhdrs.append('n')
columns['n'] = [offsetINCH - m + barINCH/2. for m in columns['m'][:lastIndex]]

for h in ['q','k','t','n']:
    hcm = h+'cm'
    columns[hcm] = [x*2.54 for x in columns[h]]
    if hcm not in colhdrs: colhdrs.append(hcm)

# horizontal center positions
h = 'horizcenter'
colhdrs.append(h)
columns[h] = [0.5*(k+q) for k,q in zip(columns['kcm'],columns['qcm'])]
h = 'backfrontcenter'
colhdrs.append(h)
a = 0.5*(columns['tcm'][0] + columns['tcm'][lastIndex-1])
b = 0.5*(columns['ncm'][0] + columns['ncm'][lastIndex-1])
columns[h] = [a,b]

## X,Y positions of right,left side measurement of bottom of vessel
hR,hL = 'newXRcm','newXLcm'
for h in [hR,hL]:
    colhdrs.append(h)
    columns[h] = []
for i in range(lastIndex):
    columns[hR].append(  0.5*(offset + columns['rcm'][i] - columns['lcm'][i]) - columns['Rcm'][i] )
    columns[hL].append( columns['Lcm'][i] - 0.5*(offset + columns['lcm'][i] - columns['rcm'][i]) )
hR,hL = 'newYRcm','newYLcm'
for h in [hR,hL]:
    colhdrs.append(h)
    columns[h] = []
for i in range(lastIndex):
    columns[hR].append( columns['tcm'][i] - columns['backfrontcenter'][0] )
    columns[hL].append( columns['ncm'][i] - columns['backfrontcenter'][1] )

h = 'newChordcm'
colhdrs.append(h)
columns[h] = [XR-XL for XR,XL in zip(columns['newXRcm'],columns['newXLcm'])]

h = 'horizBarLencm'
colhdrs.append(h)
columns[h] = [q-k for q,k in zip(columns['qcm'],columns['kcm'])]


print 'Data read from',fn,'WITH ADDITIONAL DATA AND OVERWRITING'
for h in colhdrs:print h,columns[h]
print ' '

# compare with design
h = 'horizBarLencm'
design = 109.538
lo,hi = min(columns[h])-design,max(columns[h])-design
print 'Range of',h,'wrt to design',design,'is',lo,hi
    



## additional measurements to better constrain radial position of vessel bottom
## U1,V1 and U7,V7 are positions wrt platform U(Y) is in X(Y) direction
U1 = [15.,18.,21.,24.,27.,30.,33.,36.] # inches
sU1 = [1./4.*2.54 for x in U1] # 1/4" uncertainty
U1 = [x*2.54 for x in U1] # convert to cm
V1 = [20.6,16.7,14.2,12.8,12.8,14.0,16.1,19.9] # cm, estimated uncertainty is +-0.3cm
sV1 = [0.3 for x in V1]

U7 = [15.,18.,21.,24.,27.,30.,33.,36.] # inches
sU7 = [1./4.*2.54 for x in U7]
U7 = [x*2.54 for x in U7] # convert to cm
V7 = [20.0,16.1,13.6,12.3,12.3,13.2,15.4,19.0] # cm
sV7 = [0.3 for x in U7]

## 20160204 new translations of U,V measurements to X,Y
Xuv,Yuv = [],[]
dXback = 0.5*(offset + columns['rcm'][0] - columns['lcm'][0])
dY = 0.5*(columns['tcm'][0] + columns['tcm'][lastIndex-1])
for u,v in zip(U1,V1):
    Xuv.append( dXback - u )
    Yuv.append( offset - v - dY )
sXuv,sYuv = [],[]
sXuv.extend( sU1 )
sYuv.extend( sV1 )
dXfront = 0.5*(offset + columns['rcm'][lastIndex-1] - columns['lcm'][lastIndex-1])
for u,v in zip(U7,V7):
    Xuv.append( dXfront - u )
    Yuv.append( v - dY )
sXuv.extend( sU7 )
sYuv.extend( sV7 )    

useOldUV = False
if useOldUV:

    ## calculate translation of U,V measurements to X,Y
    ## requires extrapolation of platform x-positions (called r,l in spreadsheet) to table x-positions
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
    Xq,Yq,sXq,sYq = [],[],[],[]
    for u1,v1 in zip(U1,V1):
        delta = Slope['right']*(-offset/2.) + Intercept['right']
        x = u1 - delta - offset/2.
        y = offset/2. - v1
        Xq.append(x)
        Yq.append(y)
    sXq.extend( sU1 )
    sYq.extend( sV1 )
    for u7,v7 in zip(U7,V7):
        delta = Slope['right']*(offset/2.) + Intercept['right']
        x = u7 - delta - offset/2.
        y = -offset/2. + v7
        Xq.append(x)
        Yq.append(y)
    sXq.extend( sU7 )
    sYq.extend( sV7 )

    Xqave = sum(Xq)/float(len(Xq))
    Yqave = sum(Yq)/float(len(Yq))

# put information into more convenient form
X,Y = [],[]
newX,newY = [],[]
for w in ['Lcm','Rcm']:
    c = 'X'+w
    #print 'c',c,columns[c][:lastIndex]
    X.extend(columns[c][:lastIndex])
    c = 'Y'+w
    Y.extend(columns[c][:lastIndex])
    c = 'newX'+w
    newX.extend(columns[c][:lastIndex])
    c = 'newY'+w
    newY.extend(columns[c][:lastIndex])

#### replace old with new for fit
X,Y = newX,newY
    
# estimated measurement uncertainties
sY = [1./16.*2.54 for y in Y]
sX = []
for i,x in enumerate(X):
    if i%6==0 or i%6==5:
        sX.append(1./4.*2.54)
    else:
        sX.append(1./16.*2.54)

# append additional measurements and uncertainties
useAddlMsmts = True
useNewUV = True
if useAddlMsmts:
    if useNewUV:
        X.extend(Xuv)
        Y.extend(Yuv)
        sX.extend(sXuv)
        sY.extend(sYuv)
    else:
        if useOldUV:
            X.extend(Xq)
            Y.extend(Yq)
            sX.extend(sXq)
            sY.extend(sYq)


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
    #print 'x,y,phi,residual',x,y,phi,residual
print 'chi2',chi2,'for',len(X),'points and',len(bestRxy),'parameters'

calcHorizDisp = False
if calcHorizDisp:
### calculate horizontal displacement of XR,XL derived measurements wrt 'nominal' vessel positions
    print '\n horiz displacement wrt (0,0) (dx0) and wrt uv measurements (dxq)'
    C,D = [],[]
    for x,y in zip(newX,newY):
        #print 'x,y,Rdesign,y0,y-y0',x,y,Rdesign,y0,y-y0
        s = math.sqrt(Rdesign*Rdesign - (y-0.)*(y-0.))
        xp = 0.+s
        xm = 0.-s
        dx0 = x-xm
        if abs(x-xp)<abs(x-xm):
            dx0 = x-xp
        C.append(y)
        D.append(abs(dx0))
        if useOldUV:
            s = math.sqrt(Rdesign*Rdesign - (y-Yqave)*(y-Yqave))
            xp = Xqave+s
            xm = Xqave-s
            dxq = x-xm
            if abs(x-xp)<abs(x-xm):
                dxq = x-xp
        print 'y,x',y,x,'dx0',dx0,'dxq',dxq
    print ''

### positions of bottom PMTs 20160204
### Assume PMTs in contact with bottom of vessel and bottom of vessel defines Z=0
dY = bar/2. + fpt + bar + fpt + bar + Rmushield
y = columns['tcm'][2] - 0.5*(columns['tcm'][0]+columns['tcm'][lastIndex-1]) 
Ypmt = [y - dY for i in range(4)] # all bottom PMTs at same Y position
Xpmt = []
x = columns['rcm'][2] - columns['lcm'][2] + 1.5*dQpmt
Xpmt.append(x)
for i in range(3):
    x -= dQpmt
    Xpmt.append(x)
Zpmt = [0. for i in range(len(Xpmt))]

### top PMTs
Xleft = -(offset - columns['rcm'][2] + columns['lcm'][2])/2.
x4 = x5 = Xleft + 41.6 + bar + fpt + bar/2.
Xpmt.extend([x4,x5])
Yback = offset - 0.5*(columns['ncm'][0] + columns['ncm'][lastIndex-1])
y4 = Yback - (26.+9./16.)*2.54 + bar/2. + 35.8 - bar - Rmushield
y5 = Yback - (26.+9./16.)*2.54 + bar/2  - 18.7 + bar + Rmushield
Ypmt.extend([y4,y5])
Ztop = 130.4925
Zpmt.extend([Ztop,Ztop])

### positions of bottom hodoscopes H4, H5 (give X,Y,Z of corners starting at back right)
### H5
Xright = 0.5*(offset + columns['rcm'][0] - columns['lcm'][0])
xbr = Xright - (8.+9./16.)*2.54
xbl = Xright - (20.5)*2.54
xfr = Xright - (8.5)*2.54
xfl = Xright - (20.+7./16.)*2.54
Yback = 0.5*(columns['scm'][lastIndex-1] - columns['scm'][0])
ybr = Yback + bar/2. + bar - (12.+7./8.)*2.54
ybl = Yback + bar/2. + bar - (12.+15./16.)*2.54
Yfront = 0.5*(columns['scm'][0] - columns['scm'][lastIndex-1])
yfr = Yfront - bar/2. - bar + (2.5)*2.54
yfl = Yfront - bar/2. - bar + (2.5)*2.54
Ztop = (3.5 - 27.8625)*2.54
Zbot = Ztop - (0.4)*2.54
Hodos = {}
Hodos['H5'] = [ [xbr,ybr,Ztop], [xfr,yfr,Ztop], [xfl,yfl,Ztop], [xbl,ybl,Ztop], \
                [xbr,ybr,Zbot], [xfr,yfr,Zbot], [xfl,yfl,Zbot], [xbl,ybl,Zbot] ]
    


                                                

############# draw a bunch of stuff
############# and write graphs, hists, etc. to root file
rfn = 'vP.root'
rf = TFile(rfn,'RECREATE')


if calcHorizDisp:
    title = 'absdX wrt nominal vs Ymeas'
    name = title.replace(' ','_')
    g = gU.makeTGraph(C,D,title,name)
    gU.color(g,0,0,setMarkerColor=True)
    gU.drawGraph(g,option='AP')
    rf.WriteTObject(g)

    h = gU.makeTH1D(D,'|dX|','absdX',nx=20,xmi=6.,xma=14.)
    gU.drawMultiHists([h])
    rf.WriteTObject(h)


title = 'residual vs phi'
name = title.replace(' ','_')
g = gU.makeTGraph(Phi,Residual,title,name)
gU.color(g,0,0,setMarkerColor=True)
gU.drawGraph(g,option='AP',SetLogy=True)
rf.WriteTObject(g)

tmg = gU.makeTMultiGraph('Fit_circle_and_measured_points') 

phi0 = -math.pi
nphi = 1000
dphi = 2.*math.pi/float(nphi)
Xfit,Yfit = [],[]
Xd,Yd = [],[] # design
Xduv,Yduv = [],[] # design radius, center at average of U,V measurements
for i in range(nphi):
    phi = phi0 + float(i)*dphi
    x = x0 + Rguess*math.cos(phi)
    y = y0 + Rguess*math.sin(phi)
    Xfit.append(x)
    Yfit.append(y)
    Xd.append( 0. + Rdesign*math.cos(phi) )
    Yd.append( 0. + Rdesign*math.sin(phi) )
    if useOldUV:
        Xduv.append( Xqave + Rdesign*math.cos(phi) )
        Yduv.append( Yqave + Rdesign*math.sin(phi) )
name = 'fitted_circle'
title = name.replace('_',' ')
cR = '{0:.2f}'.format(Rguess)
cx0= '{0:.2f}'.format(x0)
cy0= '{0:.2f}'.format(y0)
formula = cR + '^{2} = (x - '+cx0+')^{2} + (y - '+cy0+')^{2}'
title += ' ' + formula
g = gU.makeTGraph(Xfit,Yfit,title,name)
gU.color(g,2,2,setMarkerColor=False,setMarkerType=False)
tmg.Add(g)
rf.WriteTObject(g)

# draw  PMTs
phi0 = -math.pi
nphi = 100
dphi = 2.*math.pi/float(nphi)
for iPMT in range(len(Xpmt)):
    xc,yc = Xpmt[iPMT],Ypmt[iPMT]
    name = title = 'PMT'+str(iPMT)
    Xarc,Yarc = [],[]
    for i in range(nphi):
        phi = phi0+float(i)+dphi
        Xarc.append( xc + Rpmt*math.cos(phi) )
        Yarc.append( yc + Rpmt*math.sin(phi) )
    g = gU.makeTGraph(Xarc,Yarc,title,name)
    gU.color(g,0+iPMT,0+iPMT,setMarkerColor=False,setMarkerType=False)
    tmg.Add(g)
    rf.WriteTObject(g)

# draw hodoscopes
icol = 10
for H in Hodos:
    XYZ = Hodos[H]
    Xh,Yh = [],[]
    for xyz in XYZ:
        Xh.append(xyz[0])
        Yh.append(xyz[1])
    name = title = H
    g = gU.makeTGraph(Xh,Yh,title,name)
    icol += 1
    gU.color(g,icol,icol,setMarkerColor=True,setMarkerType=False)
    g.SetLineColor(icol) 
    tmg.Add(g)
    rf.WriteTObject(g)

name = title = 'design'
cR = '{0:.2f}'.format(Rdesign)
formula = cR + '^{2} = x^{2} + y^{2}'
title += ' ' + formula
g = gU.makeTGraph(Xd,Yd,title,name)
gU.color(g,4,4,setMarkerColor=False,setMarkerType=False)
tmg.Add(g)
rf.WriteTObject(g)

showDesignCenteredAtUV = False
if showDesignCenteredAtUV:
    name = 'design_centered_at_uv_average'
    title = name.replace('_',' ')
    g = gU.makeTGraph(Xduv,Yduv,title,name)
    gU.color(g,5,5,setMarkerColor=False,setMarkerType=False)
    tmg.Add(g)
    rf.WriteTObject(g)

name = 'measured_points'
title = name.replace('_',' ')
g = gU.makeTGraph(X,Y,title,name)
gU.color(g,0,0,setMarkerColor=True)
g.SetLineColor(10) # white lines
tmg.Add(g)

if 0: # no need, since replaced old with new above
    name = 'new_measured_points'
    title = name.replace('_',' ')
    g = gU.makeTGraph(newX,newY,title,name)
    gU.color(g,2,2,setMarkerColor=True)
    g.SetLineColor(10) # white lines
    tmg.Add(g)

show_old_UV = useOldUV
if show_old_UV:
    name = 'UV_measured_points'
    title = name.replace('_',' ')
    g = gU.makeTGraph(Xq,Yq,title,name)
    gU.color(g,3,3,setMarkerColor=True)
    g.SetLineColor(10) # white lines
    tmg.Add(g)

name = 'new_UV_measured_points'
title = name.replace('_',' ')
g = gU.makeTGraph(Xuv,Yuv,title,name)
gU.color(g,13,13,setMarkerColor=True)
g.SetLineColor(10) # white lines
tmg.Add(g)



gU.drawMultiGraph(tmg,abscissaIsTime=False,xAxisLabel='X(cm) in table coord system',yAxisLabel='Y(cm) in table coord system')
rf.WriteTObject(g)
rf.WriteTObject(tmg)

rf.Close()
print 'Closed',rfn
