#!/usr/bin/env python
'''
calibrate thermocouple based on data of 20160202
'''


import graphUtils
import sys

class calibThermocouple():
    def __init__(self,graphIt=False):
        oven = [25.,  30., 40., 35., 30., 25., 15., 20.]
        tc   = [36.1, 40.5, 49.4, 45.2, 41.2, 36.7, 27.1, 31.1]
        if graphIt:
            gU = graphUtils.graphUtils()

            tmg = gU.makeTMultiGraph('Fit_line_and_measured_points')

            g = gU.makeTGraph(tc,oven,'thermocouple readback vs oven temp (C)','tc_v_oven')
            gU.color(g,0,0,setMarkerColor=True)
            g.SetLineColor(10) # white lines
            #gU.drawGraph(g)
            tmg.Add(g)


        x = tc
        y = oven
        Sx = sum(x)
        Sy = sum(y)
        Sxy = sum([a*b for a,b, in zip(x,y)])
        Sxx = sum([a*a for a in x])
        N = float(len(x))
        m = (Sx*Sy - N*Sxy)/(Sx*Sx - N*Sxx)
        b = (Sy - m*Sx)/N
        self.slope = m
        self.intercept = b

        if graphIt:
            U,V = [],[]
            for i in range(250,500):
                X = float(i)/10.
                Y = m*X+b
                U.append(X)
                V.append(Y)

            title = 'Tcalib = {0:.2f}*Treadback + {1:.2f}'.format(m,b)
            g = gU.makeTGraph(U,V,title,'fit')
            tmg.Add(g)
            gU.drawMultiGraph(tmg,abscissaIsTime=False,xAxisLabel='Thermocouple readback (C)', yAxisLabel='Oven temp (C)')
        return
    def getCalibTC(self,TCreadback):
        return TCreadback*self.slope + self.intercept

if __name__ == '__main__' :
    graphIt = False
    if len(sys.argv)>1 and sys.argv[0]=='calibThermocouple.py': graphIt = True
    r = calibThermocouple(graphIt=graphIt)
    TCin = 35.
    TC = r.getCalibTC(TCin)
    print 'calibThermocouple.getCalibTC returns',TC,'(C) for readback',TCin,'(C)'
