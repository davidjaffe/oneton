#!/usr/bin/env python
'''
quick study of hist with time axis
20160211
'''

import numpy
import graphUtils
import ROOT
from ROOT import TFile,TH1D,TH2D,gDirectory
import sys
import datetime
import random

import graphUtils

gU = graphUtils.graphUtils()
fmt = '%H:%M:%S %m/%d/%y'

now = datetime.datetime.now()
cnow = now.strftime(fmt)
print 'now,cnow',now,cnow

tStart = '12:05:49 01/11/16'
tEnd   = '23:59:59 02/22/16'
tS = gU.getdatetime(tStart,fmt=fmt)
tE = gU.getdatetime(tStart,fmt=fmt)
print 'from getdatetime',tS,tE
tS = datetime.datetime.strptime(tStart,fmt)
tE = datetime.datetime.strptime(tEnd,fmt)
print 'from datetime',tS,tE
tS = gU.getTDatime(tStart,fmt=fmt)
tE = gU.getTDatime(tEnd,fmt=fmt)
print 'from getTDatime',tS,tE

h = TH1D('h','h',100,tS,tE)
for i in range(1000):
    t = tS + (tE-tS)*random.random()
    h.Fill(t)

gU.drawMultiHists([h],abscissaIsTime=True)
