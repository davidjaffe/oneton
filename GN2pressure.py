#!/usr/bin/env python
'''
plot GN2 pressure vs day to estimate usage
data gleaned from ELOG 20160928 - 20161121
20161121
'''
import os,sys
import matplotlib.pyplot as plt
import math,numpy

pressure = [250,350,400,400,400, 400,500,550,500,600, 650,750,850,900,900, 1000,1000,1100,1100,1150, 1500,1700,1750,1800,2000]
dayofyear= [326,323,322,321,320, 319,315,314,313,312, 309,305,302,301,300, 299,298,295,294,293,      287,284,280,279,272]

x = numpy.array(dayofyear)
y = numpy.array(pressure)
order = 1

plt.clf()
plt.subplot(2,1,1)
plt.xlabel('Day of year')
plt.ylabel('GN2 Pressure (psi)')
plt.plot(x,y,'ro')
plt.ylim(0., 1.05*max(y) )
par = numpy.polyfit(x,y,order)
xf = numpy.linspace(min(x)-1.,max(x)+20.,1000) 
yf = numpy.polyval(par,xf)
plt.plot(xf,yf,'b-') 
plt.grid()

x = x-min(x)
plt.subplot(2,1,2)
plt.xlabel('Days elapsed')
plt.ylabel('GN2 Pressure (psi)')
plt.plot(x,y,'ro')
plt.ylim(0., 1.05*max(y) )

par = numpy.polyfit(x,y,order)
xf = numpy.linspace(0.,60.,1000)
yf = numpy.polyval(par,xf)
plt.plot(xf,yf,'b-') 

plt.grid()


pdf = 'Figures/GN2pressure_vs_day.pdf'
plt.savefig(pdf)
print 'GN2pressure wrote',pdf
