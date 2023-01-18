#!/usr/bin/env python
'''
plot attenuation length measurements
20230118
results from 1508.07029
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit

class atten():
    def __init__(self):

        self.debug = 0
        
        # 1508.07029 contents of Figure 1 (dubbed 'LB' for Lindsey Bignell)
        # wavelength : 770 pixels = 205 nm starting at 350 nm
        # attenuation length : 499,399,299,200,99,0 pixels = 1e-3,1e-2,1e-1,1e0,1e1,1e2 meters
        # z = 10^(m*y + b) where z = attlen in meters and y is in pixels
        y1,z1 = 0.,1.e-3
        y2,z2 = 499.,1.e2
        m = math.log10(z2/z1)/(y2-y1)
        b = math.log10(z1) - m*y1
        if self.debug > 0 : 
            print('atten.__init__ m {:.3f} b {:.3f}'.format(m,b))
            for y in [0.,99.,200.,299.,399.,499.]:
                print('atten.__init__ pixels {:.1f} z(m) {:6.3f}'.format(y,math.pow(10.,m*y+b)))
        # figure contains attenuation for 2 different WbLS formulations
        self.LBpixel_wl = numpy.array([  0,  20,  40,  60,  80, 100, 130, 160, 200, 300, 450, 500, 550, 600, 750],dtype=float)
        self.LBpixel_al1= y2 - numpy.array([140, 113,  94,  86,  80,  76,  71,  63,  56,  44,  35,  42,  44,  48,  81],dtype=float)
        self.LBpixel_al2= y2 - numpy.array([443, 377, 344, 325, 310, 298, 274, 253, 243, 227, 201, 189, 176, 165, 138],dtype=float)
        self.LB_wl = 350. + self.LBpixel_wl*205./770.
        self.LB_al1 = numpy.power(10.,m*self.LBpixel_al1+b)
        self.LB_al2 = numpy.power(10.,m*self.LBpixel_al2+b)
        print('{:6} {:7} {:7} from Figure 1 of  1508.07029 atten.__init__'.format(' nm ','L1(m)','L2(m)'))
        for i,wl in enumerate(self.LB_wl):
            al1,al2 = self.LB_al1[i],self.LB_al2[i]
            px1,px2 = self.LBpixel_al1[i],self.LBpixel_al2[i]
            print('{:6.1f} {:7.4f} {:7.4f}'.format(wl,al1,al2))


        ### Aiwu's 20180426 results
        ### 'scattering correction'
        ### wavelength (nm) = 350 nm + pixels*(400nm/882)
        ### attlen(m) = (750 - pixels)*12m
        self.Apixel_wl1 = numpy.array([000,  36,  52, 112, 222, 317, 398, 486, 518, 539, 567, 598, 663],dtype=float)
        self.Apixel_al1 = numpy.array([750, 637, 563, 460, 269, 146, 180, 217, 346, 434, 498, 527, 545],dtype=float)
        self.A_wl1 = 350. + (400./882.)*self.Apixel_wl1
        self.A_al1 = (750. - self.Apixel_al1)*12./750.
        print('{:6} {:7} from page 5 of Aiwu 20180426 results for `scattering correction` atten.__init__'.format(' nm ','L(m)'))
        for wl,al in zip(self.A_wl1,self.A_al1):
            print('{:6.1f} {:7.4f}'.format(wl,al))

        ### Aiwu's 20180426 results from 1cm vs 10cm measurements
        ### wavelength (nm) = 340 nm + pixels*(300 nm/740)
        ### attlen(m) = (628 - pixels)* 10m/628
        self.Apixel_wl2 = numpy.array([ 26,  55,  75, 140, 239, 309, 354, 375, 392, 405, 427, 429, 457, 463, 538, 583, 618, 730, 845],dtype=float)
        self.Apixel_al2 = numpy.array([628, 565, 462, 354, 188,  87,  55,  47,  50,  90,  90,  79,  78, 105, 150, 298, 432, 461, 528],dtype=float)
        self.A_wl2 = 340. + self.Apixel_wl2*(300./740.)
        self.A_al2 = (628. - self.Apixel_al2)*(10./628.)
        print('{:6} {:7} from page 8 lower right of Aiwu 20180426 results for 1cm v 10cm cell msmt atten.__init__'.format(' nm ','L(m)'))
        for wl,al in zip(self.A_wl2,self.A_al2):
            print('{:6.1f} {:7.4f}'.format(wl,al))
            
        return


    def main(self):


        for yax in ['linear','log']:
            ## set wavelength and atten length ranges
            x1,x2 = 300.,700.
            y1,y2 = 0., 13.5
            if yax=='log' : y1,y2 = 1.e-3, 1.e2

            x = self.LB_wl
            y = self.LB_al1
            yname = 'arXiv:1508.07029 improved'
            plt.plot(x,y,marker='o',linestyle='-',color='black',label=yname)
            y = self.LB_al2
            yname = 'arXiv:1508.07029 original'
            plt.plot(x,y,marker='o',linestyle='-',color='red',label=yname)

            x = self.A_wl1
            y = self.A_al1
            yname = 'Aiwu 20180426 scat corr'
            plt.plot(x,y,marker='o',linestyle='-',color='blue',label=yname)
            x = self.A_wl2
            y = self.A_al2
            yname = 'Aiwu 20180426 1cm v 10cm cell'
            plt.plot(x,y,marker='o',linestyle='-',color='green',label=yname)


            plt.title('WbLS attenuation lengths')
            plt.grid()
            plt.legend(loc='best')
            plt.xlim(x1,x2)
            plt.ylim(y1,y2)
            plt.yscale(yax)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Attenuation length (m)')
            if 1 : 
                png = 'ATTEN_FIGURES/atten_WbLS_'+yax+'.png'
                plt.savefig(png)
                print('atten.main Wrote',png)
            plt.show()

        
        return
if __name__ == '__main__' :
    P = atten()
    P.main()
    
