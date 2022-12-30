#!/usr/bin/env python
'''
plot and integrate  time profile distribution from 
Cherenkov and scintillation separation in water-based liquid scintillator using an LAPPD Figure12 arXiv:2110.13222
Left figure 12 corresponding to 1% WbLS

20221230
'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import numpy

class timeprofile():
    def __init__(self):

        # 1750 counts corresponds to 224 mm
        # I measured depth of each point, so subtract measured values from 224 and multiply by 1750/224 to get number of counts.
        # time bins are 0.1 ns, low edge is (arbitrary) t=0ns
        self.Cmeas = [223, 223, 222, 223, 220, 217, 211, 203, 186, 172,\
                      150, 136, 114, 90, 62.8, 54.7,59.7,117, 155, 190,\
                      198, 202, 210, 211, 213, 214, 218, 220, 222, 221, 224]
        self.Smeas = [223, 223, 223, 224, 224, 224, 223, 221, 222, 219,\
                      216, 211, 210, 205, 196, 183, 179, 165, 165, 153,\
                      162, 154, 157, 159, 157, 158, 165, 161, 166, 165,\
                      172, 177, 172, 169, 181, 182, 184, 181, 185, 189,\
                      191, 197, 194, 199, 195, 197, 200, 198, 199, 201,\
                      202, 205, 208, 208, 207, 205, 210, 211, 209, 212,\
                      215, 212, 214, 211, 214, 213, 215, 213, 218, 217,\
                      215, 216, 217, 216, 217, 218, 220, 220, 218, 220, \
                      218, 219, 218, 219, 220, 220, 221, 221, 220, 220]
        q = 1750./224.
        self.Ccts = [(224-x)*q for x in self.Cmeas]
        self.Scts = [(224-x)*q for x in self.Smeas]
        tmi = 0.1/2.
        dt = 0.1
        self.Ctime = []
        for i in range(len(self.Ccts)):
            self.Ctime.append( tmi + float(i)*dt)
        self.Stime = []
        for i in range(len(self.Scts)):
            self.Stime.append( tmi + float(i)*dt)

        if len(self.Ccts)!=len(self.Ctime): sys.exit('ERROR lengths not equal for Ccts and Ctime')

            
        return

    def main(self,model='Caravaca'):

        # Cerenkov = blue, Scintillation = orange
        Ct,Cc = numpy.array(self.Ctime), numpy.array(self.Ccts)
        St,Sc = numpy.array(self.Stime), numpy.array(self.Scts)
        eCc = numpy.sqrt(Cc)
        eSc = numpy.sqrt(Sc)

        Tc = []
        Tt = []
        for i,t in enumerate(St):
            S = Sc[i]
            C = 0.
            if i<len(Cc) : C = Cc[i]
            Tt.append(t)
            Tc.append(S+C)

        Tt,Tc = numpy.array(Tt),numpy.array(Tc)
        eTc = numpy.sqrt(Tc)

        sumC = numpy.sum(Cc)
        sumS = numpy.sum(Sc)
        sumT = numpy.sum(Tc)
        print('timeprofile.main Sum of counts Cerenkov {:.0f} Scintillation {:.0f} Total {:.0f}'.format(sumC,sumS,sumT))
        print('timeprofile.main Fractions Cerenkov {:.3f} Scintillation {:.3f}'.format(sumC/sumT,sumS/sumT))
        fC = ', fraction {:.3f}'.format(sumC/sumT)
        fS = ', fraction {:.3f}'.format(sumS/sumT)
            
        plt.errorbar(Tt,Tc,fmt='o',yerr=eTc,color='red',label='MC')
        plt.errorbar(Ct,Cc,fmt='o',yerr=eCc,color='blue',label='MC Cerenkov'+fC)
        plt.errorbar(St,Sc,fmt='o',yerr=eSc,color='orange',label='MC Scintillation'+fS)

        plt.title('arXiv:2110.13222 Cherenkov and scintillation separation in \n1% WbLS using an LAPPD Figure12')
        plt.grid()
        plt.xlim(0, 10.)
        plt.legend(loc='best')
        png = 'TIMEPROFILE_FIGURES/timeprofile_arXiv2110.13222.png'
        plt.savefig(png)
        print('timeprofile.main Wrote',png)
        plt.show()

        
        return
if __name__ == '__main__' :
    P = timeprofile()
    P.main()
    
