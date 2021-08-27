#!/usr/bin/env python
'''
return positions of acrylic vessel, pmts, hodoscopes, leds in
platform coordinate system see vesselPosition.py
units are centimeters
20160527
20210826 output useful info for final technote
'''

import sys


class getPosition():
    def __init__(self):
        self.data_file = 'Positions.dat'
        self.data = {}
        return
    def getData(self,name):
        '''
        return data as list for data corresponding to name in data_file
        on first call, fill dict with contents of data_file
        '''
        debug = False
        if len(self.data)==0: 
            f = open(self.data_file,'r')
            N = 0
            for line in f:
                if line[0]!='*':
                    nme = line.split(' ')[0]
                    if nme in self.data:
                        sys.exit('getPositions.getData ERROR Duplicate name ' + nme + ' in ' + self.data_file)
                    self.data[nme] = line[:-1].split(' ') # strip control character
                    if debug: print line[:-1]
                    N += 1
            f.close()
            print 'getPositions.getData Read',N,'data items from',self.data_file
            if debug:
                for nme in sorted(self.data): print self.data[nme]

        if name in self.data:

            d = [float(x) for x in self.data[name][1:]]
            return d
        else:
            sys.exit('getPositions.getData ERROR name ' + name + ' is not present in data file')
if __name__ == '__main__':
    name = ['S0']
    # for H0-H5, from Positions.dat
    # *Box name, x,y,z of center, half-height,-width,-length. Units cm
    name = ['H0', 'H1', 'H2', 'H3', 'H4', 'H5']

    if len(sys.argv)>1: name = sys.argv[1:]
    gP = getPosition()
    for n in name:
        d = gP.getData(n)
        print n,d
        # convert half-dimensions to full dim and convert to inch
        cmD   = []
        inchD = []
        for x in d[3:]:
            z = 2.*x
            cmD.append(z)
            z = 2.*x/2.54
            inchD.append(z)
        print n,cmD,inchD
        s = '{0} '.format(n)
        s+= '{0:.2f} {1:.2f} {2:.2f} cm '.format(*cmD)
        s+= '{0:.2f} {1:.2f} {2:.3f} inch '.format(*inchD)
        print s
