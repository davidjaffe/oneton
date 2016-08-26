#!/usr/bin/env python
'''
return positions of acrylic vessel, pmts, hodoscopes, leds in
platform coordinate system see vesselPosition.py
units are centimeters
20160527
'''

import sys
import os

class getPosition():
    def __init__(self):
        modulepath = os.path.dirname(__file__)
        self.data_file = os.path.join(modulepath, 'Positions.dat')
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
                    if debug: print(line[:-1])
                    N += 1
            f.close()
            print('getPositions.getData Read',N,'data items from',self.data_file)
            if debug:
                for nme in sorted(self.data): print(self.data[nme])

        if name in self.data:

            d = [float(x) for x in self.data[name][1:]]
            return d
        else:
            sys.exit('getPositions.getData ERROR name ' + name + ' is not present in data file')
if __name__ == '__main__':
    name = ['S0']

    if len(sys.argv)>1: name = sys.argv[1:]
    gP = getPosition()
    for n in name:
        d = gP.getData(n)
        print(n,d)
