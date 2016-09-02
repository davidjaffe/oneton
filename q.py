#!/usr/bin/env python
'''
estimate offset of emission M/C from pre-calibration emission scan 
/Users/djaffe/work/GIT/QY/Henry/Emission/Empty/empty_box_frosted_ex350nm_em330-370nm_2sec_20160901.txt
20160902
'''
import os,sys
import graphUtils
import ROOT
import math

# these data taken from file given above by cut and paste
a = [348,  3994.47,      348,  26470.9648,
    348.5,9051.292,     348.5,59976.2656,
    349,  25886.9258,   349,  171519.75,
    349.5,63846.3672,   349.5,423386.219,
    350,  152923.969,   350,  1013791.5,
    350.5,286784.6,     350.5,1901791.75,
    351,  466242.6,     351,  3090916.5,
    351.5,638382.1,     351.5,4232679,
    352,  759341.7,     352,  5033610.5,
    352.5,778043.563,   352.5,5156592.5,
    353,  699284.9,     353,  4637469,
    353.5,562220.063,   353.5,3724572,
    354,  359466.875,   354,  2381754.5,
    354.5,210926.734,   354.5,1396357.75,
    355,  98464.57,     355,  652106.25,
    355.5,40051.5469,   355.5,265176.344,
    356,  12847.6436,   356,  85048.0547]
l = len(a)
wl,raw,corr = [],[],[]
for i in range(0,l,4):
    wl.append(a[i])
    raw.append(a[i+1])
    corr.append(a[i+3])

print wl
print raw
print corr

wlmin = 352.4-2.
wlmax = 352.4+2.

srw,sr,scw,sc = 0.,0.,0.,0.
for i,x in enumerate(wl):
    if wlmin<=x and x<=wlmax:
        r = raw[i]
        c = corr[i]
        srw += r*x
        sr  += r
        scw += c*x
        sc  += c
print 'weighted means. raw',srw/sr,'corr', scw/sc,'wlmin,wlmax=',wlmin,wlmax
