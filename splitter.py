#!/usr/bin/env python
'''
splitter and amplifier measurements
20151228

splitter msmts from book#2,p39
amp gain msmts from book#2,p81 and book#3,p64
'''
import math

print 'Measured splitter fractions'
splitter = {'s1': [860.,104.],
            's2': [840.,112.],
            's3': [870.,112.],
            's4': [860.,120.],
            's5': [840.,104.],
            's6': [850.,108.],
            's7': [870.,104.],
            's8': [860.,104.],
            's9': [850.,108.]
            }
S,rms = 0.,0.
for s in sorted(splitter):
    r,q = splitter[s]
    t = r+q
    fr,fq = r/t,q/t
    S += fr
    rms += fr*fr
    #print s,'fractions r,q',fr,fq
    print '{0:} fractions r,q {1:.3f} {2:.3f}'.format(s,fr,fq)
N = float(len(splitter))
S = S/N
rms = rms/N
rms = math.sqrt( (rms-S*S)/(N-1.) )
print 'mean fraction r {0:.3f} +- {1:.3f}'.format(S,rms)

print '\n Measured amplification'
ampinput = [0.5*(62.0 + 61.5) for i in range(8)]
ampinput.extend( [0.5*(36.0+36.8) for i in range(8)] )
amp = [123.]
amp.extend( [124. for i in range(5)] )
amp.append( 125. )
amp.append( 124. )
amp.extend( [1830., 1820., 1840., 1850., 1860., 1830., 1860., 1830.] )
j = 0
for i,a in enumerate(amp):
    denom = ampinput[i]
    if i>7:
        j += 1
        s = 's'+str(j)
        r,q = splitter['s'+str(j)]
        fr = r/(r+q)
        q = q/(r+q)*a/denom
        ratio = q/fr
    else:
        s = 'NA'
        fr = 1.
        q = a/denom
        ratio = 0.
    print 'chan{0:2d} ampl {1:.2f} frac*gain({2:}) {3:.2f} 1-frac({2:}) {4:.2f} ratio {5:.3f}'.format(i,a/denom,s,q,fr,ratio) 
