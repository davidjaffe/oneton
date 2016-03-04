#!/usr/bin/env python
'''
test compression algorithms for h5py
test use of ZipFile, gzip
20160209
'''

import h5py
import numpy
import random
from uuid import uuid4
import zipfile
import gzip, shutil

testZip = False
testgzip = True

data = {}
Idata = {}
Sdata = {}
iRange = [10000,100001]
if testZip or testgzip: iRange = [10000]
for i in iRange:
    data[i] = numpy.random.sample(i)
    Idata[i] = numpy.random.random_integers(99999999,high=None,size=i)
    Sdata[i] = []
    for j in range(i): Sdata[i].append( str(uuid4()) )

itypes = ['F','I','S']
dtypes = [data,Idata,Sdata]
CompAlgs = [None,"gzip","lzf"]
if testZip or testgzip: # only one file with no compression
    itypes = ['S']
    dtypes = [Sdata]
    CompAlgs = [None]
        
    
for what,DATA in zip(itypes, dtypes):
    
    for compalg in CompAlgs:
        for i in DATA:
            label = 'none'
            if compalg is not None: label = compalg
            label += '_'+what+str(i)
            print label
            fn = 'TESTH5/'+label+'.h5'
            f = h5py.File(fn,'w')
            if compalg=='gzip':
                f.create_dataset(label,data=DATA[i],compression=compalg,compression_opts=9)
            else:
                f.create_dataset(label,data=DATA[i],compression=compalg)

            f.close()

            if testZip:
                zfn = fn.replace('.h5','.zip',zipfile.ZIP_DEFLATED)
                zf = zipfile.ZipFile(zfn,'w')
                zf.write(fn)
                zf.close()
                print 'wrote',zfn,'from',fn

            if testgzip:
                gzfn = fn+'.gz'
                with open(fn,'rb') as fin, gzip.open(gzfn,'wb') as fout:
                    shutil.copyfileobj(fin,fout)
                print 'wrote',gzfn,'from',fn
                newfn = fn.replace(label,'new'+label)
                with gzip.open(gzfn,'rb') as fin, open(newfn,'wb') as fout:
                    shutil.copyfileobj(fin,fout)
                print 'wrote',newfn,'from',gzfn

