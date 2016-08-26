# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 08:49:02 2016

@author: lbignell
"""

from rootpy.tree import Tree
from rootpy.io import root_open
from random import gauss
import oneton.SaveAllEvts
import imp
import numpy as np
import ROOT

def simpletree():
    f = root_open("testrootpy.root", "recreate")

    tree = Tree("test")
    tree.create_branches(
        {'x': 'F',
         'y': 'F',
         'z': 'F',
         'i': 'I',
         'QDC': 'F'})

    for i in xrange(10000):
        tree.x = gauss(.5, 1.)
        tree.y = gauss(.3, 2.)
        tree.z = gauss(13., 42.)
        tree.i = i
        tree.QDC = [gauss(1.,1.), gauss(2.,1.), gauss(3.,1.), gauss(4.,1.)]
        tree.fill()
    tree.write()

    f.close()

def cleararrs(obj):
    for k in range(8):
        for j in range(32):
            prefix = 'Ch{0}_'.format(k)
            ch = '[{0}]'.format(j)
            setattr(obj, prefix + 'area' + ch, 0)
            setattr(obj, prefix + 'time' + ch, 0)
            setattr(obj, prefix + 'ped' + ch, 0)
            setattr(obj, prefix + 'pedsd' + ch, 0)
    
def realtree():
    imp.reload(oneton.SaveAllEvts)
    obj = oneton.SaveAllEvts.rootpyEvts('testrootpy_realtree_3.root')
    #area = ROOT.std.vector(float)()
    #time = ROOT.std.vector(float)()
    #ped = ROOT.std.vector(float)()
    #pedsd = ROOT.std.vector(float)()    
    trigtypes = ['CT ', 'M ', 'LED ', 'CT M ', 'CT LED ', 'M LED ', 'CT M LED ']
    for i in range(100):
        obj.tree.rnum = i
        obj.tree.evtnum = i
        obj.tree.temp = np.random.normal(25., 0.2)
        obj.tree.time = i
        obj.tree.QDC1 = [np.random.normal(j, 0.1) for j in range(8)]
        obj.tree.QDC2 = [np.random.normal(j, 0.1) for j in range(8)]
        obj.tree.scaler = [i for j in range(16)]
        nch = 200
        while nch>16: nch = np.random.poisson(5) 
        obj.tree.TDCnumch = nch
        #obj.tree.WFDnpulses = [nch]*8
        #obj.tree.WFDnsubp = [1]*8
        for j in range(nch):
            obj.tree.TDCch[j] = j
            obj.tree.TDCval[j] = np.random.normal(j, 0.1)
            #for k in range(8):
            #    prefix = 'Ch{0}_'.format(k)
            #obj.tree.Ch0_area[j] = np.random.normal(j, 0.1)
            #obj.tree.Ch0_time[j] = np.random.normal(j, 0.1)
        for k in range(8):
            prefix = 'Ch{0}_'.format(k)
            #setattr(obj.tree, prefix + 'area', [np.random.normal(j, 0.1) for j in range(32)])            
            #setattr(obj.tree, prefix + 'time', [np.random.normal(j, 0.1) for j in range(32)])            
            #setattr(obj.tree, prefix + 'ped', [np.random.normal(100*j, 0.1) for j in range(32)])            
            #setattr(obj.tree, prefix + 'pedsd', [np.random.normal(j, 0.1) for j in range(32)])            
            for j in range(nch):
                pass
                #prefix = 'Ch{0}_'.format(k)
                #ch = '[{0}]'.format(j)
                #setattr(obj.tree, prefix + 'area' + ch, np.random.normal(j, 0.1))
                #setattr(obj.tree, prefix + 'time' + ch, np.random.normal(j, 0.1))
                #setattr(obj.tree, prefix + 'ped' + ch, np.random.normal(100*j, 0.1))
                #setattr(obj.tree, prefix + 'pedsd' + ch, np.random.normal(j, 0.1))
                #val = getattr(obj.tree, prefix + 'area')
                #print(prefix + 'area' + ' = {0}'.format(val))
        #        area.push_back(np.random.normal(k, 0.1))
        #        time.push_back(np.random.normal(k, 0.1))
        #        ped.push_back(np.random.normal(100*k, 0.1))
        #        pedsd.push_back(np.random.normal(k, 0.1))
        #    obj.tree.WFDarea.push_back(area)
        #    obj.tree.WFDtime.push_back(time)
        #    obj.tree.WFDped.push_back(ped)
        #    obj.tree.WFDpedsd.push_back(pedsd)
        obj.tree.trigtype = trigtypes[np.random.choice(len(trigtypes))].ljust(9)
        obj.tree.Fill()
        cleararrs(obj)
        
    obj.rfile.Write()
    obj.rfile.Close()
    return