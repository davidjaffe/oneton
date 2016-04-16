# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 09:47:35 2016

@author: lbignell
"""
import numpy as np
import ROOT
import rootpy
from rootpy.tree import IntCol, FloatCol, FloatArrayCol, IntArrayCol, CharCol, CharArrayCol
from rootpy.tree import Tree
from rootpy import stl

class TDC(rootpy.tree.TreeModel):
    TDCnumch = IntCol()
    TDCch = IntArrayCol(16, length_name='TDCnumch')
    TDCval = FloatArrayCol(16, length_name='TDCnumch')

class WFDproc(rootpy.tree.TreeModel):
    area = FloatArrayCol(32)
    time = FloatArrayCol(32)

class WFD(rootpy.tree.TreeModel):
    WFDped = FloatArrayCol(8)
    WFDpedsd = FloatArrayCol(8)
    
class Trigger(TDC, WFD,
              WFDproc.prefix('S0_'), WFDproc.prefix('S1_'), WFDproc.prefix('S2_'),
              WFDproc.prefix('S3_'), WFDproc.prefix('S4_'), WFDproc.prefix('S5_'),
              WFDproc.prefix('S6_'), WFDproc.prefix('S7_')):
    '''
    This class defines the rootpy event tree model.
    '''
    rnum = IntCol()
    evtnum = IntCol()
    temp = FloatCol()
    time = FloatCol()
    QDC1 = FloatArrayCol(16)
    QDC2 = FloatArrayCol(16)
    scaler = IntArrayCol(16)
    
class rootpyEvts():
    def __init__(self, fname):
        '''
        Initialize the object, creating the ROOT file fname and the tree.
        '''
        self.AllTrigsfname = fname
        self.rfile = rootpy.io.root_open(self.AllTrigsfname, 'recreate')
        self.tree = Tree('evttree', model=Trigger)
        return
