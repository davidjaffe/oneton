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

#import sys,pybindgen
#from pybindgen import ReturnValue, Parameter, Module, Function, FileCodeSink

#def vector_gen(file_out):
#    '''
#    this function creates bindings to c++ STL vectors using pybindgen.
#    '''
#    mod=pybindgen.Module('vectors')
#    mod.add_include('"loadvec.h"')
#    veci = mod.add_container('std::vector<int>', 'int',
#                             'vector',custom_name="VecI")
#    vecd = mod.add_container('std::vector<double>', 'double',
#                             'vector',custom_name="VecD")
#    vecveci = mod.add_container('std::vector< std::vector<int> >',
#                                'std::vector<int>',
#                                'vector', custom_name='VecVecI')
#    vecvecf = mod.add_container('std::vector< std::vector<double> >',
#                                'std::vector<double>',
#                                'vector', custom_name='VecVecI')
#    mod.generate(file_out)

class TDC(rootpy.tree.TreeModel):
    TDCnumch = IntCol()
    TDCch = IntArrayCol(16, length_name='TDCnumch')#ROOT.std.vector('int')#stl.vector('int')
    TDCval = FloatArrayCol(16, length_name='TDCnumch')#ROOT.std.vector('float')#stl.vector('float')

class WFDproc(rootpy.tree.TreeModel):
    area = FloatArrayCol(32)#ROOT.std.vector(float)()
    time = FloatArrayCol(32)#ROOT.std.vector(float)()

class WFD(rootpy.tree.TreeModel):
    WFDped = FloatArrayCol(8)#ROOT.std.vector(float)()
    WFDpedsd = FloatArrayCol(8)#ROOT.std.vector(float)()
    
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
        #vector_gen(sys.stdout)
        self.AllTrigsfname = fname
        self.rfile = rootpy.io.root_open(self.AllTrigsfname, 'recreate')
        self.tree = Tree('evttree'+self.AllTrigsfname, model=Trigger)
        return

class ROOTevts():
    '''
    This class will create and handle the ROOT event tree object.
    '''
    def __init__(self, fname):
        '''
        Initialize the object, creating the ROOT file fname and the tree.
        '''
        self.maxlenTDC = 16   
        self.maxNPulse = 16
        #The code below is burdensome, maybe I should use rootpy?
        ROOT.gROOT.ProcessLine("\
        struct digitizer{\
            Int_t \
        };")

        self.AllTrigsfname = fname
        self.rfile = ROOT.TFile(self.AllTrigsfname, 'recreate')
        self.tree = ROOT.TTree('evttree'+self.AllTrigsfname,
                                    'Tree containing all triggers for ' + self.AllTrigsfname)
        self.rnum = np.array([0], dtype=int)
        self.evtnum = np.array([0], dtype=int)
        self.temp = np.array([0], dtype=float)
        self.time = np.array([0], dtype=float)
        self.QDC2 = np.array([0]*8, dtype=float)
        self.QDC1 = np.array([0]*8, dtype=float)
        self.scaler = np.array([0]*16, dtype=int)
        self.lenTDC = np.array([0], dtype=int)
        #self.TDC = np.zeros(self.maxlenTDC, dtype='int, float')
        #self.TDC = np.array([[0.]*self.maxlenTDC, [0.]*self.maxlenTDC], dtype=float)
        self.TDCch = np.array([0]*self.maxlenTDC, dtype=int)
        self.TDCval = np.array([0.]*self.maxlenTDC, dtype=float)
        self.NPulse = np.array([0], dtype=int)
        self.QPulse = np.array([0]*self.maxNPulse, dtype=float)
        self.TPulse = np.array([0]*self.maxNPulse, dtype=float)
        #create the trees        
        self.tree.Branch('rnum', self.rnum, 'rnum/I')
        self.tree.Branch('evtnum', self.evtnum, 'evtnum/I')
        self.tree.Branch('temp', self.temp, 'temp/D')
        self.tree.Branch('time', self.time, 'time/D')
        self.tree.Branch('QDC1', self.QDC1, 'QDC1[8]/D')
        self.tree.Branch('QDC2', self.QDC2, 'QDC2[8]/D')
        self.tree.Branch('scaler', self.scaler, 'scaler[16]/I')
        self.tree.Branch('lenTDC', self.lenTDC, 'lenTDC/I')
        self.tree.Branch('TDCch', self.TDCch, 'TDCch[lenTDC]/I')
        self.tree.Branch('TDCval', self.TDCval, 'TDCval[lenTDC]/D')
        self.tree.Branch('NPulse', self.NPulse, 'NPulse/I')
        self.tree.Branch('QPulse', self.QPulse, 'QPulse[NPulse]/D')
        self.tree.Branch('TPulse', self.TPulse, 'TPulse[NPulse]/D')
        return