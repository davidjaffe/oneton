#!/usr/bin/env python
'''
20230103 compare, plot, etc. a number of results on relative and absolute light yield

Ref1 - 2006.00173 : WbLS aLY (LAB-based) and LABPPO(2 g/L) aLY
Ref2 - 1508.07029 :  WbLS aLY (PC-based) and  LABPPO(2 g/L) aLY
Ref3 - NIM A 701 (2013) 133–144 : Relative LY of LAB+PPO, DIN+PPO to PC+PPO
Ref4 - 2205.15046 : LAB+4g/L PPO, DIN+4g/L PPO, p-Xylene LY, EJ301 relative to anthracene
Ref5 - 2003.10491 : LAB + g/L PPO relative to EJ301

'''
import numpy

import os

import sys
import math
import matplotlib.pyplot as plt
import numpy
from scipy.optimize import curve_fit

class LS_yield():
    def __init__(self):

        self.figDir = 'LS_YIELD_FIGURES/'

        # 1508.07029 WbLS LY contents of table 4
        self.BNLabsLY = {'0.4%': (19.9, 2.3),
                         '1%'  : (108.9, 10.9),
                        'LABPPO':(9156,917)}
        
        
        # 2006.00173 contents of table 1 absLY[concentration] = (central value, unc)
        self.absLY = {'1%': (234., 30.),
                      '5%': (770., 72.),
                      '10%':(1357,125.),
                      'LABPPO':(11076, 1004)}
        self.conc = {'0.4%':0.4,
                     '1%': 1., 
                     '5%':5., 
                     '10%':10., 
                     'LABPPO':100.}
        # calculate LY relative to LAPPPO assuming uncertainties are independent (add in quad)
        self.relLY = {}
        keyden = 'LABPPO'
        den,eden = self.absLY[keyden]
        for key in self.absLY:
            if key!=keyden:
                num,enum = self.absLY[key]
                r = num/den
                er = r*math.sqrt(enum*enum/num/num + eden/eden/den/den)
                self.relLY[key] = (r,er)
                print('LS_yield.__init__ LY('+key+') relative to keyden {:.4f}({:.4f})'.format(r,er))

        self.relBNLLY = {}
        for key in self.BNLabsLY:
            if key!=keyden:
                num,enum = self.BNLabsLY[key]
                r = num/den
                er = r*math.sqrt(enum*enum/num/num + eden/eden/den/den)
                self.relBNLLY[key] = (r,er)
                print('LS_yield.__init__ BNL LY('+key+') relative to keyden {:.4f}({:.4f})'.format(r,er))

                
        
        # section 3.1.3 reported results of linear fit to 3 wbls points
        self.lfit = {'slope': (127.9,17.0),
                     'intercept': (108.3, 51.0)}

        ### 1105.2100 = Ref6, H. Wan Chan Tseung, J.Kasper, N.Tolich
        ### measurements of LY vs electron energy (Ee) taken from Figure 7
        ### for LAB + 3 g/L PPO and EJ301
        self.Ref6_name = 'LAB30'
        self.Ref6_conc = '3.0 g/L PPO'
        self.Ref6_relTo= 'EJ301'
        ### LY and Ee in pixel units
        ### LAB CALIBRATION
        ### LY = (800-pLY) * (4e3/(800-163)) = optical photons
        ### Ee = pEe * 3.0 MeV / 456
        self.LAB_pLY = numpy.array([301, 365, 434, 474, 516, 567, 610, 663, 704, 722, 747, 756, 764, 773, 783])
        self.LAB_pEe = numpy.array([420, 373, 314, 282, 249, 207, 169, 130,  93,  79,  60,  51,  44,  33,  23])
        self.LAB_LY = (800.-self.LAB_pLY) * (4.e3/(800.-163.))
        self.LAB_Ee = self.LAB_pEe * 3.0 /456.
        words = 'LS_yield.__init__ LAB (Ee,LY) ='
        for Ee,LY in zip(self.LAB_Ee,self.LAB_LY):
            words += ' ({:.2f},{:.0f})'.format(Ee,LY)
        print(words)

        ### EJ301 CALIBRATION
        ### LY = (800-pLY) * (3.5e3/(800-91)) = optical photons
        ### Ee = pEe * 2.5 MeV / 653.
        self.EJ301_pLY = numpy.array([ 82, 235, 428, 611, 645, 679, 692, 706, 717, 769, 785])
        self.EJ301_pEe = numpy.array([541, 422, 290, 157, 132, 107,  98,  88,  78,  36,  21])
        self.EJ301_LY = (800.-self.EJ301_pLY) * (3.5e3/(800.-91.))
        self.EJ301_Ee = self.EJ301_pEe * 2.5 / 653.
        words = 'LS_yield.__init__ EJ301 (Ee,LY) ='
        for Ee,LY in zip(self.EJ301_Ee,self.EJ301_LY):
            words += ' ({:.2f},{:.0f})'.format(Ee,LY)
        print(words)



        # 2205.15046 = Ref4 dict[solvent] = [WLS, (rLY,drLY), LY]
        # results from table 2
        self.LY4 = {'EJ301'     : [ 'NA', (78./100.,0.), 13570],
                    'p-Xylene'  : [ '4 g/L PPO', (77.2/100.,1.9/100.), 13430],
                    'LAB40'     : [ '4 g/L PPO', (60.2/100.,1.4/100.), 10470],
                    'DIN40'     : [ '4 g/L PPO', (70.0/100.,1.7/100.), 12190],
                    'relativeTo': ['anthracene']}

        # NIM A701 (2013) 133
        # results from tables 1, 6 with uncertainties taken from
        #  statements in the text
        self.LY3 = {'LAB15'  : ['1.5 g/L PPO', (1.16, 1.16*.03)],
                    'LAB20'  : ['2.0 g/L PPO', (1.22, 1.22*.03)],
                    'LAB30'  : ['3.0 g/L PPO', (1.26, 1.26*.03)],
                    'LAB60'  : ['6.0 g/L PPO', (1.33, 1.33*.03)],
                    'DIN15'  : ['1.5 g/L PPO', (0.96, 0.96*.03)],
                    'DIN20'  : ['2.0 g/L PPO', (0.99, 0.99*.03)],
                    'DIN30'  : ['3.0 g/L PPO', (1.07, 1.07*.03)],
                    'DIN60'  : ['6.0 g/L PPO', (1.12, 1.12*.03)],
                    'relativeTo': ['PC 1.5 g/L PPO']}

        # 2003.10491
        # values and uncertainies estimated from figure 6
        self.LY5 = {'LAB05'  : ['0.5 g/L PPO', (0.37, 0.01)],
                    'LAB10'  : ['1.0 g/L PPO', (0.46, 0.01)],
                    'LAB20'  : ['2.0 g/L PPO', (0.64, 0.01)],
                    'LAB30'  : ['3.0 g/L PPO', (0.67, 0.01)],
#                    'LAB40est':['4.0 g/L PPO', (0.70, 0.01)],
                    'LAB50'  : ['5.0 g/L PPO', (0.72, 0.01)],
                    'relativeTo': ['EJ301']}

        # 1812.02998 Production and Properties of the Liquid Scintillators used in the Stereo Reactor Neutrino Experiment
        # values taken from text and Table 3
        # Table 3 lists LY of 8700 and 9000 ph/MeV for LAB30 and LAB70, resp., but text states 8400 ph/MeV for LAB30.
        # Text states Anthracene produces about 17400 ph/MeV as reference
        # Ref7
        self.LY7 = {'LAB03' : ['3.0 g/L PPO',(8700./17400., 0.05*8700./17400.), 8700.],
                    'LAB07' : ['7.0 g/L PPO',(9000./17400., 0.05*9000./17400.), 9000.],
                    'relativeTo': ['anthracene']}


        ### 2011.12924 Development, characterisation, and deployment of the SNO+ liquid scintillator
        ### Ref8
        ### LAB20x = measured rLY 0.96 of LY of PC + 2 g/L PPO, inferred aLY 11900(1100)
        ### LAB20a = measured aLY 10830(570) based on data and simulation. Not de-oxygenated.
        ### LAB20b = measured aLY 11920(630) based on data and simulation. De-oxygenated.
        x,dx = 11900.,1100.
        ref = x/0.96
        a,da = 10830.,570.
        b,db = 11920.,630.
        self.LY8 = {'LAB20x': ['2.0 g/L PPO',(x/ref,dx/ref), x],
                    'LAB20a': ['2.0 g/L PPO',(a/ref,da/ref), a],
                    'LAB20b': ['2.0 g/L PPO',(b/ref,db/ref), b],
                    'relativeTo': ['PC 2.0 g/L PPO']}

        self.LY = {'Ref3' : self.LY3,
                   'Ref4' : self.LY4,
                   'Ref5' : self.LY5,
                   'Ref7' : self.LY7,
                   'Ref8' : self.LY8}

        
        
        return
    def plotLYvEe(self):
        '''
        plot Ref6 results and set Ref6 values
        '''
        refTitle = 'Measurement of the dependence of the LY of LAB-based\n and EJ-301 scintillators on electron energy'
        LABEe, LABLY = self.LAB_Ee, self.LAB_LY
        EJEe,  EJLY  = self.EJ301_Ee, self.EJ301_LY

        fig, axs = plt.subplots(2,1,sharex=True)
        
        axs[0].plot(LABEe,LABLY,'o',color='black',linestyle='dotted',label='LAB')
        axs[0].plot(EJEe ,EJLY ,'o',color='blue',linestyle='dotted',label='EJ301')
        axs[0].set_title('Ref6 ' + refTitle)
        axs[0].legend(loc='best')
        axs[0].grid()
        axs[0].set_ylabel('Light yield')
        #plt.show()

        delMax = 0.02
        Emin = 0.30
        cuts = 'Require $\Delta E <${:.2f}, Emin {:.2f} MeV'.format(delMax,Emin)
        rLY,rE,dE = [],[],[]
        for E,LY in zip(EJEe,EJLY):
            idx = numpy.argmin(numpy.abs(LABEe-E))
            iLY,iE = LABLY[idx],LABEe[idx]
            delta = abs(E-iE)/2.
            if delta<delMax and E>Emin and iE>Emin: 
                rLY.append( iLY/LY )
                rE.append( (E+iE)/2.)
                dE.append( delta)
        rLY,rE,dE = self.toNPA(rLY), self.toNPA(rE), self.toNPA(dE)
        print('dE',dE)
        axs[1].errorbar(rE,rLY,xerr=dE,marker='o',color='black',label='LAB/EJ301')

        m,s = numpy.mean(rLY),numpy.std(rLY)
        rmi = m - 3.*s
        rma = m + 3.*s
        words = 'mean {:.3f} std.dev. {:.3f}'.format(m,s)
        for i,ls in zip([-1.,0.,1.],['dotted','dashed','dotted']):
            axs[1].plot([min(rE),max(rE)],[m+i*s,m+i*s],linestyle=ls,color='grey')
        
        axs[1].set_ylabel('LY(LAB)/LY(EJ301)')
        axs[1].set_xlabel('Electron energy (MeV)')
        axs[1].grid()
        axs[1].set_ylim(rmi,rma)
        axs[1].sharex(axs[0])
        #axs[1].set_xscale('log')
        axs[1].text(min(rE)+.5,m-2.*s,words)
        axs[1].text(min(rE),m+2.*s,cuts)

        fig.subplots_adjust(hspace=0)

        png = self.figDir + 'Ref6_LAB_relative_to_EJ301.png'
        plt.savefig(png)
        print('LS_yield.plotLYvEe Wrote',png)

        plt.show()
        
        ### set Ref6 values
        self.LY6 = {self.Ref6_name : [self.Ref6_conc, (m,s)],
                    'relativeTo':[self.Ref6_relTo ]}
        self.LY['Ref6'] = self.LY6
        
        return
    def getConc(self,words):
        '''
        return concentration in percent by parsing input words
        '''
        s = words.split(' ')[0]
        if s=='NA' : return 0.
        return float(s)
        
    def unpackRef(self,ref,solvent='LAB',report=False):
        '''
        return relTo, X,Y,dY,Z  = reference LS, WLS conc, relative LY, uncertainty on relative light yield, absolute LY
        relTo is a string
        X, Y, dY, Z are numpy arrays. Z can be zero length. 
        given ref
        '''
        if ref not in self.LY :
            sys.exit('LS_yield.unpackRef ERROR Invalid ref '+ref)

        x,y,dy,z = [],[],[],[]
        relTo = None
        block = self.LY[ref]
        for LS in block:
            
            if LS=='relativeTo' :
                relTo = block[LS][0]
            elif solvent in LS:
                conc = block[LS][0]
                v,dv = block[LS][1]
                if len(block[LS])> 2: z.append( block[LS][2] )
                x.append( self.getConc(conc) )
                y.append( v )
                dy.append( dv )
        if relTo is None: sys.exit('LS_yield.unpackRef ERROR No relativeTo for ref '+ ref)
            
        X,Y,dY,Z = self.toNPA(x),self.toNPA(y),self.toNPA(dy),self.toNPA(z)
        
        if report : self.reportUnp(ref,solvent,relTo, X,Y,dY,Z)
        return relTo, X,Y,dY,Z 
    def toNPA(self,a):
        return numpy.array(a)
    def linear(self, x, m, b):
        return m*x + b
    def reportUnp(self,ref,solvent,relTo,X,Y,dY,Z):
        '''
        report results of unpacking
        '''
        words = ''
        if len(X)==0: words = 'NO ENTRIES'
        print('LS_yield.reportUnp ref',ref,'solvent',solvent,'relTo',relTo,words)
        if len(X)>0:
            print('LS_yield.reportUnp conc',X)
            print('LS_yield.reportUnp rLY ',Y)
            print('LS_yield.reportUnp drLY',dY)
        if len(Z)>0 : print('LS_yield.reportUnp aLY',Z)
        return
    def makeRatio(self,ref,sol1,sol2,report=False):
        '''
        for ref return ratio of response of solvent sol1 to solvent sol2
        same variables returned as unpackRef
        '''
        Num = relTo, X,Y,dY,Z = self.unpackRef(ref,solvent=sol1,report=report)
        Den = relTo, X,Y,dY,Z = self.unpackRef(ref,solvent=sol2,report=report)

        relTo = sol2
        x,y,dy,z = [],[],[],[]
        Xden,Yden,dYden = Den[1],Den[2],Den[3]
        Xnum,Ynum,dYnum = Num[1],Num[2],Num[3]
        
        for i,pair in enumerate(zip(Ynum,dYnum)):
            if len(Xden)>1 and abs(Xden[i]-Xnum[i])>1.e-6 :
                print('LS_yield.makeRatio ERROR i',i,'Xden',Xden[i],'Xnum[i]',Xnum[i])
                sys.exit('LS_yield.makeRatio ERROR concentrations not equal')
            v,dv = Yden[i],dYden[i]
            u,du = pair
            r = u/v
            dr= r*math.sqrt(du*du/u/u + dv*dv/v/v)
            y.append( r )
            dy.append(dr)
            x.append( Xnum[i] )
            
        X,Y,dY,Z = self.toNPA(x),self.toNPA(y),self.toNPA(dy),self.toNPA(z)

        ratio = sol1 
        if report : self.reportUnp(ref,ratio,relTo, X,Y,dY,Z)

        return relTo, X,Y,dY,Z        
    def makeName(self,a,b,c):
        s = a
        if len(b)>0 : s+= '_' + b
        if len(c)>0 : s+= '_' + c
        return s
    def main(self):

        ### plots Ref6 data and sets Ref6 relative LY
        self.plotLYvEe()
        
        
        # ref3 = LAB, DIN relative to PC + 1.5 g/L PPO
        # ref4 = LAB, DIN, EJ301 relative to anthracene
        # ref5 = LAB relative to EJ301
        # ref6 = LAB relative to EJ301
        # ref7 = LAB relative to anthracene
        # ref8 = LAB relative to PC + 2.0 g/L PPO

        # put unpacked results into dict C
        C = {}
        report = True
        for ref in ['Ref3','Ref4','Ref5','Ref6','Ref7','Ref8']:
            Solvents = ['EJ301','LAB','DIN']
            if ref=='Ref8' : Solvents = ['LAB20x', 'LAB20a', 'LAB20b']
            for solvent in Solvents:
                C[self.makeName(ref,solvent,'')] = relTo, X,Y,dY,Z = self.unpackRef(ref,solvent=solvent,report=report)

        ##### plot LAB relative to DIN vs concentration
        ref = 'Ref4'
        sol1,sol2 = 'LAB','DIN'
        label = self.makeName(ref,sol1,sol2)
        C[label] = self.makeRatio(ref,sol1,sol2,report=report)
        relToA,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='black',linestyle='dotted',label=label)

        ref = 'Ref3'
        sol1,sol2 = 'LAB','DIN'
        label = self.makeName(ref,sol1,sol2)
        C[label] = self.makeRatio(ref,sol1,sol2,report=report)
        relToB,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='blue',linestyle='dotted',label=label)

        if relToA!=relToB : sys.exit('LS_yield.main ERROR relToA '+relToA+' not equal to '+relToB)
        title = 'LAB light yield relative to '+relToA
        plt.title(title)
        plt.xlabel('PPO concentration in g/L')
        plt.ylabel('Relative light yield')
        plt.ylim(0.,1.3)
        plt.grid()
        plt.legend(loc='best')
        pdf = self.figDir + title.replace(' ','_')
        plt.savefig(pdf)
        print('LS_yield.main Wrote',pdf)
        plt.show()
        
        ##### plot LABPPO relative to EJ301
        ref = 'Ref4'
        sol1,sol2 = 'LAB','EJ301'
        C[self.makeName(ref,sol1,sol2)] = self.makeRatio(ref,sol1,sol2,report=report)
        
        label = self.makeName('Ref5','LAB','')
        relToA,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='black',linestyle='dotted',label=label)
        
        label = self.makeName('Ref6','LAB','')
        relToA,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='brown',linestyle='dotted',label=label)

        label = self.makeName('Ref4','LAB','EJ301')
        relToB,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='blue',linestyle='dotted',label=label)

        if relToA!=relToB : sys.exit('LS_yield.main ERROR relToA '+relToA+' not equal to '+relToB)

        title = 'LAB light yield relative to '+relToA
        plt.title(title)
        plt.xlabel('PPO concentration in g/L')
        plt.ylabel('Relative light yield')
        plt.ylim(0.19,.81)
        plt.grid()
        plt.legend(loc='best')
        pdf = self.figDir + title.replace(' ','_')
        plt.savefig(pdf)
        print('LS_yield.main Wrote',pdf)
        plt.show()

        ##### plot LABPPO relative to anthracene
        label = self.makeName('Ref4','LAB','')
        relToA,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='black',linestyle='dotted',label=label)
        
        label = self.makeName('Ref7','LAB','')
        relToB,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='brown',linestyle='dotted',label=label)
        
        if relToA!=relToB : sys.exit('LS_yield.main ERROR relToA '+relToA+' not equal to '+relToB)

        title = 'LAB light yield relative to '+relToA
        plt.title(title)
        plt.xlabel('PPO concentration in g/L')
        plt.ylabel('Relative light yield')
        plt.ylim(0.19,.81)
        plt.grid()
        plt.legend(loc='best')
        pdf = self.figDir + title.replace(' ','_')
        plt.savefig(pdf)
        print('LS_yield.main Wrote',pdf)
        plt.show()


        ##### plot LABPPO relative to PC
        label = self.makeName('Ref3','LAB','')
        relToA,X,Y,dY,Z = C[label]
        plt.errorbar(X,Y,yerr=dY,marker='o',color='black',linestyle='dotted',label=label)

        markers = ['o','v','s']
        colors  = ['brown','green','blue']
        for i,solvent in enumerate(['LAB20a','LAB20b','LAB20x']):
            offset = float(i-1)*.1
            label = self.makeName('Ref8',solvent,'')
            relToB,X,Y,dY,Z = C[label]
            plt.errorbar(X+offset,Y,yerr=dY,marker=markers[i],color=colors[i],linestyle='dotted',label=label)
        
        #if relToA!=relToB : sys.exit('LS_yield.main ERROR relToA '+relToA+' not equal to '+relToB)

        title = 'LAB light yield relative to '+relToA+'(Ref3) and '+relToB+'(Ref8)'
        
        plt.title(title)
        plt.xlabel('PPO concentration in g/L')
        plt.ylabel('Relative light yield')
        #plt.ylim(0.19,.81)
        plt.grid()
        plt.legend(loc='best')
        pdf = self.figDir + 'LAB LY relative to PC'.replace(' ','_')
        plt.savefig(pdf)
        print('LS_yield.main Wrote',pdf)
        plt.show()

        

        return

    def OLDmain(self):

        LY,dLY,x = [],[],[]
        for c in self.absLY:
            y,dy = self.absLY[c]
            LY.append(y)
            dLY.append(dy)
            x.append( self.conc[c] )
        LY, dLY, x = self.toNPA(LY), self.toNPA(dLY), self.toNPA(x)

        # BNL data
        BNLLY, dBNLLY, BNLx = [],[],[]
        for c in self.BNLabsLY:
            y,dy = self.BNLabsLY[c]
            BNLLY.append(y)
            dBNLLY.append(dy)
            BNLx.append( self.conc[c] )
        BNLLY, dBNLLY, BNLx = self.toNPA(BNLLY), self.toNPA(dBNLLY), self.toNPA(BNLx)

        ## do linear fit to the first 3 points
        param, cov = curve_fit(self.linear, x[:3], LY[:3], sigma=dLY[:3])
        print('LS_yield.main results of linear fit, param',param,'cov',cov)

        ## generate a family of fit parameters consistent with
        ## the fit results
        size = 500
        p = numpy.random.multivariate_normal(param,cov,size=size)

            
        
        X = numpy.linspace(min(x)/10.,max(x),1000)
        m, b = self.lfit['slope'][0],self.lfit['intercept'][0]
        Y = m*X+b

        for xma,rng in zip([101.,10.1],['Xfull','Xpartial']):
            for yscale in ['linear','log']:
                yw = 'Y'+yscale
                for xscale in ['linear','log']:
                    words = rng+yw+'X'+xscale 

                    ### plot data
                    plt.errorbar(x,LY,fmt='o',yerr=dLY,color='black')

                    ### plot BNL data
                    plt.errorbar(BNLx,BNLLY,fmt='o',yerr=dBNLLY,color='blue')

                    ### draw family of fit results as grey band
                    for pair in p:
                        m,b = pair
                        YY = m*X+b
                        plt.plot(X,YY,color='grey',linestyle='solid',alpha=0.1)
                    ### main fit result
                    plt.plot(X,Y,'r-')

                    plt.title('absolute LY arXiv:2006.00173 table 1')
                    plt.grid()

                    plt.xlabel('% Concentration')
                    plt.ylabel('Light yield (ph/MeV/% conc)')
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    if xscale=='linear': plt.xlim(-1.,xma)
                    if yscale=='linear' and xma<100. : plt.ylim(0.,2000.) 
                    png = 'LY_YIELD_FIGURES/LS_yield_arXiv2006.00173_'+words+'.png'
                    plt.savefig(png)
                    print('LS_yield.main Wrote',png)
                    plt.show()

        
        return
if __name__ == '__main__' :
    P = LS_yield()
    P.main()
    
