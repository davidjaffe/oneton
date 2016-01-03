# -*- coding: utf-8 -*-
import h5py
import numpy as np

f = h5py.File('E:\\1TonData\\SignalPMT_HV\\S0-S4_2000V_Bright_150918.h5', 'r')
CalData1 = f['Calibration']['Digitizer_1']['Pedestal']
CalData2 = f['Calibration']['Digitizer_2']['Pedestal']

EvtDir = f['Events']
Dig1Evts = []
cbaseline = []
cdummy = []
tdummy = []
PulseCharge = [[], [], [], []]
PulseTime = [[],[],[],[]]
PulseWidth = [[],[],[],[]]
for evtkey in EvtDir.keys():
    #Subtract BG from digitizers, do some analysis, save a ROOT file.
    thisdig1 = EvtDir[evtkey]['Digitizer_1']
    theDig1 = []
    theDig1 = np.subtract(thisdig1, CalData1)
    Dig1Evts.append(theDig1)
    cbaseline.append([sum((theDig1[0][160:170], theDig1[0][240:250])), sum((theDig1[1][160:170], theDig1[1][240:250])), sum((theDig1[2][160:170], theDig1[2][240:250])), sum((theDig1[3][160:170], theDig1[3][240:250]))])
    cdummy.append([sum(theDig1[0][170:240]), sum(theDig1[1][170:240]), sum(theDig1[2][170:240]), sum(theDig1[3][170:240])])
    tdummy.append([np.argmin(theDig1[0][170:240]), np.argmin(theDig1[1][170:240]), np.argmin(theDig1[2][170:240]), np.argmin(theDig1[3][170:240])])
    #Step along wfm and analyse for peaks
    threshold = -10
    TailLength = 10
    #print("Event = ", evtkey)
    for i in np.linspace(0, 3, 4):
        chan = int(i)
        MvAvgPrev = 0
        ispulse = 0
        BLstart = 0
        BLend = 0
        Qpulse = 0
        Wpulse = 0
        NumAfter = 0
        BL = 0
        theTime = 0
        #print("Channel = ", chan)
        for j in np.linspace(41, 2558, 2518):
            #if np.mod(int(j), 500)==0:
            #    print(j)
            if ((theDig1[chan][int(j)]-MvAvgPrev)<threshold)&(not(ispulse)):
                #process a new pulse
                #print("Found a pulse!")
                ispulse=1
                BLstart = MvAvgPrev
                Qpulse += theDig1[chan][int(j)]
                Wpulse += 1
                theTime = int(j)
            elif ispulse:
                #still in a pulse
                if (theDig1[chan][int(j)]-MvAvgPrev)<threshold:
                    Qpulse += theDig1[chan][int(j)]
                    Wpulse += 1
                    NumAfter = 0
                else:
                    Qpulse += theDig1[chan][int(j)]
                    Wpulse += 1
                    NumAfter += 1
                    if NumAfter == TailLength:
                        #process the previous pulse
                        ispulse == 0
                        #print("Pulse Ended!!")
                        BLend = (theDig1[chan][int(j-1)] + theDig1[chan][int(j)] + theDig1[chan][int(j+1)])/3#np.mean(theDig1[chan][int(j-1):int(j+1)])
                        BL = np.mean((BLstart, BLend))
                        Qpulse -= BL*Wpulse
                        #Add pulse charge to dictionary/array
                        PulseCharge[chan] = np.append(PulseCharge[chan], Qpulse)
                        PulseTime[chan] = np.append(PulseTime[chan], theTime)
                        PulseWidth[chan] = np.append(PulseWidth[chan], Wpulse)
                        #print("Charge = ", Qpulse, " Time = ", theTime, " Width = ", Wpulse, "Chan = ", chan)
                        NumAfter = 0
                        ispulse = 0
                        Wpulse = 0
                        MvAvgPrev = BLend
            else:
                MvAvgPrev = (theDig1[chan][int(j-1)] + theDig1[chan][int(j)] + theDig1[chan][int(j+1)])/3
                
        
netcharge = np.subtract(cdummy, np.multiply(np.mean(cbaseline),70/20))
Charge = np.transpose(cdummy)
MinTime = np.transpose(tdummy)

f.close()
