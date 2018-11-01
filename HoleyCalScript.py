import math
import random
from scipy.optimize import least_squares
import kinematics
import HoleyCalibration

cal=HoleyCalibration.HoleyCalibration()
cal.kin.isQuadKinematics=False

# Starting machine calibration
cal.SP_D=3601.2
cal.SP_motorOffsetY=468.4
cal.SP_rotationDiskRadius=139.1
cal.SP_chainSagCorrection=31.865887
cal.SP_leftChainTolerance=0
cal.SP_rightChainTolerance=0
cal.SP_chainOverSprocket=1

# adjust based upon machine settings
workspaceHeight = 1219.2
workspaceWidth = 2438.4
gearTeeth = 10
chainPitch = 6.35

# Ideal x,y coordinates
aH0x = 0
aH0y = 0
aH1x = (workspaceWidth/2.0-254.0)*-1.0
aH1y = (workspaceHeight/2.0-254.0)*-1.0
aH2x = aH1x
aH2y = aH1y*-1.0
aH3x = aH1x*-1.0
aH3y = aH2y
aH4x = aH3x
aH4y = aH1y
IdealLengthArray=cal.SetIdealXyCoordinates(aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y)

#measured distances of hole pattern
##---CHANGE THESE TO WHAT YOU MEASURED---##
##---USE MILLIMETERS ONLY---##
##---My tape measure was off by 101 mm so the -101.0 adjust for it---##
##---CHANGE IT BECAUSE YOURS IS LIKELY DIFFERENT---###
dH0H1 = 1133.0-101.0
dH0H2 = 1133.0-101.0
dH0H3 = 1129.0-101.0
dH0H4 = 1133.0-101.0
dH1H2 = 814.0-101.0
dH1H4 = 2037.0-101.0
dH2H3 = 2036.0-101.0
dH3H4 = 812.0-101.0
dH0M5 = 466.0-101.0
dH2M5 = 1070.0-101.0
dH2Top=254
dH3Top=255
cal.SetMeasurements(dH0H1,dH0H2,dH0H3,dH0H4,dH1H2,dH1H4,dH2H3,dH3H4,dH2Top,dH3Top)
       
cal.Calibrate()

CalResults=cal.OptimizationOutput.x

CalibratedError=cal.LengthDeltaFromIdeal(CalResults)

print("Calibrated Error:")
print(CalibratedError)
print("")
print("Calibration Deltas:")
print(CalResults)
print("")
print("Distance Between Motors:")
print(cal.Opt_D)
print("")
print("Motor Y Offset:")
print(cal.Opt_motorOffsetY)
print("")
print("Rotational Radius:")
print(cal.Opt_rotationDiskRadius)
print("")
print("Chain Sag Correction Constant:")
print(cal.Opt_chainSagCorrection)
print("")
print("Left Chain Tolerance:")
print(cal.Opt_leftChainTolerance)
print("")
print("Right Chain Tolerance:")
print(cal.Opt_rightChainTolerance)
print("")
print("Errors:")
print(cal.OptimizationOutput.fun)