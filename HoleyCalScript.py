import math
import random
from scipy.optimize import least_squares
import kinematics
import HoleyCalibration
def GeometricLength(x1,y1,x2,y2):
    return math.sqrt(math.pow(x1-x2,2) + math.pow(y1-y2,2))

cal=HoleyCalibration.HoleyCalibration()
cal.kin.isQuadKinematics=False

# Starting machine calibration
cal.SP_D=3048
cal.SP_motorOffsetY=711
cal.SP_rotationDiskRadius=139.1
cal.SP_sledWeight=20
cal.SP_leftChainTolerance=0
cal.SP_rightChainTolerance=0
cal.SP_chainOverSprocket=0

# adjust based upon machine settings
workspaceHeight = 1219.2
workspaceWidth = 2438.4
gearTeeth = 10
chainPitch = 6.35

# Ideal x,y coordinates
aH1x = -(workspaceWidth/2.0-254.0)
aH1y = (workspaceHeight/2.0-254.0)
aH2x = 0
aH2y = aH1y
aH3x = -aH1x
aH3y = aH1y
aH4x = aH1x
aH4y = -aH1y
aH5x = 0
aH5y = aH4y
aH6x = aH3x
aH6y = aH4y
IdealLengthArray=cal.SetIdealXyCoordinates(aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y,aH5x,aH5y,aH6x,aH6y)

#measured distances of hole pattern
##---CHANGE THESE TO WHAT YOU MEASURED---##
##---USE MILLIMETERS ONLY---##
##---My tape measure was off by 101 mm so the -101.0 adjust for it---##
##---CHANGE IT BECAUSE YOURS IS LIKELY DIFFERENT---###

#dH0H1 = 1133.0-101.0
#dH0H2 = 1133.0-101.0
#dH0H3 = 1129.0-101.0
#dH0H4 = 1133.0-101.0
#dH1H2 = 814.0-101.0
#dH1H4 = 2037.0-101.0
#dH2H3 = 2036.0-101.0
#dH3H4 = 812.0-101.0
#dH0M5 = 466.0-101.0
#dH2M5 = 1070.0-101.0
#dH2Top=254
#dH3Top=255

#Simulate Measurement.  Modify machine parameters, use kin.forward to determine x,y coordinates
cal.SP_D=3040
cal.SP_motorOffsetY=700
cal.SP_rotationDiskRadius=139.1
cal.SP_sledWeight=20
cal.SP_leftChainTolerance=.5
cal.SP_rightChainTolerance=-.5
cal.kin.D=cal.SP_D
cal.kin.motorOffsetY=cal.SP_motorOffsetY
cal.kin.rotationDiskRadius=cal.SP_rotationDiskRadius
cal.kin.sledWeight=cal.SP_sledWeight
cal.kin.leftChainTolerance=cal.SP_leftChainTolerance
cal.kin.rightChainTolerance=cal.SP_rightChainTolerance
cal.kin.recomputeGeometry()

mH1x,mH1y=cal.kin.forward(cal.LC01,cal.RC01)
mH2x,mH2y=cal.kin.forward(cal.LC02,cal.RC02)
mH3x,mH3y=cal.kin.forward(cal.LC03,cal.RC03)
mH4x,mH4y=cal.kin.forward(cal.LC04,cal.RC04)
mH5x,mH5y=cal.kin.forward(cal.LC05,cal.RC05)
mH6x,mH6y=cal.kin.forward(cal.LC06,cal.RC06)

dH1H2=GeometricLength(mH1x,mH1y,mH2x,mH2y)
dH2H3=GeometricLength(mH2x,mH2y,mH3x,mH3y)
dH4H5=GeometricLength(mH4x,mH4y,mH5x,mH5y)
dH5H6=GeometricLength(mH5x,mH5y,mH6x,mH6y)
dH1H4=GeometricLength(mH1x,mH1y,mH4x,mH4y)
dH2H5=GeometricLength(mH2x,mH2y,mH5x,mH5y)
dH3H6=GeometricLength(mH3x,mH3y,mH6x,mH6y)
dH2H4=GeometricLength(mH2x,mH2y,mH4x,mH4y)
dH1H5=GeometricLength(mH1x,mH1y,mH5x,mH5y)
dH3H5=GeometricLength(mH3x,mH3y,mH5x,mH5y)
dH2H6=GeometricLength(mH2x,mH2y,mH6x,mH6y)
dH2Top=GeometricLength(mH2x,mH2y,mH2x,workspaceHeight/2)

#Reset the IdealLengthArray because this is what the machine thought the cal was when it ran the test-holes
IdealLengthArray=cal.SetIdealXyCoordinates(aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y,aH5x,aH5y,aH6x,aH6y)

cal.SetMeasurements(dH1H2,dH2H3,dH4H5,dH5H6,dH1H4,dH2H5,dH3H6,dH2H4,dH1H5,dH3H5,dH2H6,dH2Top)
                  # [dH1H2,dH2H3,dH4H5,dH5H6,dH1H4,dH2H5,dH3H6,dH2H4,dH1H5,dH3H5,dH2H6,dH2Top]
                  #  dH1H2,dH2H3,dH4H5,dH5H6,dH1H4,dH2H5,dH3H6,dH2H4,dH1H5,dH3H5,dH2H6,dH2Top

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
#print("Rotational Radius:")
#print(cal.Opt_rotationDiskRadius)
#print("")
#print("Sled Weight:")
#print(cal.Opt_sledWeight)
#print("")
print("Left Chain Tolerance:")
print(cal.Opt_leftChainTolerance)
print("")
print("Right Chain Tolerance:")
print(cal.Opt_rightChainTolerance)
print("")
print("Errors:")
print(cal.OptimizationOutput.fun)