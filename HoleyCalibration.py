import kinematics
import math
from scipy.optimize import curve_fit
import numpy
class HoleyCalibration():
    SP_D=3601.2
    SP_motorOffsetY=468.4
    SP_rotationDiskRadius=139.1
    SP_chainSagCorrection=31.865887
    
    Opt_D=3601.2
    Opt_motorOffsetY=468.4
    Opt_rotationDiskRadius=139.1
    Opt_chainSagCorrection=31.865887
    
    #Chain lengths @ each hole
    LC00=0
    RC00=0
    LC01=0
    RC01=0
    LC02=0
    RC02=0
    LC03=0
    RC03=0
    LC04=0
    RC04=0
    
    IdealLengthArray=0
    MeasuredLengthArray=0
    DesiredLengthDeltaArray=0
    
    OptimalMachineParameterDeltas=0
    CoefficientOfVariance=0
    
    kin=kinematics.Kinematics()
    #Define function with input of (ideal lengths and) machine parameters (delta) and output of length error
    def LengthChangeFromStartingPoint(self,MeasuredLengthArray,Del_D,Del_motorOffsetY,Del_rotationDiskRadius,Del_chainSagCorrection):
        self.kin.D=self.SP_D+Del_D
        self.kin.motorOffsetY=self.SP_motorOffsetY+Del_motorOffsetY
        self.kin.rotationDiskRadius=self.SP_rotationDiskRadius+Del_rotationDiskRadius
        self.kin.chainSagCorrection=self.SP_chainSagCorrection+Del_chainSagCorrection
        
        aH0x,aH0y=self.kin.forward(self.LC00,self.RC00)
        aH1x,aH1y=self.kin.forward(self.LC01,self.RC01)
        aH2x,aH2y=self.kin.forward(self.LC02,self.RC02)
        aH3x,aH3y=self.kin.forward(self.LC03,self.RC03)
        aH4x,aH4y=self.kin.forward(self.LC04,self.RC04)
        
        return self.CalculateLengthArray(aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y)-MeasuredLengthArray
    def LengthChangePlusOne(self,MeasuredLengthArray,Del_D,Del_motorOffsetY,Del_rotationDiskRadius,Del_chainSagCorrection):
        LenChange=self.LengthChangeFromStartingPoint(MeasuredLengthArray,Del_D,Del_motorOffsetY,Del_rotationDiskRadius,Del_chainSagCorrection)
        return LenChange+1
    
    def CalculateLengthArray(self,aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y):
        LenH00H01=GeometricLength(aH0x,aH0y,aH1x,aH1y)
        LenH00H02=GeometricLength(aH0x,aH0y,aH2x,aH2y)
        LenH00H03=GeometricLength(aH0x,aH0y,aH3x,aH3y)
        LenH00H04=GeometricLength(aH0x,aH0y,aH4x,aH4y)
        LenH01H02=GeometricLength(aH1x,aH1y,aH2x,aH2y)
        LenH01H04=GeometricLength(aH1x,aH1y,aH4x,aH4y)
        LenH02H03=GeometricLength(aH2x,aH2y,aH3x,aH3y)
        LenH03H04=GeometricLength(aH3x,aH3y,aH4x,aH4y)
        return numpy.array([LenH00H01,LenH00H02,LenH00H03,LenH00H04,LenH01H02,LenH01H04,LenH02H03,LenH03H04])
    
    def SetIdealXyCoordinates(self,aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y):
        self.kin.D=self.SP_D
        self.kin.motorOffsetY=self.SP_motorOffsetY
        self.kin.rotationDiskRadius=self.SP_rotationDiskRadius
        self.kin.chainSagCorrection=self.SP_chainSagCorrection
        self.IdealLengthArray=self.CalculateLengthArray(aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y)
        
        self.LC00,self.RC00=self.kin.inverse(aH0x,aH0y)
        self.LC01,self.RC01=self.kin.inverse(aH1x,aH1y)
        self.LC02,self.RC02=self.kin.inverse(aH2x,aH2y)
        self.LC03,self.RC03=self.kin.inverse(aH3x,aH3y)
        self.LC04,self.RC04=self.kin.inverse(aH4x,aH4y)
        return self.IdealLengthArray    
    
    def SetMeasurements(self,dH0H1,dH0H2,dH0H3,dH0H4,dH1H2,dH1H4,dH2H3,dH3H4):
        self.MeasuredLengthArray=numpy.array([dH0H1,dH0H2,dH0H3,dH0H4,dH1H2,dH1H4,dH2H3,dH3H4])
    
    def Calibrate(self):
        TargetLengthChange=self.IdealLengthArray-self.MeasuredLengthArray+1
        self.OptimalMachineParameterDeltas,self.CoefficientOfVariance=curve_fit(self.LengthChangePlusOne,self.MeasuredLengthArray,TargetLengthChange,p0=numpy.array([1,1,1,1]),absolute_sigma=True,sigma=numpy.array([.5,.5,.5,.5,.5,.5,.5,.5]))
        self.Opt_D=self.OptimalMachineParameterDeltas[0]+self.SP_D
        self.Opt_motorOffsetY=self.OptimalMachineParameterDeltas[1]+self.SP_motorOffsetY
        self.Opt_rotationDiskRadius=self.OptimalMachineParameterDeltas[2]+self.SP_rotationDiskRadius
        self.Opt_chainSagCorrection=self.OptimalMachineParameterDeltas[3]+self.SP_chainSagCorrection
        self.kin.D=self.Opt_D
        self.kin.motorOffsetY=self.Opt_motorOffsetY
        self.kin.rotationDiskRadius=self.Opt_rotationDiskRadius
        self.kin.chainSagCorrection=self.Opt_chainSagCorrection

def GeometricLength(x1,y1,x2,y2):
    return math.sqrt(math.pow(x1-x2,2) + math.pow(y1-y2,2))
    