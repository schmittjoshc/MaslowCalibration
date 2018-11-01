import kinematics
import math
from scipy.optimize import least_squares
import numpy
class HoleyCalibration():
    SP_D=3601.2
    SP_motorOffsetY=468.4
    SP_rotationDiskRadius=139.1
    SP_chainSagCorrection=31.865887
    SP_leftChainTolerance=0
    SP_rightChainTolerance=0
    
    Opt_D=3601.2
    Opt_motorOffsetY=468.4
    Opt_rotationDiskRadius=139.1
    Opt_chainSagCorrection=31.865887
    Opt_leftChainTolerance=0
    Opt_rightChainTolerance=0
    
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
    
    OptimizationOutput=0
    CoefficientOfVariance=0
    
    kin=kinematics.Kinematics()
    #Define function with input of (ideal lengths and) machine parameters (delta) and output of length error
    def LengthDeltaFromIdeal(self,DeltaArray): #Del_D,Del_motorOffsetY,Del_rotationDiskRadius,Del_chainSagCorrection):
        self.kin.D=self.SP_D+DeltaArray[0]
        self.kin.motorOffsetY=self.SP_motorOffsetY+DeltaArray[1]
        self.kin.rotationDiskRadius=self.SP_rotationDiskRadius+DeltaArray[2]
        self.kin.chainSagCorrection=self.SP_chainSagCorrection+DeltaArray[3]
        self.kin.leftChainTolerance=self.SP_leftChainTolerance+DeltaArray[4]
        self.kin.rightChainTolerance=self.SP_rightChainTolerance+DeltaArray[5]
        self.kin.recomputeGeometry()
        aH0x,aH0y=self.kin.forward(self.LC00,self.RC00)
        aH1x,aH1y=self.kin.forward(self.LC01,self.RC01)
        aH2x,aH2y=self.kin.forward(self.LC02,self.RC02)
        aH3x,aH3y=self.kin.forward(self.LC03,self.RC03)
        aH4x,aH4y=self.kin.forward(self.LC04,self.RC04)
        
        return self.CalculateLengthArray(aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y)+self.MeasuredLengthArray-self.IdealLengthArray-self.IdealLengthArray
    
    def CalculateLengthArray(self,aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y):
        LenH00H01=GeometricLength(aH0x,aH0y,aH1x,aH1y)
        LenH00H02=GeometricLength(aH0x,aH0y,aH2x,aH2y)
        LenH00H03=GeometricLength(aH0x,aH0y,aH3x,aH3y)
        LenH00H04=GeometricLength(aH0x,aH0y,aH4x,aH4y)
        LenH01H02=GeometricLength(aH1x,aH1y,aH2x,aH2y)
        LenH01H04=GeometricLength(aH1x,aH1y,aH4x,aH4y)
        LenH02H03=GeometricLength(aH2x,aH2y,aH3x,aH3y)
        LenH03H04=GeometricLength(aH3x,aH3y,aH4x,aH4y)
        LenH01Top=self.kin.machineHeight/2-aH2y
        LenH04Top=self.kin.machineHeight/2-aH3y
        return numpy.array([LenH00H01,LenH00H02,LenH00H03,LenH00H04,LenH01H02,LenH01H04,LenH02H03,LenH03H04,LenH01Top,LenH04Top])
    
    def SetIdealXyCoordinates(self,aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y):
        self.kin.D=self.SP_D
        self.kin.motorOffsetY=self.SP_motorOffsetY
        self.kin.rotationDiskRadius=self.SP_rotationDiskRadius
        self.kin.chainSagCorrection=self.SP_chainSagCorrection
        self.kin.leftChainTolerance=self.SP_leftChainTolerance
        self.kin.rightChainTolerance=self.SP_rightChainTolerance
        self.kin.recomputeGeometry()
        self.IdealLengthArray=self.CalculateLengthArray(aH0x,aH0y,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y)
        
        self.LC00,self.RC00=self.kin.inverse(aH0x,aH0y)
        self.LC01,self.RC01=self.kin.inverse(aH1x,aH1y)
        self.LC02,self.RC02=self.kin.inverse(aH2x,aH2y)
        self.LC03,self.RC03=self.kin.inverse(aH3x,aH3y)
        self.LC04,self.RC04=self.kin.inverse(aH4x,aH4y)
        return self.IdealLengthArray    
    
    def SetMeasurements(self,dH0H1,dH0H2,dH0H3,dH0H4,dH1H2,dH1H4,dH2H3,dH3H4,d1top,d4top):
        self.MeasuredLengthArray=numpy.array([dH0H1,dH0H2,dH0H3,dH0H4,dH1H2,dH1H4,dH2H3,dH3H4,d1top,d4top])
    
    def Calibrate(self):
        TargetLengthChange=self.IdealLengthArray-self.MeasuredLengthArray+1
        self.OptimizationOutput=least_squares(self.LengthDeltaFromIdeal,numpy.array([20,100,-10,100,0,0]),jac='2-point',diff_step=.1)
        Deltas=self.OptimizationOutput.x
        self.Opt_D=Deltas[0]+self.SP_D
        self.Opt_motorOffsetY=Deltas[1]+self.SP_motorOffsetY
        self.Opt_rotationDiskRadius=Deltas[2]+self.SP_rotationDiskRadius
        self.Opt_chainSagCorrection=Deltas[3]+self.SP_chainSagCorrection
        self.Opt_leftChainTolerance=Deltas[4]+self.SP_leftChainTolerance
        self.Opt_rightChainTolerance=Deltas[5]+self.SP_rightChainTolerance
        self.kin.D=self.Opt_D
        self.kin.motorOffsetY=self.Opt_motorOffsetY
        self.kin.rotationDiskRadius=self.Opt_rotationDiskRadius
        self.kin.chainSagCorrection=self.Opt_chainSagCorrection
        self.kin.leftChainTolerance=self.Opt_leftChainTolerance
        self.kin.rightChainTolerance=self.Opt_rightChainTolerance
        self.kin.recomputeGeometry()

def GeometricLength(x1,y1,x2,y2):
    return math.sqrt(math.pow(x1-x2,2) + math.pow(y1-y2,2))
    