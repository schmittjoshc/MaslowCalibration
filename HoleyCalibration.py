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
    SP_sledWeight=22
    
    Opt_D=3601.2
    Opt_motorOffsetY=468.4
    Opt_rotationDiskRadius=139.1
    Opt_chainSagCorrection=31.865887
    Opt_leftChainTolerance=0
    Opt_rightChainTolerance=0
    Opt_sledWeight=22
    
    #Chain lengths @ each hole
    LC01=0
    RC01=0
    LC02=0
    RC02=0
    LC03=0
    RC03=0
    LC04=0
    RC04=0
    LC05=0
    RC05=0
    LC06=0
    RC06=0
    
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
        #self.kin.rotationDiskRadius=self.SP_rotationDiskRadius+DeltaArray[2]
        #self.kin.sledWeight=abs(self.SP_sledWeight+DeltaArray[3])
        self.kin.leftChainTolerance=self.SP_leftChainTolerance+DeltaArray[4]
        self.kin.rightChainTolerance=self.SP_rightChainTolerance+DeltaArray[5]
        self.kin.recomputeGeometry()
        aH1x,aH1y=self.kin.forward(self.LC01,self.RC01)
        aH2x,aH2y=self.kin.forward(self.LC02,self.RC02)
        aH3x,aH3y=self.kin.forward(self.LC03,self.RC03)
        aH4x,aH4y=self.kin.forward(self.LC04,self.RC04)
        aH5x,aH5y=self.kin.forward(self.LC05,self.RC05)
        aH6x,aH6y=self.kin.forward(self.LC06,self.RC06)
        
        return (self.CalculateLengthArray(aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y,aH5x,aH5y,aH6x,aH6y)+self.MeasuredLengthArray-self.IdealLengthArray-self.IdealLengthArray)
    
    def CalculateLengthArray(self,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y,aH5x,aH5y,aH6x,aH6y):
        dH1H2=GeometricLength(aH1x,aH1y,aH2x,aH2y)
        dH2H3=GeometricLength(aH2x,aH2y,aH3x,aH3y)
        dH4H5=GeometricLength(aH4x,aH4y,aH5x,aH5y)
        dH5H6=GeometricLength(aH5x,aH5y,aH6x,aH6y)
        dH1H4=GeometricLength(aH1x,aH1y,aH4x,aH4y)
        dH2H5=GeometricLength(aH2x,aH2y,aH5x,aH5y)
        dH3H6=GeometricLength(aH3x,aH3y,aH6x,aH6y)
        dH2H4=GeometricLength(aH2x,aH2y,aH4x,aH4y)
        dH1H5=GeometricLength(aH1x,aH1y,aH5x,aH5y)
        dH3H5=GeometricLength(aH3x,aH3y,aH5x,aH5y)
        dH2H6=GeometricLength(aH2x,aH2y,aH6x,aH6y)
        dH2Top=GeometricLength(aH2x,aH2y,aH2x,self.kin.machineHeight/2)
        return numpy.array([dH1H2,dH2H3,dH4H5,dH5H6,dH1H4,dH2H5,dH3H6,dH2H4,dH1H5,dH3H5,dH2H6,dH2Top])
    
    def SetIdealXyCoordinates(self,aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y,aH5x,aH5y,aH6x,aH6y):
        self.kin.D=self.SP_D
        self.kin.motorOffsetY=self.SP_motorOffsetY
        self.kin.rotationDiskRadius=self.SP_rotationDiskRadius
        self.kin.sledWeight=self.SP_sledWeight
        self.kin.leftChainTolerance=self.SP_leftChainTolerance
        self.kin.rightChainTolerance=self.SP_rightChainTolerance
        self.kin.chainOverSprocket=self.SP_chainOverSprocket
        self.kin.recomputeGeometry()
        self.IdealLengthArray=self.CalculateLengthArray(aH1x,aH1y,aH2x,aH2y,aH3x,aH3y,aH4x,aH4y,aH5x,aH5y,aH6x,aH6y)
        
        #DeleteMe=self.kin.inverse(1219.2,609.6)
        self.LC01,self.RC01=self.kin.inverse(aH1x,aH1y)
        self.LC02,self.RC02=self.kin.inverse(aH2x,aH2y)
        self.LC03,self.RC03=self.kin.inverse(aH3x,aH3y)
        self.LC04,self.RC04=self.kin.inverse(aH4x,aH4y)
        self.LC05,self.RC05=self.kin.inverse(aH5x,aH5y)
        self.LC06,self.RC06=self.kin.inverse(aH6x,aH6y)
        return self.IdealLengthArray
    
    def SetMeasurements(self,dH1H2,dH2H3,dH4H5,dH5H6,dH1H4,dH2H5,dH3H6,dH2H4,dH1H5,dH3H5,dH2H6,dH2Top):
        self.MeasuredLengthArray=numpy.array([dH1H2,dH2H3,dH4H5,dH5H6,dH1H4,dH2H5,dH3H6,dH2H4,dH1H5,dH3H5,dH2H6,dH2Top])
    
    def Calibrate(self):
        TargetLengthChange=self.IdealLengthArray-self.MeasuredLengthArray+1
        self.OptimizationOutput=least_squares(self.LengthDeltaFromIdeal,numpy.array([0,0,0,0,0,0]),jac='2-point',diff_step=.1,ftol=1e-11)
        Deltas=self.OptimizationOutput.x
        self.Opt_D=Deltas[0]+self.SP_D
        self.Opt_motorOffsetY=Deltas[1]+self.SP_motorOffsetY
        #self.Opt_rotationDiskRadius=Deltas[2]+self.SP_rotationDiskRadius
        #self.Opt_sledWeight=Deltas[3]+self.SP_sledWeight
        self.Opt_leftChainTolerance=Deltas[4]+self.SP_leftChainTolerance
        self.Opt_rightChainTolerance=Deltas[5]+self.SP_rightChainTolerance
        self.kin.D=self.Opt_D
        self.kin.motorOffsetY=self.Opt_motorOffsetY
        #self.kin.rotationDiskRadius=self.Opt_rotationDiskRadius
        #self.kin.sledWeight=self.Opt_sledWeight
        self.kin.leftChainTolerance=self.Opt_leftChainTolerance
        self.kin.rightChainTolerance=self.Opt_rightChainTolerance
        self.kin.recomputeGeometry()

def GeometricLength(x1,y1,x2,y2):
    return math.sqrt(math.pow(x1-x2,2) + math.pow(y1-y2,2))
    