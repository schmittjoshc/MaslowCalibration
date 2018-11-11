#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 15:19:22 2018

@author: john
"""

import kinematics
import numpy
import matplotlib.pyplot as plt
from matplotlib import ticker
import math

kin=kinematics.Kinematics()
# Starting machine calibration
kin.D=3048
kin.motorOffsetY=711
kin.rotationDiskRadius=139.1
kin.sledWeight=20
kin.leftChainTolerance=0
kin.rightChainTolerance=0
kin.chainOverSprocket=0
kin.isQuadKinematics=False

workspaceHeight = 1219.2
workspaceWidth = 2438.4
NumPoints=30
xPosition=numpy.linspace(-workspaceWidth/2+10,workspaceWidth/2-10,NumPoints)
yPosition=numpy.linspace(-workspaceHeight/2+10,workspaceHeight/2-10,NumPoints)

leftChainStraight=numpy.zeros((NumPoints,NumPoints))
rightChainStraight=numpy.zeros((NumPoints,NumPoints))
leftChainCatenary=numpy.zeros((NumPoints,NumPoints))
rightChainCatenary=numpy.zeros((NumPoints,NumPoints))
xPositionCatenary=numpy.zeros((NumPoints,NumPoints))
yPositionCatenary=numpy.zeros((NumPoints,NumPoints))
xPositionElasticity=numpy.zeros((NumPoints,NumPoints))
yPositionElasticity=numpy.zeros((NumPoints,NumPoints))
xPositionTotal=numpy.zeros((NumPoints,NumPoints))
yPositionTotal=numpy.zeros((NumPoints,NumPoints))
CatenaryDelta=numpy.zeros((NumPoints,NumPoints))
ElasticityDelta=numpy.zeros((NumPoints,NumPoints))
TotalDelta=numpy.zeros((NumPoints,NumPoints))

def GeometricLength(x1,y1,x2,y2):
    return math.sqrt(math.pow(x1-x2,2) + math.pow(y1-y2,2))

for i in range(0,NumPoints):
    for j in range(0,NumPoints):
        kin.chainWeight=1E-6
        leftChainStraight[j,i],rightChainStraight[j,i]=kin.inverse(xPosition[i],yPosition[j])
        kin.chainWeight=.09/304.8
        xPositionCatenary[j,i],yPositionCatenary[j,i]=kin.forward(leftChainStraight[j,i],rightChainStraight[j,i])
        CatenaryDelta[j,i]=GeometricLength(xPosition[i],yPosition[j],xPositionCatenary[j,i],yPositionCatenary[j,i])
        kin.chainElasticity=1E-7
        leftChainStraight[j,i],rightChainStraight[j,i]=kin.inverse(xPosition[i],yPosition[j])
        kin.chainElasticity=.000023
        xPositionElasticity[j,i],yPositionElasticity[j,i]=kin.forward(leftChainStraight[j,i],rightChainStraight[j,i])
        ElasticityDelta[j,i]=GeometricLength(xPosition[i],yPosition[j],xPositionElasticity[j,i],yPositionElasticity[j,i])
        kin.chainWeight=1E-6
        kin.chainElasticity=1E-7
        leftChainStraight[j,i],rightChainStraight[j,i]=kin.inverse(xPosition[i],yPosition[j])
        kin.chainWeight=.09/304.8
        kin.chainElasticity=.000023
        xPositionTotal[j,i],yPositionTotal[j,i]=kin.forward(leftChainStraight[j,i],rightChainStraight[j,i])
        TotalDelta[j,i]=GeometricLength(xPosition[i],yPosition[j],xPositionTotal[j,i],yPositionTotal[j,i])
fig,ax=plt.subplots()
levels=numpy.linspace(numpy.min(CatenaryDelta),numpy.max(CatenaryDelta),10)
locator = ticker.LinearLocator()
plot=ax.contourf(xPosition,yPosition,CatenaryDelta,levels=levels)
ax.set_xlabel("X-Position [mm]")
ax.set_ylabel("Y-Position [mm]")
ax.set_title("Chain-Sag Displacement vs. Position on Work Surface")
clb  = plt.colorbar(plot, format='%.e', ticks=locator)
plt.show()

fig,ax=plt.subplots()
levels=numpy.linspace(numpy.min(ElasticityDelta),numpy.max(ElasticityDelta),10)
locator = ticker.LinearLocator()
plot=ax.contourf(xPosition,yPosition,ElasticityDelta,levels=levels)
ax.set_xlabel("X-Position [mm]")
ax.set_ylabel("Y-Position [mm]")
ax.set_title("Chain Elasticity Displacement vs. Position on Work Surface")
clb  = plt.colorbar(plot, format='%.e', ticks=locator)
plt.show()

fig,ax=plt.subplots()
levels=numpy.linspace(numpy.min(TotalDelta),numpy.max(TotalDelta),10)
locator = ticker.LinearLocator()
plot=ax.contourf(xPosition,yPosition,TotalDelta,levels=levels)
ax.set_xlabel("X-Position [mm]")
ax.set_ylabel("Y-Position [mm]")
ax.set_title("Total, Sum of Chain Sag and Chain Stretch, Displacement vs. Position on Work Surface")
clb  = plt.colorbar(plot, format='%.e', ticks=locator)
plt.show()