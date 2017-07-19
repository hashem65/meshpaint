#!/usr/bin/env python

#> \file 
#> \author David Ladd, Reused: Hashem Yousefi 
#> \brief This is an example to use linear fitting to fit the beginning of linear heart tube.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>
### to be used for segmenting embryonic heart and fitting with an initial meshes.
#<

import sys, os
import exfile
import numpy
from numpy import linalg
import math
import random
from random import randint 
# running optimization parameters ... 
# Intialise OpenCMISS/iron 
from opencmiss.iron import iron

# defining the output file to be written in the ExDataFile
def writeExdataFile(filename,dataPointLocations,dataErrorVector,dataErrorDistance,offset):
    "Writes data points to an exdata file"

    numberOfDimensions = dataPointLocations[1].shape[0]
    try:
        f = open(filename,"w")
        if numberOfDimensions == 1:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=2, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=3, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 2:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=3, #Derivatives=0, #Versions=1
  y.  Value index=4, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=5, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 3:
             header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
  x.  Value index=3, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=4, #Derivatives=0, #Versions=1
  y.  Value index=5, #Derivatives=0, #Versions=1
  z.  Value index=6, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=7, #Derivatives=0, #Versions=1
'''
        f.write(header)

        numberOfDataPoints = len(dataPointLocations)
        for i in range(numberOfDataPoints):
            line = " Node: " + str(offset+i+1) + '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataPointLocations[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataErrorVector[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            line = ' ' + str(dataErrorDistance[i])
            f.write(line)
            line = '\n'
            f.write(line)
        f.close()
            
    except IOError:
        print ('Could not open file: ' + filename)



# initial mesh 
# the location of  nodes for the mesh  
numberOfDimensions = 3
numberOfGaussXi = 3 
numberOfCircumfrentialElementsPerQuarter = 2
numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfLengthElements = 8
numberOfLengthNodes = numberOfLengthElements+1
numberOfWallElements = 1
numberOfWallNodes = numberOfWallElements+1
origin = [0.0,0.0,0.0]
meshOrigin = [0.0,0.0,0.0]
print "mesh resolution and parameters fixed"
#Epi = True
Epi = False 
 
#Endo = True
Endo = False
  
Myo = True
#Myo = False
   
#CJ = True
CJ = False
   
if (Epi):    
    numberOfDataPoints = 1024
elif (Endo):
    numberOfDataPoints = 473
elif (Myo):
    numberOfDataPoints = 920
else: 
    numberOfDataPoints = 654

rmsInitial = 30 
print "rmsInitial = ", rmsInitial
kappa = 0.001
tau = 0.001
ratioRefin = 0.5
criteriaDistance = 60 
startpoint = 0
fixInterior = True
zeroTolerance = 0.00001
hermite = True
startIteration = 0
if startIteration > 1:
    exfileMesh = True
    exnode = exfile.Exnode("DeformedGeometry" + str(startIteration-1) + ".part0.exnode")
    exelem = exfile.Exelem("UndeformedGeometry.part0.exelem")
else:
    exfileMesh = False
print "other parameters setted up "

manualNodePoints = numpy.zeros((numberOfLengthNodes,8,3,2))
manualNodePoints[0,0,:,0] = [528,318,-5]
manualNodePoints[0,1,:,0] = [521,321,-5]
manualNodePoints[0,2,:,0] = [510,324,-5]
manualNodePoints[0,3,:,0] = [502,322,-5]
manualNodePoints[0,4,:,0] = [497,317,-5]
manualNodePoints[0,5,:,0] = [500,311,-5]
manualNodePoints[0,6,:,0] = [508,305,-5]
manualNodePoints[0,7,:,0] = [522,308,-5]

manualNodePoints[1,0,:,0] = [529,314,25]
manualNodePoints[1,1,:,0] = [519,319,25]
manualNodePoints[1,2,:,0] = [509,322,25]
manualNodePoints[1,3,:,0] = [497,319,25]
manualNodePoints[1,4,:,0] = [494,311,25]
manualNodePoints[1,5,:,0] = [501,306,25]
manualNodePoints[1,6,:,0] = [511,306,25]
manualNodePoints[1,7,:,0] = [523,307,25]

manualNodePoints[2,0,:,0] = [548,300,130]
manualNodePoints[2,1,:,0] = [541,312,127]
manualNodePoints[2,2,:,0] = [526,330,131]
manualNodePoints[2,3,:,0] = [512,311,128]
manualNodePoints[2,4,:,0] = [505,296,130]
manualNodePoints[2,5,:,0] = [520,290,130]
manualNodePoints[2,6,:,0] = [531,287,130]
manualNodePoints[2,7,:,0] = [544,285,130]

manualNodePoints[3,0,:,0] = [596,270,255]
manualNodePoints[3,1,:,0] = [573,302,245]
manualNodePoints[3,2,:,0] = [546,332,217]
manualNodePoints[3,3,:,0] = [498,296,242]
manualNodePoints[3,4,:,0] = [467,262,255]
manualNodePoints[3,5,:,0] = [491,246,255]
manualNodePoints[3,6,:,0] = [534,253,255]
manualNodePoints[3,7,:,0] = [564,252,255]

manualNodePoints[4,0,:,0] = [625,230,380]
manualNodePoints[4,1,:,0] = [591,250,380]
manualNodePoints[4,2,:,0] = [538,236,386]
manualNodePoints[4,3,:,0] = [488,267,380]
manualNodePoints[4,4,:,0] = [413,265,380]
manualNodePoints[4,5,:,0] = [455,230,380]
manualNodePoints[4,6,:,0] = [547,194,380]
manualNodePoints[4,7,:,0] = [587,201,380]

manualNodePoints[5,0,:,0] = [672,250,490]
manualNodePoints[5,1,:,0] = [610,262,490]
manualNodePoints[5,2,:,0] = [578,270,496]
manualNodePoints[5,3,:,0] = [551,285,496]
manualNodePoints[5,4,:,0] = [471,296,490]
manualNodePoints[5,5,:,0] = [495,245,490]
manualNodePoints[5,6,:,0] = [523,213,490]
manualNodePoints[5,7,:,0] = [590,232,490]

manualNodePoints[6,0,:,0] = [637,249,580]
manualNodePoints[6,1,:,0] = [627,270,580]
manualNodePoints[6,2,:,0] = [600,294,581]
manualNodePoints[6,3,:,0] = [587,301,582]
manualNodePoints[6,4,:,0] = [562,302,580]
manualNodePoints[6,5,:,0] = [557,271,580]
manualNodePoints[6,6,:,0] = [582,260,580]
manualNodePoints[6,7,:,0] = [599,231,580]

manualNodePoints[7,0,:,0] = [609,307,610]
manualNodePoints[7,1,:,0] = [608,328,610]
manualNodePoints[7,2,:,0] = [588,328,610]
manualNodePoints[7,3,:,0] = [570,308,610]
manualNodePoints[7,4,:,0] = [563,282,610]
manualNodePoints[7,5,:,0] = [574,273,610]
manualNodePoints[7,6,:,0] = [585,272,610]
manualNodePoints[7,7,:,0] = [604,285,610]

manualNodePoints[8,0,:,0] = [610,304,635]
manualNodePoints[8,1,:,0] = [609,319,635]
manualNodePoints[8,2,:,0] = [586,323,635]
manualNodePoints[8,3,:,0] = [569,307,635]
manualNodePoints[8,4,:,0] = [563,282,635]
manualNodePoints[8,5,:,0] = [574,273,635]
manualNodePoints[8,6,:,0] = [585,272,635]
manualNodePoints[8,7,:,0] = [602,284,635]

# node positions of the outer surface ... 
manualNodePoints[0,0,:,1] = [600,311,-5]
manualNodePoints[0,1,:,1] = [572,336,-5]
manualNodePoints[0,2,:,1] = [518,345,-5]
manualNodePoints[0,3,:,1] = [461,340,-5]
manualNodePoints[0,4,:,1] = [436,312,-5]
manualNodePoints[0,5,:,1] = [455,288,-5]
manualNodePoints[0,6,:,1] = [520,278,-5]
manualNodePoints[0,7,:,1] = [573,287,-5]

manualNodePoints[1,0,:,1] = [611,285,25]
manualNodePoints[1,1,:,1] = [580,334,25]
manualNodePoints[1,2,:,1] = [520,350,25]
manualNodePoints[1,3,:,1] = [435,355,25]
manualNodePoints[1,4,:,1] = [395,325,25]
manualNodePoints[1,5,:,1] = [430,285,25]
manualNodePoints[1,6,:,1] = [513,252,25]
manualNodePoints[1,7,:,1] = [575,260,25]

manualNodePoints[2,0,:,1] = [661,283,137]
manualNodePoints[2,1,:,1] = [656,367,158]
manualNodePoints[2,2,:,1] = [534,397,158]
manualNodePoints[2,3,:,1] = [391,373,159]
manualNodePoints[2,4,:,1] = [331,300,125]
manualNodePoints[2,5,:,1] = [390,245,135]
manualNodePoints[2,6,:,1] = [495,220,135]
manualNodePoints[2,7,:,1] = [612,217,136]

manualNodePoints[3,0,:,1] = [702,232,259]
manualNodePoints[3,1,:,1] = [657,309,269]
manualNodePoints[3,2,:,1] = [530,366,252]
manualNodePoints[3,3,:,1] = [392,347,260]
manualNodePoints[3,4,:,1] = [315,290,260]
manualNodePoints[3,5,:,1] = [376,203,264]
manualNodePoints[3,6,:,1] = [505,181,257]
manualNodePoints[3,7,:,1] = [656,176,260]

manualNodePoints[4,0,:,1] = [775,210,380]
manualNodePoints[4,1,:,1] = [719,295,376]
manualNodePoints[4,2,:,1] = [572,312,380]
manualNodePoints[4,3,:,1] = [415,328,378]
manualNodePoints[4,4,:,1] = [313,270,380]
manualNodePoints[4,5,:,1] = [350,190,380]
manualNodePoints[4,6,:,1] = [500,140,380]
manualNodePoints[4,7,:,1] = [681,145,377]

manualNodePoints[5,0,:,1] = [819,227,490]
manualNodePoints[5,1,:,1] = [755,323,499]
manualNodePoints[5,2,:,1] = [590,320,490]
manualNodePoints[5,3,:,1] = [453,377,487]
manualNodePoints[5,4,:,1] = [305,305,490]
manualNodePoints[5,5,:,1] = [324,205,491]
manualNodePoints[5,6,:,1] = [514,164,497]
manualNodePoints[5,7,:,1] = [706,145,490]

manualNodePoints[6,0,:,1] = [839,239,586]
manualNodePoints[6,1,:,1] = [778,332,585]
manualNodePoints[6,2,:,1] = [606,369,590]
manualNodePoints[6,3,:,1] = [461,383,577]
manualNodePoints[6,4,:,1] = [335,305,586]
manualNodePoints[6,5,:,1] = [385,205,590]
manualNodePoints[6,6,:,1] = [521,165,587]
manualNodePoints[6,7,:,1] = [716,168,590]

manualNodePoints[7,0,:,1] = [832,251,622]
manualNodePoints[7,1,:,1] = [785,345,618]
manualNodePoints[7,2,:,1] = [605,367,616]
manualNodePoints[7,3,:,1] = [445,393,608]
manualNodePoints[7,4,:,1] = [342,318,620]
manualNodePoints[7,5,:,1] = [385,205,615]
manualNodePoints[7,6,:,1] = [536,166,620]
manualNodePoints[7,7,:,1] = [711,172,621]

manualNodePoints[8,0,:,1] = [814,260,645]
manualNodePoints[8,1,:,1] = [785,361,644]
manualNodePoints[8,2,:,1] = [603,363,640]
manualNodePoints[8,3,:,1] = [432,413,633]
manualNodePoints[8,4,:,1] = [323,346,643]
manualNodePoints[8,5,:,1] = [391,211,644]
manualNodePoints[8,6,:,1] = [541,188,650]
manualNodePoints[8,7,:,1] = [703,177,653]

if (Myo):
    for i in range (9):
        for j in range (8):
            for k in range(3):
                manualNodePoints[i,j,k,0]=ratioRefin*manualNodePoints[i,j,k,0]+(1-ratioRefin)*manualNodePoints[i,j,k,1]

#calculating the derivatives 
difference = numpy.zeros((numberOfLengthNodes,8,3,2))
differenceAverage = numpy.zeros((numberOfLengthNodes,8,3,2))
circumDeriv = numpy.zeros((numberOfLengthNodes,8,3,2))
directDeriv = numpy.zeros((numberOfLengthNodes,8,3,2))
lengthDeriv = numpy.zeros((numberOfLengthNodes,8,3,2))
#circumferential derivative to be calculated 
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            if (i<7):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i+1,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,0,m,k]-manualNodePoints[j,7,m,k]
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            if (i<7):
                for m in range (3):
                    differenceAverage[j,i+1,m,k]=(difference[j,i+1,m,k]+difference[j,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,0,m,k]=(difference[j,0,m,k]+difference[j,7,m,k])/2
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            for m in range (3):
                circumDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# derivative of the length direction
for k in range (2):
    for i in range (8):
        for j in range (numberOfLengthNodes):
            if (j<numberOfLengthNodes-1):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j+1,i,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i,m,k]-manualNodePoints[j-1,i,m,k]
for k in range (2):
    for i in range (8):
        for j in range (numberOfLengthNodes):
            if (j == 0):
                for m in range (3): 
                    differenceAverage[j,i,m,k]=difference[j,i,m,k]
            if (j<numberOfLengthNodes-1):
                for m in range (3):
                    differenceAverage[j+1,i,m,k]=(difference[j,i,m,k]+difference[j+1,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,i,m,k]=difference[j-1,i,m,k]
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (8):
            for m in range (3):
                lengthDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# the derivatives of the wall direction is defined in the below lines ... 
for i in range (8):
    for j in range (numberOfLengthNodes):
        for m in range (3):
            for k in range (2):
                difference[j,i,m,k] = manualNodePoints[j,i,m,1] - manualNodePoints[j,i,m,0]
for i in range (8):
    for j in range (numberOfLengthNodes):

        for k in range (2):
            for m in range (3):
                directDeriv[j,i,m,k] = difference[j,i,m,k]/math.sqrt(math.pow(difference[j,i,0,k],2) + math.pow(difference[j,i,1,k],2) + math.pow(difference[j,i,2,k],2))

# setting up optimization parameters and loops 


(coordinateSystemUserNumber, regionUserNumber, basisUserNumber, generatedMeshUserNumber, meshUserNumber, decompositionUserNumber, geometricFieldUserNumber,
        equationsSetFieldUserNumber,dependentFieldUserNumber,independentFieldUserNumber,dataPointFieldUserNumber,materialFieldUserNumber,
        analyticFieldUserNumber,dependentDataFieldUserNumber,dataPointsUserNumber,dataProjectionUserNumber,equationsSetUserNumber,problemUserNumber) = range(1,19)
    # Get the computational nodes information
    #print dir(iron),'\n\n'
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()
    # Create a RC CS
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()
    # Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()
   # define a basis 
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfDimensions
if hermite:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
else:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()
print "CS, Region and basis setted up"
# ====================================================
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
if (numberOfWallElements == 0): 
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements
elif (numberOfWallElements == 1): 
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements
else:
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements
if (exfileMesh):
    # Read previous mesh
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, numberOfDimensions)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(exelem.num_elements)
    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, exnode.num_nodes)
    nodes.CreateFinish()
    # Define elements for the mesh
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    for elem in exelem.elements:
        elements.NodesSet(elem.number, elem.nodes)
    elements.CreateFinish()
    mesh.CreateFinish()
else:
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber,region,3)
    mesh.origin = meshOrigin
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(numberOfElements)
# Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region,numberOfNodes)
    nodes.CreateFinish()
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    elementNumber = 0
    for wallElementIdx in range(1,numberOfWallElements+1):
       for lengthElementIdx in range(1,numberOfLengthElements+1):
            for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                elementNumber = elementNumber + 1
                localNode1 = circumfrentialElementIdx + (lengthElementIdx - 1)*numberOfCircumfrentialElements + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                if circumfrentialElementIdx == numberOfCircumfrentialElements:
                    localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                        (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                else: 
                    localNode2 = localNode1 + 1
                localNode3 = localNode1 + numberOfCircumfrentialNodes
                localNode4 = localNode2 + numberOfCircumfrentialNodes
                localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
#		print "Element Number = ",elementNumber
#         	print "Node numbers of the element", localNode1, localNode2, localNode3, localNode4, localNode5, localNode6, localNode7, localNode8 
                elements.NodesSet(elementNumber,localNodes)  
    elements.CreateFinish()
    mesh.CreateFinish() 
# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()
print "mesh decomposition finished"

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
for dimension in range(3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimension+1,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()
# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)    
numberOfNodes = nodes.numberOfNodes

# Get or calculate geometric parameters
if (exfileMesh):
    # Read the geometric field from the exnode file
    for node_num in range(1, exnode.num_nodes + 1):
        for derivative in range(1,9):
            version = 1
            for component in range(1, numberOfDimensions + 1):
                component_name = ["x", "y", "z"][component - 1]
                value = exnode.node_value("Coordinate", component_name, node_num, derivative)
                geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                      version, derivative, node_num, component, value)
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
else:
    # Create the geometric field
    for wallNodeIdx in range(1,numberOfWallNodes+1):
        for lengthNodeIdx in range(1,numberOfLengthNodes+1):
            for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
                x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                xtangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                ytangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                ztangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                xnormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                ynormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                znormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                zxnormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                zynormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                zznormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,1,x)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,2,y)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,3,z)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,ztangent)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,zxnormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,zynormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,zznormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
                geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,znormal)
    # Update the geometric field
    geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    # Export undeformed mesh geometry
######    print("Writing undeformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("UndeformedGeometry","FORTRAN")
    fields.ElementsExport("UndeformedGeometry","FORTRAN")
    fields.Finalise()

#=================================================================
# Data Points
#=================================================================
# Create the data points
dataPoints = iron.DataPoints()
dataPoints.CreateStart(dataPointsUserNumber,region,numberOfDataPoints)
dataPointLocations = numpy.zeros((numberOfDataPoints,3))
######print("Number of data points: " + str(numberOfDataPoints))
# reading from a text file containing the point clouds   
if (Epi):
    with open("EpiDataPoints.txt", "r") as ins:
	    arrayOfInputData = []
	    for line in ins:
		    arrayOfInputData.append(line)
elif (Endo): 
    with open("EndoDataPoints.txt", "r") as ins:
	    arrayOfInputData = []
	    for line in ins:
		    arrayOfInputData.append(line)
elif (Myo): 
    with open("MyoDataPoints.txt", "r") as ins:
	    arrayOfInputData = []
	    for line in ins:
		    arrayOfInputData.append(line)
else: 
    with open("CJDataPoints.txt", "r") as ins:
	    arrayOfInputData = []
	    for line in ins:
		    arrayOfInputData.append(line)

x = 0.0
y = 0.0
z = 0.0
#for i in range (startpoint, numberOfDataPoints + startpoint):
for i in range (numberOfDataPoints):
#for i in range (200):
	for j in range (5):
		sample = arrayOfInputData[i*5 + j]
#		sample = arrayOfInputData[i*25 + j]
		if (math.fmod(j,5) == 1):
			x = float (sample[12:25])				
		elif (math.fmod(j,5) == 2):
			y = float (sample[12:25])
		elif (math.fmod(j,5) == 3):
			z = float (sample[12:17])
#		dataPointLocations[i - startpoint,:] = [x,y,z]
		dataPointLocations[i,:] = [x,y,z]


# Set up data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.PositionSet(dataPointId,dataList)
dataPoints.CreateFinish()
 
#=================================================================
# Data Projection on Geometric Field
#=================================================================
######print("Projecting data points onto geometric field")
candidateElements = range(1,numberOfElements+1)

#    candidateElements = range(1,numberOfElements+1)
if (Epi): 
    candidateFaceNormals = iron.ElementNormalXiDirections.PLUS_XI3*numpy.ones(numberOfElements,dtype=numpy.int32)
else: 
    candidateFaceNormals = iron.ElementNormalXiDirections.MINUS_XI3*numpy.ones(numberOfElements,dtype=numpy.int32)

# Set up data projection
dataProjection = iron.DataProjection()
dataProjection.CreateStart(dataProjectionUserNumber,dataPoints,geometricField,iron.FieldVariableTypes.U)
#dataProjection.projectionType = iron.DataProjectionProjectionTypes.ALL_ELEMENTS
dataProjection.projectionType = iron.DataProjectionProjectionTypes.BOUNDARY_FACES
dataProjection.ProjectionCandidateFacesSet(candidateElements,candidateFaceNormals)
#dataProjection.ProjectionDataCandidateFacesSet([1,2,3],[1,2],[iron.ElementNormalXiDirections.PLUS_XI3,iron.ElementNormalXiDirections.PLUS_XI3])
dataProjection.CreateFinish()

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition topology for data projection
decomposition.TopologyDataProjectionCalculate()

# Cancel some projections
dataProjection.ProjectionCancelByDistance(iron.DataProjectionDistanceRelations.GREATER_EQUAL,criteriaDistance)

# Output data projection results
dataProjection.ResultAnalysisOutput("ProjectionAnalysis")

rmsError=dataProjection.ResultRMSErrorGet()
######print("RMS error = "+ str(rmsError))
rmsError1 = rmsError
# Output the .exdata file.                                           
dataErrorVector = numpy.zeros((numberOfDataPoints,3))
dataErrorDistance = numpy.zeros(numberOfDataPoints)
for elementIdx in range(1,numberOfElements+1):
    numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementIdx)
    for dataPointIdx in range(1,numberOfProjectedDataPoints+1):
        dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementIdx,dataPointIdx)
        errorVector = dataProjection.ResultProjectionVectorGet(dataPointNumber,3)
        dataErrorVector[dataPointNumber-1,0]=errorVector[0]
        dataErrorVector[dataPointNumber-1,1]=errorVector[1]
        dataErrorVector[dataPointNumber-1,2]=errorVector[2]
        errorDistance = dataProjection.ResultDistanceGet(dataPointNumber)
        dataErrorDistance[dataPointNumber-1]=errorDistance
 
# write data points to exdata file for CMGUI
offset = 0
######writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,dataErrorVector,dataErrorDistance,offset)

######print("Projection complete")
#exit(0)

#=================================================================
# Equations Set
#=================================================================
# Create vector fitting equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.DATA_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.DATA_POINT_FITTING, 
 			     iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

#=================================================================
# Dependent Field
#=================================================================
# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,numberOfDimensions)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions)
equationsSet.DependentCreateFinish()
# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

# Initialise dependent field to undeformed geometric field
for component in range (1,numberOfDimensions+1):
    geometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, component)

#=================================================================
# Independent Field
#=================================================================
# Create data point field (independent field, with vector values stored at the data points)
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U,"data point vector")
independentField.VariableLabelSet(iron.FieldVariableTypes.V,"data point weight")
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()
# Initialise data point vector field to 0
#independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
# Initialise data point weight field to 1
#independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)
# loop over each element's data points and set independent field values to data point locations on surface of the sphere
for element in range(numberOfElements):
    elementId = element + 1
    elementDomain = decomposition.ElementDomainGet(elementId)
    if (elementDomain == computationalNodeNumber):
        numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementId)
        for dataPoint in range(numberOfProjectedDataPoints):
            dataPointId = dataPoint + 1
            dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementId,dataPointId)
            dataList = dataPoints.PositionGet(dataPointNumber,3)
            # set data point field values
            for component in range(numberOfDimensions):
                componentId = component + 1
                dataPointNumberIndex = dataPointNumber - 1
                value = dataList[component]
                independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementId,dataPointId,componentId,value)

#=================================================================
# Material Field
#=================================================================
# Create material field (Sobolev parameters)
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Smoothing Parameters")
equationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================
# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.FULL
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================
# Create fitting problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.STATIC_FITTING]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE # NONE / MATRIX
#solver.outputType = iron.SolverOutputTypes.MATRIX # NONE / MATRIX
solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.linearType = iron.LinearSolverTypes.DIRECT
#solver.LibraryTypeSet(iron.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
#solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
solver.linearIterativeAbsoluteTolerance = 1.0E-10
solver.linearIterativeRelativeTolerance = 1.0E-05
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
#solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================


# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
'''
for nodeIdx in range(145,145+8):
    for componentIdx in range(1,4):
        for derivativeIdx in range(1,9):
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
for nodeIdx in range(217-16,217):
    for componentIdx in range(1,4):
        for derivativeIdx in range(1,9):
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
'''
if (Epi):
    for nodeIdx in range(73,73+8):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                 boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
    for nodeIdx in range(144-16,145):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
else: 
    for nodeIdx in range(1,1+8):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
    for nodeIdx in range(73-8,73):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)

# for nodeIdx in range(numberOfLengthNodes*numberOfCircumfrentialNodes+1,numberOfLengthNodes*numberOfCircumfrentialNodes*2+1):
if (Epi):
    for nodeIdx in range(1,numberOfLengthNodes*numberOfCircumfrentialNodes*(numberOfWallNodes-1)+1):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)
else:
    for nodeIdx in range(numberOfLengthNodes*numberOfCircumfrentialNodes+1,numberOfLengthNodes*numberOfCircumfrentialNodes*numberOfWallNodes+1):
        for componentIdx in range(1,4):
            for derivativeIdx in range(1,9):
                boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
                                       1,derivativeIdx,nodeIdx,componentIdx,
                                       iron.BoundaryConditionsTypes.FIXED,0.0)

#for nodeIdx in range(5,9):
#    for componentIdx in range(1,4):
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                                   1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,nodeIdx,componentIdx,
#                                   iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeIdx,2,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeIdx,1,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeIdx,2,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
#    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,
#                               1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeIdx,3,
#                               iron.BoundaryConditionsTypes.FIXED,0.0)
solverEquations.BoundaryConditionsCreateFinish()

derivativeVector=[0.0,0.0,0.0,0.0]
numberOfIterations = 6
for iteration in range (startIteration,startIteration+numberOfIterations+1):
    # Solve the problem
    print("Solving fitting problem, iteration: " + str(iteration))
    problem.Solve()
    # Normalise derivatives
    for nodeIdx in range(1,numberOfNodes+1):
      for derivativeIdx in [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
          length=0.0
          for componentIdx in range(1,4):
              derivativeVector[componentIdx]=dependentField.ParameterSetGetNode(iron.FieldVariableTypes.U,
                                                                                iron.FieldParameterSetTypes.VALUES,
                                                                                1,derivativeIdx,nodeIdx,componentIdx)
              length=length + derivativeVector[componentIdx]*derivativeVector[componentIdx]
          length=math.sqrt(length)
          for componentIdx in range(1,4):
              value=derivativeVector[componentIdx]/length
              dependentField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,derivativeIdx,nodeIdx,componentIdx,value)
    # Copy dependent field to geometric 
    for componentIdx in range(1,numberOfDimensions+1):
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx,geometricField,
                                                                iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx)
    # Reproject
    dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
    rmsError=dataProjection.ResultRMSErrorGet()
    print("RMS error = "+ str(rmsError))
    if (rmsError < rmsInitial):
        print("RMSFitting = "+ str(rmsError))
        print "criteriaDistance = ", criteriaDistance, "kappa =", kappa, "tau =", tau, "ratioRefin=", ratioRefin
        print "iterationNumber = ", iteration
        rmsInitial = rmsError
    # Export fields
    print("Writing deformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# iterating through a random variable of tau and kappa ... 
# different values of criteria Distance and Refin Ratio will be considered later
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def findRMS(optimin,**kwargs):
    cmissobject=kwargs['cmiss']
    geometricField = cmissobject['geometricField']
    materialField = cmissobject['materialField']
    dependentField= cmissobject['dependentField']
    dataProjection= cmissobject['dataProjection']
    numberOfWallNodes= cmissobject['numberOfWallNodes']
    numberOfCircumfrentialNodes= cmissobject['numberOfCircumfrentialNodes']
    g = []
    try:
        kappa, tau ,ratioRefin, criteriaDistance = optimin[0],optimin[1],optimin[2],optimin[3]
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
                    nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
                    x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                    y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                    z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                    xtangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                    ytangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                    ztangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                    xnormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                    ynormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                    znormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                    zxnormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
                    zynormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
                    zznormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,1,nodeNumber,1,x)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,1,nodeNumber,2,y)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,1,nodeNumber,3,z)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,ztangent)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,zxnormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,zynormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,zznormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
                    geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,znormal)
        # Update the geometric field
        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
        materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)
        derivativeVector=[0.0,0.0,0.0,0.0]
        numberOfIterations = 6
        rmsError = 10e16
        for iteration in range (startIteration,startIteration+numberOfIterations+1):
            # Solve the problem
            #print("Solving fitting problem, iteration: " + str(iteration))
            problem.Solve()
            # Normalise derivatives
            for nodeIdx in range(1,numberOfNodes+1):
              for derivativeIdx in [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
                  length=0.0
                  for componentIdx in range(1,4):
                      derivativeVector[componentIdx]=dependentField.ParameterSetGetNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,derivativeIdx,nodeIdx,componentIdx)
                      length=length + derivativeVector[componentIdx]*derivativeVector[componentIdx]
                  length=math.sqrt(length)
                  for componentIdx in range(1,4):
                      value=derivativeVector[componentIdx]/length
                      dependentField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,derivativeIdx,nodeIdx,componentIdx,value)
            # Copy dependent field to geometric 
            for componentIdx in range(1,numberOfDimensions+1):
                dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, componentIdx,geometricField,
                                                                        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, componentIdx)
            # Reproject
            dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
    #        mesh.TopologyDataPointsCalculateProjection(dataProjection)
            # Create decomposition topology for data projection
    #        decomposition.TopologyDataProjectionCalculate()
            # Cancel some projections
            dataProjection.ProjectionCancelByDistance(iron.DataProjectionDistanceRelations.GREATER_EQUAL,criteriaDistance)
            # Output data projection results
            dataProjection.ResultAnalysisOutput("ProjectionAnalysis")
            rmsError=min([dataProjection.ResultRMSErrorGet(),rmsError])
        f=rmsError
        fail = 0
        return f,g,fail
    except:
        fail = 1
        f = 1e16
        return f,g,fail


cmissobject = dict()
cmissobject['geometricField'] =    geometricField 
cmissobject['materialField'] =    materialField 
cmissobject['dependentField'] =    dependentField
cmissobject['dataProjection']=    dataProjection
cmissobject['numberOfWallNodes'] =  numberOfWallNodes
cmissobject['numberOfCircumfrentialNodes'] =  numberOfCircumfrentialNodes

x0 = [0.1,0.1,0.8,50.0]

from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP

fkeys = ['kappa','tau','wallRefineXi','Cuttoff']
bnds = [(1e-7,1.0),(1e-7,1.0),(0.1,0.9),(20.0,100.0)]
opt_prob = Optimization('Mesh Fitting',findRMS)
for i,v in enumerate(bnds):
    opt_prob.addVar(fkeys[i],'c',lower=v[0],upper=v[1],value=x0[i])
opt_prob.addObj('f')
print opt_prob

# Global Optimization
nsga2 = NSGA2()
nsga2(opt_prob,cmiss=cmissobject)
print opt_prob.solution(0)

# Local Optimization Refinement
slsqp = SLSQP()
slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
print opt_prob.solution(0).solution(0)


#-----------------------------------------------------------------
iron.Finalise()
