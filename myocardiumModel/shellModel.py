# -*- coding: utf-8 -*-
"""
Created on Wed July 29 13:16:10 2020

@author: bieZY (e0348827@u.nus.edu)
Objective: ideal left ventricle myocardium ideal model based on the middle layer shell
Method: 
Modules: numpy, scipy, trimesh, meshplex

History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: jorry.zhengyu@gmail.com         29July2020           -V1.0.0 Created, test version
  Author: jorry.zhengyu@gmail.com         26AUGU2020           -V1.0.1, test version, improve innerSurface_ES function to solve inner_ES STL volume error, point_ES error
"""
print('shellModel test version 1.0.1')

#import os
import sys
#sys.path.insert(0, "/home/yaplab/Programming/python3")
#sys.path.insert(0,"C:\\Users\\bieZY\\Documents\\Python3")
sys.path.insert(0,"C:\\Users\\bieZY\\GoogleDrive\\Work - NUS\\2.PhD Project\\2.PyCode")
sys.path.insert(0,"C:\\zhengyu\\PyCode")
import numpy as np
#import scipy as sp
import scipy as sp
import scipy.spatial as spatial
import trimesh
import meshplex

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

class shellModel:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self):
        '''
        Initialize shellModel class
        '''
        # shell part
        self.sliceNum = None
        self.slicePoint = None      # points number of each slice,last slice must be 1 (apex point)
        self.sliceRadius = None     # radius of each slice,last slice must be 0 (apex point)
        self.sliceInterval = None   # interval between nearby slices
        self.shellCoord = None      # all the points
        self.shellPoint = None      # one point each slice to calculate strain
        self.connectivity_shell = None
        
        self.sliceInterval_ES = None
        self.sliceRadius_ES = None
        self.shellCoord_ES = None      # all the points
        self.shellPoint_ES = None      # one point each slice to calculate strain
        self.shellStretch = None
        self.shellStrain = None
        
        # inner part
        self.thickness_inner = None
        self.sliceNum_inner = None
        self.slicePoint_inner = None
        self.sliceRadius_inner = None
        self.sliceInterval_inner = None
        self.innerCoord = None
        self.innerPoint = None
        self.connectivity_inner = None
        
#        self.thickness_inner_ES = None
#        self.sliceNum_inner_ES = None
#        self.slicePoint_inner_ES = None
        self.sliceRadius_inner_ES = None
        self.sliceInterval_inner_ES = None
        self.innerCoord_ES = None
        self.innerPoint_ES = None
        
        # outer part
        self.thickness_outer = None
        self.sliceNum_outer = None
        self.slicePoint_outer = None
        self.sliceRadius_outer = None
        self.sliceInterval_outer = None
        self.outerCoord = None
        self.outerPoint = None
        self.connectivity_outer = None
        
        self.sliceRadius_outer_ES = None
        self.sliceInterval_outer_ES = None
        self.outerCoord_ES = None
        self.outerPoint_ES = None
        
        self.connectivity = None
        
        print('shellModel_init_ done')
        
    def surfaceGeneration(self,sliceNum=None,slicePoint=None,sliceRadius=None,sliceInterval=None):
        temp = []       # all the points
        temp3 = []      # one point each slice
        for i in range(sliceNum):
            if i==sliceNum-1:
                temp.append([0,0,0])
                temp3.append([0,0,0])
                continue
            angle = np.pi*2/slicePoint[i]
            radius = sliceRadius[i]
            #interval = self.sliceInterval[i]
            z = np.sum(sliceInterval[i:])
            
            for j in range(int(slicePoint[i])):
                x = radius*np.cos(j*angle)
                y = radius*np.sin(j*angle)
                #z = np.sum(self.sliceInterval[i:])
                #z = interval*(sliceNum - i -1)
                temp.append([x,y,z])
            temp3.append([radius,0,z])
            
        temp2 = np.array(temp)
        temp4 = np.array(temp3)
        return temp2,temp4
        
    def shellSurface(self,sliceNum=None,slicePoint=None,sliceRadius=None,sliceInterval=None):
        self.sliceNum = sliceNum
        if type(slicePoint) == int and sliceNum>1:
            temp = np.ones(sliceNum)*slicePoint
            temp[-1] = 1
            self.slicePoint = temp
        else:
            self.slicePoint = slicePoint
            self.slicePoint[-1] = 1
        for i in range(sliceNum-1):
            if self.slicePoint[i]<3:
                self.slicePoint[i] = 3
        if type(sliceRadius) == int and sliceNum>1:
            temp = np.ones(sliceNum)*sliceRadius
            temp[-1] = 0
            self.sliceRadius = temp
        else:
            self.sliceRadius = sliceRadius
        if type(sliceInterval) == int and sliceNum>1:
            temp = np.ones(sliceNum-1)*sliceInterval
            self.sliceInterval = temp
        else:
            self.sliceInterval = sliceInterval
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum, slicePoint=self.slicePoint, sliceRadius=self.sliceRadius, sliceInterval=self.sliceInterval)
        self.shellCoord = coord.copy()
        self.shellPoint = point.copy()
        '''
        temp = []
        for i in range(sliceNum):
            if i==sliceNum-1:
                temp.append([0,0,0])
                continue
            angle = np.pi*2/self.slicePoint[i]
            radius = self.sliceRadius[i]
            #interval = self.sliceInterval[i]
            z = np.sum(self.sliceInterval[i:])
            
            for j in range(int(self.slicePoint[i])):
                x = radius*np.cos(j*angle)
                y = radius*np.sin(j*angle)
                #z = np.sum(self.sliceInterval[i:])
                #z = interval*(sliceNum - i -1)
                temp.append([x,y,z])
        temp2 = np.array(temp)
        self.shellCoord = temp2.copy()
        '''
        print('shellSurface completed >_<')
        
    def innerSurface(self,sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, thickness=None):
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        if type(sliceRadius)==type(None):
            sliceRadius = self.sliceRadius.copy()
        if type(sliceInterval)==type(None):
            sliceInterval = self.sliceInterval.copy()
        
        if type(thickness) == int or type(thickness) == float and self.sliceNum>1:
#        if len(thickness) == 1 and self.sliceNum>1:
            temp = np.ones(self.sliceNum)*thickness
            self.thickness_inner = temp
        else:
            self.thickness_inner = thickness
        
        # use coordGeneration function
        radiusInner = []
        lastSlice = 0
        for i in range(sliceNum-1):
            if (sliceRadius[i] - self.thickness_inner[i])>=0 and (sliceRadius[i+1] - self.thickness_inner[i+1])<0:
                z = np.sum(sliceInterval[i:])
                sliceNum_inner = i+1
                radiusInner.append(0)
                lastSlice = 1
                break
            radius = sliceRadius[i] - self.thickness_inner[i]
            radiusInner.append(radius)
        if lastSlice==0:
            sliceNum_inner = sliceNum
            radiusInner.append(0)
            
        self.sliceNum_inner = sliceNum_inner
        temp = slicePoint[0:self.sliceNum_inner]
        self.slicePoint_inner = temp.copy()
        self.slicePoint_inner[-1] = 1
        self.sliceRadius_inner = np.array(radiusInner).copy()
        temp = sliceInterval[0:self.sliceNum_inner-1]
        if lastSlice==0:
            temp2=temp[-1]/2
            temp[-1]=temp2
        self.sliceInterval_inner = temp.copy()
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum_inner, slicePoint=self.slicePoint_inner, sliceRadius=self.sliceRadius_inner, sliceInterval=self.sliceInterval_inner)
        coord[:,2] = coord[:,2] + z
        point[:,2] = point[:,2] + z
        self.innerCoord = coord.copy()
        self.innerPoint = point.copy()
        
        '''
        temp = []
        radiusInner = []
        for i in range(sliceNum):
            if (sliceRadius[i] - self.thickness_inner[i])>=0 and (sliceRadius[i+1] - self.thickness[i+1])<0:
                z = np.sum(sliceInterval[i:])
                temp.append([0,0,z])
                self.sliceNum_inner = i+1
                radiusInner.append(0)
                break
            
            angle = np.pi*2/slicePoint[i]
            radius = sliceRadius[i] - self.thickness_inner[i]
            radiusInner.append(radius)
            z = np.sum(sliceInterval[i:])
            
            for j in range(int(slicePoint[i])):
                x = radius*np.cos(j*angle)
                y = radius*np.sin(j*angle)
                temp.append([x,y,z])
        
        temp2 = np.array(temp)
        self.innerCoord = temp2.copy()
        temp = slicePoint[0:self.sliceNum_inner]
        self.slicePoint_inner = temp.copy()
        self.slicePoint_inner[-1] = 1
        self.sliceRadius_inner = np.array(radiusInner).copy()
        temp = sliceInterval[0:self.sliceNum_inner-1]
        self.sliceInterval_inner = temp.copy()
        '''
        print('innerSurface completed >_<')
    
    def innerSurface_ES(self,sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, stretch_radial=None):
        '''
        simplified inner surface motion, use radial stretch to calculate radius (not same direction, especially apical part)
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum_inner
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint_inner.copy()
        if type(stretch_radial)==type(None):
            stretch_radial = self.shellStretch[:,2].copy()
        
        sliceInterval = self.sliceInterval_ES[:(sliceNum-1)]
        sliceRadius = self.sliceRadius_ES
        sliceRadius_inner = []
        flag = 0    # if sliceRadius_inner becomes negative
        for i in range(sliceNum-1):
            temp = sliceRadius[i]-stretch_radial[i]*(self.sliceRadius[i]-self.sliceRadius_inner[i])
            if temp<=0 and flag == 0:
                if i == 0:      # inner surface can't be smaller
                    return 1    # special return
                flag = 1
                negativeSlice = i
                lastRadius = sliceRadius_inner[-1]
                increment = lastRadius/(sliceNum-negativeSlice)     # modified radius
            if flag == 0:
                sliceRadius_inner.append(temp)
            else:
                temp = lastRadius-increment*(i-negativeSlice+1)
                sliceRadius_inner.append(temp)
                
                temp = sliceInterval[i-2]*(sliceRadius_inner[-2]-sliceRadius_inner[-1])/(sliceRadius_inner[-3]-sliceRadius_inner[-2])
                sliceInterval[i-1] = temp
        sliceRadius_inner.append(0)
        
        #if flag == 1:
            #temp = sliceInterval[-2]*(sliceRadius_inner[-2]-sliceRadius_inner[-1])/(sliceRadius_inner[-3]-sliceRadius_inner[-2])
            #sliceInterval[-1] = temp
        temp = sliceInterval[-2]*(sliceRadius_inner[-2]-sliceRadius_inner[-1])/(sliceRadius_inner[-3]-sliceRadius_inner[-2])
        sliceInterval[-1] = temp
        
        coord,point = self.surfaceGeneration(sliceNum=sliceNum, slicePoint=slicePoint, sliceRadius=sliceRadius_inner, sliceInterval=sliceInterval)
        
        temp = self.shellCoord_ES[0,2]-coord[0,2]
        coord[:,2] = coord[:,2] + temp
        point[:,2] = point[:,2] + temp
#        coord[:,2] = coord[:,2] + self.innerCoord[-1,2]
#        point[:,2] = point[:,2] + self.innerCoord[-1,2]
        self.sliceRadius_inner_ES = sliceRadius_inner.copy()
        self.sliceInterval_inner_ES = sliceInterval.copy()
        self.innerCoord_ES = coord.copy()
        self.innerPoint_ES = point.copy()
        return 0    # normal return
    
    def outerSurface(self,sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, thickness=None):
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        if type(sliceRadius)==type(None):
            sliceRadius = self.sliceRadius.copy()
        if type(sliceInterval)==type(None):
            sliceInterval = self.sliceInterval.copy()
        
        if type(thickness) == int or type(thickness) == float and self.sliceNum>1:
#        if len(thickness) == 1 and self.sliceNum>1:
            temp = np.ones(self.sliceNum)*thickness
            self.thickness_outer = temp
        else:
            self.thickness_outer = thickness
        
        # use coordGeneration function
        radiusOuter = []
        for i in range(sliceNum):
            radius = sliceRadius[i] + self.thickness_outer[i]
            radiusOuter.append(radius)
        sliceNum_outer = sliceNum+1
        radiusOuter.append(0)
            
        self.sliceNum_outer = sliceNum_outer
        temp = slicePoint[0:self.sliceNum_outer-1]
        temp[-1] = temp[-2]
        temp.append(1)
        self.slicePoint_outer = temp.copy()
        #self.slicePoint_outer[-2] = self.slicePoint_outer[-3]
        self.sliceRadius_outer = np.array(radiusOuter).copy()
        temp = sliceInterval
        lastInterval = (radiusOuter[-2]*temp[-1]/(radiusOuter[-3]-radiusOuter[-2]))*2/3
        temp.append(lastInterval)
        self.sliceInterval_outer = temp.copy()
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum_outer, slicePoint=self.slicePoint_outer, sliceRadius=self.sliceRadius_outer, sliceInterval=self.sliceInterval_outer)
        z = 0-self.sliceInterval_outer[-1]
        coord[:,2] = coord[:,2] + z
        point[:,2] = point[:,2] + z
        self.outerCoord = coord.copy()
        self.outerPoint = point.copy()
        print('outerSurface completed >_<')
    
    def outerSurface_ES(self,sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, stretch_radial=None):
        '''
        simplified outer surface motion, use radial stretch to calculate radius (not same direction, especially apical part)
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum_outer
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint_outer.copy()
        if type(stretch_radial)==type(None):
            stretch_radial = self.shellStretch[:,2].copy()
        
        sliceInterval = self.sliceInterval_ES.copy()
        #temp = self.sliceInterval_ES[-1]/self.sliceInterval[-1]*self.sliceInterval_outer[-1]
        temp = self.sliceInterval_outer[-1]
        sliceInterval.append(temp)
        sliceRadius = self.sliceRadius_ES
        sliceRadius_outer = []
        for i in range(sliceNum-2):
            sliceRadius_outer.append(sliceRadius[i]-stretch_radial[i]*(self.sliceRadius[i]-self.sliceRadius_outer[i]))
        temp = self.sliceRadius_outer[-2]*sliceRadius_outer[-1]/self.sliceRadius_outer[-3]
        sliceRadius_outer.append(temp)
        sliceRadius_outer.append(0)
        
        coord,point = self.surfaceGeneration(sliceNum=sliceNum, slicePoint=slicePoint, sliceRadius=sliceRadius_outer, sliceInterval=sliceInterval)
        
        temp = self.shellCoord_ES[0,2]-coord[0,2]
        coord[:,2] = coord[:,2] + temp
        point[:,2] = point[:,2] + temp
#        coord[:,2] = coord[:,2] + self.innerCoord[-1,2]
#        point[:,2] = point[:,2] + self.innerCoord[-1,2]
        self.sliceRadius_outer_ES = sliceRadius_outer.copy()
        self.sliceInterval_outer_ES = sliceInterval.copy()
        self.outerCoord_ES = coord.copy()
        self.outerPoint_ES = point.copy()
    
    def distCalc(self,coordA=None, coordB=None, mode='min'):
        '''
        calculate the distance between coordA and coordB: min (to do), max (to do)
        '''
        
        temp = spatial.distance_matrix(coordA,coordB)
        temp = temp.transpose()
        dist=temp.tolist()
        distValue=[]
        distIndex=[]
        if mode=='min':
            for i in range(len(coordB)):
                distValue.append(min(dist[i]))
                distIndex.append(dist[i].index(distValue[i]))
        elif mode=='max':
            for i in range(len(coordB)):
                distValue.append(max(dist[i]))
                distIndex.append(dist[i].index(distValue[i]))
        return distValue,distIndex
    
    def STLgeneration(self,coord=None, sliceNum=None, slicePoint=None, stlName=None):
        '''
        obtain point connection to generation STL
        saveName: name of STL file
        '''
        
        connect = []
        startPoint = 0
        endPoint = 0
        
        for i in range(sliceNum-1):
            endPoint += slicePoint[i]
            slicePoint_1 = coord[startPoint:endPoint,:]
            slicePoint_2 = coord[endPoint:(endPoint+slicePoint[i+1]),:]
            
            if slicePoint[i+1]==1:
                for j in range(slicePoint[i]):
                    if j+1 == slicePoint[i]:
                        connect.append([startPoint+j, endPoint, startPoint])
                    else:
                        connect.append([startPoint+j, endPoint, startPoint+j+1])
                break
            
            distMatrix = spatial.distance_matrix(slicePoint_1,slicePoint_2)
            #temp = distMatrix.copy()
            distMatrix_2 = distMatrix.copy()
            distMatrix_2[0:-1,:] = distMatrix[1:,:]
            distMatrix_2[-1,:] = distMatrix[0,:]
            distSum = distMatrix + distMatrix_2
            #np.around(distSum, decimals=6)
            
            for j in range(distSum.shape[0]):
                # upper slice, use first minDist to form a triangle
                temp = np.around(distSum[j,:],decimals=6)
                index = np.where(temp == np.min(temp))
                index_first = index[0][0]
                if j+1 == distSum.shape[0]:
                    connect.append([startPoint+j, endPoint+index_first, startPoint])
                else:
                    connect.append([startPoint+j, endPoint+index_first, startPoint+j+1])
                
            distMatrix_3 = distMatrix.copy()
            distMatrix_3[:, 0:-1] = distMatrix[:, 1:]
            distMatrix_3[:, -1] = distMatrix[:, 0]
            distSum = distMatrix + distMatrix_3
            
            #connect2 = []
            for j in range(distSum.shape[1]):
                # upper slice, use first minDist to form a triangle
                temp = np.around(distSum[:,j],decimals=6)
                index = np.where(temp == np.min(temp))
                index_last = index[0][-1]
                if j+1 == distSum.shape[1]:
                    connect.append([endPoint+j, endPoint, startPoint+index_last])
                    #connect2.append([endPoint+j, endPoint, startPoint+index_last])
                else:
                    connect.append([endPoint+j, endPoint+j+1, startPoint+index_last])
                    #connect2.append([endPoint+j, endPoint+j+1, startPoint+index_last])
            
            startPoint += slicePoint[i]
            
        self.connectivity = connect
        mesh = trimesh.Trimesh(vertices=coord,faces=connect)
        #mesh = trimesh.Trimesh(vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],faces=[[0, 1, 2]])
        #trimesh.io.export.export_mesh(mesh,saveName)
        if type(stlName)!=None:
            mesh.export(stlName)
        return coord,connect
            
    def STLvolume(self,coord=None,connect=None):
        '''
        calculate the volume of STL surface
        meshplex module for tetrahedron mesh
        '''
        centerPoint = [0,0,np.max(coord[:,2])]
        points = np.concatenate((coord, [centerPoint]))
        #coord=np.append(coord,[0,0,49])    # reshape is required, coord becomes a new (n,) array
        length = len(points)
        cells = np.concatenate((connect,np.ones((len(connect),1))*(length-1)),axis=-1).astype(np.int32)
        
        mesh = meshplex.MeshTetra(points, cells)    # all vertices must be used by the cells, or assertion error happens
        volume=mesh.cell_volumes
        volumeTotal = np.sum(volume)
        return volumeTotal
    
    def simpleStretch(self, sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, sliceStretch=None):
        '''
        simple stretch: change interval (first) and radius to achieve movement
        sliceStretch: stretch values for each slice or each point, [n, 2]
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        if type(sliceRadius)==type(None):
            sliceRadius = self.sliceRadius.copy()
        if type(sliceInterval)==type(None):
            sliceInterval = self.sliceInterval.copy()
            
        sliceRadius_new = sliceRadius.copy()
        sliceInterval_new = sliceInterval.copy()
        
        for i in range(sliceNum-1):
            sliceInterval_new[i] = sliceInterval[i]*sliceStretch[i,0]
            sliceRadius_new[i] = sliceRadius[i]*sliceStretch[i,1]
            
        coord,point = self.surfaceGeneration(sliceNum=sliceNum, slicePoint=slicePoint, sliceRadius=sliceRadius_new, sliceInterval=sliceInterval_new)
#        self.shellCoord_ES = coord.copy()
#        self.shellPoint_ES = point.copy()
        return sliceInterval_new, sliceRadius_new, coord, point
            
    def stretchMotion(self,coord=None, sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, stretchMode='lcr', sliceStretch=None):
        '''
        unfinished
        stretch motion along XYZ axis, or LCR direction
        stretch mode: 'xyz' or 'XYZ' along regular axis, 'lcr' or 'LCR' along clinical strain direction (default): stretch for each ideal slice; 'lcr_point' or 'LCR_point': stretch for each point
        sliceStretch: stretch values for each slice or each point, [n, 3]
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        if type(sliceRadius)==type(None):
            sliceRadius = self.sliceRadius.copy()
        if type(sliceInterval)==type(None):
            sliceInterval = self.sliceInterval.copy()
            
        sliceRadius_new = sliceRadius.copy()
        sliceInterval_new = sliceInterval.copy()
        
        if stretchMode=='lcr' or 'LCR':
            for i in range(sliceNum-1):
                sliceInterval_new[i] = sliceInterval[i]*sliceStretch[i,0]
            
    def strainCalculation(self,point_ED=None,point_ES=None, strain_border = False):
        '''
        calculate stretch and strain
        lengthStrain: longitudinal strain calculated by total length change
        '''
        if type(point_ED)==type(None):
            point_ED = self.shellPoint.copy()
        if type(point_ES)==type(None):
            point_ES = self.shellPoint_ES.copy()
        
        stretch_circum = point_ES[:-1,0]/point_ED[:-1,0]
        
        temp = spatial.distance_matrix(point_ED,point_ED)
        temp2 = []
        for i in range(len(point_ED)-1):
            temp2.append(temp[i,i+1])
        dist_ED = temp2.copy()
        dist_ED = np.array(dist_ED)
        
        temp = spatial.distance_matrix(point_ES,point_ES)
        temp2 = []
        for i in range(len(point_ES)-1):
            temp2.append(temp[i,i+1])
        dist_ES = temp2.copy()
        dist_ES = np.array(dist_ES)
        
        if strain_border==True:
            stretch = np.sum(dist_ES)/np.sum(dist_ED)
            strain = 0.5*(stretch**2-1)
            return stretch, strain
        
        stretch_longit = dist_ES/dist_ED
        stretch_radial = 1/stretch_longit/stretch_circum
        
        strain_longit = 0.5*(stretch_longit**2-1)
        strain_circum = 0.5*(stretch_circum**2-1)
        strain_radial = 0.5*(stretch_radial**2-1)
        
        temp = np.concatenate((stretch_longit,stretch_circum,stretch_radial))
        temp2 = temp.reshape((-1,3), order='F')
        stretch = temp2.copy()
        
        temp = np.concatenate((strain_longit,strain_circum,strain_radial))
        temp2 = temp.reshape((-1,3), order='F')
        strain = temp2.copy()
        
        print('strain calculated >_<')
        return stretch, strain
    
    def shellEndSystoleSolver(self, strain=None, sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, sliceStretch=None):
        '''
        control longit. and radial stretch to obtain specific longitudinal and circumferential strain [A,B]
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        if type(sliceRadius)==type(None):
            sliceRadius = self.sliceRadius.copy()
        if type(sliceInterval)==type(None):
            sliceInterval = self.sliceInterval.copy()
        
        strain_circum = strain[1]
        stretch_circum = (2*strain_circum+1)**0.5
        
        strain_longit = strain[0]
        stretch_longit = (2*strain_longit+1)**0.5
        sliceStretch = np.zeros((sliceNum-1,2))
        sliceStretch[:,1] = stretch_circum
        
        point_ED = self.shellPoint
        for i in range(sliceNum-1):
            # test
#            if i==27:
#                a=1
            temp1 = point_ED[i,:]
            temp2 = point_ED[i+1,:]
            
            point1 = [temp1[0],temp1[1],0]
            point2 = [temp2[0],temp2[1],temp2[2]-temp1[2]]
            temp = (point1[0]-point2[0])**2 + point2[2]**2      # distance**2
            temp2 = (stretch_longit**2)*temp-(stretch_circum*point1[0]-stretch_circum*point2[0])**2     # (z direction disctance) **2
            
            if temp2<=0:
                sliceStretch[i,0] = sliceStretch[i-1,0]
            else:
                sliceStretch[i,0] = (temp2**0.5)/abs(point2[2])
        sliceInterval_new, sliceRadius_new, coord, point = self.simpleStretch(sliceStretch=sliceStretch)
        return sliceInterval_new, sliceRadius_new, coord, point
    
    def borderEndSystoleSolver(self, surface='inner', volume=None, stretch_radial=None, error=0.005, inc=0.1, minVolume=500, stlName=None, sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None):
        '''
        achieve myocardium incompressibility by reducing the volume error
        adjust stretch_radial to achieve
        surface: inner or outer
        volume [shell_ED volume, shell_ES volume, inner_ED volume/outer_ED volume]
        minVolume (mm3): min volume of ES inner surface, avoid STL issue
        '''
        volume_myoc = abs(volume[0]-volume[2])
        if type(stretch_radial)==type(None):
            stretch_radial = self.shellStretch[:,2].copy()     # bigger stretch_radial, bigger volume_myoc
        
        volume_ES = float('inf')    # avoid accident like 0
        error_ite = (volume_myoc-abs(volume[1]-volume_ES))/volume_myoc
        while abs(error_ite)>error:
            if surface=='inner':
                index = self.innerSurface_ES(stretch_radial = stretch_radial)
                if index == 1:      # special return, stop iteration
                    return volume_ES
                coord,connect = self.STLgeneration(coord=self.innerCoord_ES, sliceNum=self.sliceNum_inner, slicePoint=self.slicePoint_inner,stlName = stlName)
                volume_ES = self.STLvolume(coord=coord,connect=connect)
                if volume_ES<=minVolume:
                    return volume_ES
            elif surface =='outer':
                self.outerSurface_ES(stretch_radial = stretch_radial)
                coord,connect = self.STLgeneration(coord=self.outerCoord_ES, sliceNum=self.sliceNum_outer, slicePoint=self.slicePoint_outer,stlName = stlName)
                volume_ES = self.STLvolume(coord=coord,connect=connect)
            else:
                print('input surface error: inner or outer!')
                sys.exit()
            error_ite = (volume_myoc-abs(volume[1]-volume_ES))/volume_myoc
            stretch_radial_backup = stretch_radial
            
            if volume_myoc-abs(volume[1]-volume_ES)>0:
                positive = True
                try:
                    if positive==positive_backup:
                        pass
                    else:
                        inc = inc*0.8
                except:
                    pass
                stretch_radial = stretch_radial*(1.0+inc)
                
            else:
                positive = False
                try:
                    if positive==positive_backup:
                        pass
                    else:
                        inc = inc*0.8
                except:
                    pass
                stretch_radial = stretch_radial*(1.0-inc)
            positive_backup = positive
        
        #stretch, strain = self.strainCalculation(point_ED= self.innerPoint, point_ES = self.innerPoint_ES)
        return volume_ES
    
    def txtStringRead(self,txtName):
        '''
        read string txt file and return float data
        '''
        try:
            with open(txtName, 'r') as f:
                data = f.readlines()
        except:
            txtName=txtName.replace('\\','/')
            #verticesFilePath = eval(repr(verticesFilePath).replace('\\', '@'))
            with open(txtName, 'r') as f:
                data = f.readlines()
        temp=[]
        for value in data:
            temp.append([])
            valueSplit=value.split('\t')
            for valueSingle in valueSplit:
                #if valueSingle.isdigit():
                try:
                    temp[-1].append(float(valueSingle))
                except:
                    pass
        return temp
    









