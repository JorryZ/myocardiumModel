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
  Author: jorry.zhengyu@gmail.com         29July2020           -V1.0.0 Created, test version, horizontal radial stretch
  Author: jorry.zhengyu@gmail.com         26AUGU2020           -V1.0.1, release version, improve innerSurface_ES function to solve inner_ES STL volume error, point_ES error
  Author: jorry.zhengyu@gmail.com         29AUGU2020           -V1.0.2, release version, add thickRatio for inner and outer ED surface for thickness ajustment, different strain calculation
  Author: jorry.zhengyu@gmail.com         05Sept2020           -V2.0.0 Algorithm improvement, test version, radial direction using shell surface normal
  Author: jorry.zhengyu@gmail.com         11Sept2020           -V2.0.1 inner ES smoothing improvement
"""
print('shellModel release version 2.0.1')

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

def normalize_v2(a):
    #a = np.array(a)
    b = np.zeros(np.shape(a))
    for i in range(len(b)):
        b[i,:] = a[i,:]/np.linalg.norm(a[i,:])
    return b

def perpendicular( a ) :
    b = np.zeros(np.shape(a))
    for i in range(len(b)):
        b[i,0] = a[i,1]
        b[i,1] = -a[i,0]
    return b

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
        self.radialDirection_inner = None
        
        self.thickness_inner_ES = None
#        self.sliceNum_inner_ES = None
#        self.slicePoint_inner_ES = None
        self.sliceRadius_inner_ES = None
        self.sliceInterval_inner_ES = None
        self.innerCoord_ES = None
        self.innerPoint_ES = None
        self.index = 0
        
        # outer part
        self.thickness_outer = None
        self.sliceNum_outer = None
        self.slicePoint_outer = None
        self.sliceRadius_outer = None
        self.sliceInterval_outer = None
        self.outerCoord = None
        self.outerPoint = None
        self.connectivity_outer = None
        self.radialDirection_outer = None
        
        self.sliceRadius_outer_ES = None
        self.sliceInterval_outer_ES = None
        self.outerCoord_ES = None
        self.outerPoint_ES = None
        
        self.connectivity = None
        
        print('shellModel_init_ done')
    
    def radialDirection(self, shellPoint = None):
        '''
        calculate radial direction based on the longitudinal direction
        shellPoint [n,3], but y coordinates are 0
        '''
        longitDirection = np.zeros((len(shellPoint)-1,2))
        longitDirection[:,0] = shellPoint[1:,0] - shellPoint[:-1,0]
        longitDirection[:,1] = shellPoint[1:,2] - shellPoint[:-1,2]
        temp = normalize_v2(longitDirection)
        radialDirection = perpendicular(temp)
        return radialDirection
    
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
        
    
    def simpleStretch(self, sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None, sliceStretch=None):
        '''
        simple stretch: change interval (first) and radius to achieve movement, shell ES surface
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
        
        
    def innerSurface(self,sliceNum=None, slicePoint=None, thickness=None, thickRatio=1.):
        '''
        use thickness to obtain inner ED surface
        thickRatio = thick_apex/thick_basal
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        
        if type(thickness) == int or type(thickness) == float and self.sliceNum>1:
#        if len(thickness) == 1 and self.sliceNum>1:
            temp = np.ones(self.sliceNum)*thickness
            self.thickness_inner = temp
        else:
            self.thickness_inner = thickness
        if thickRatio > 1:
            pass
        elif abs(thickRatio-1.)!=0.:
            decrement = (thickRatio-1)/(sliceNum-1)
            ratio = np.arange(1.,thickRatio+decrement,decrement)
            self.thickness_inner = self.thickness_inner*ratio
        
        shellPoint = self.shellPoint.copy()
        radialDirection = self.radialDirection(shellPoint)
        innerPoint = np.zeros(np.shape(shellPoint))
        innerPoint[:-1,0] = self.thickness_inner[:-1]*radialDirection[:,0] + shellPoint[:-1,0]
        innerPoint[:-1,2] = self.thickness_inner[:-1]*radialDirection[:,1] + shellPoint[:-1,2]
        innerPoint[-1,2] = self.thickness_inner[-1] + shellPoint[-1,2]
        
        self.sliceNum_inner = sliceNum
        self.slicePoint_inner = slicePoint.copy()
        self.sliceRadius_inner = np.array(innerPoint[:,0]).copy()
        self.radialDirection_inner = radialDirection.copy()
        
        z = innerPoint[-1,2]
        sliceInterval_inner = innerPoint[:-1,2] - innerPoint[1:,2]
        self.sliceInterval_inner = sliceInterval_inner.copy()
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum_inner, slicePoint=self.slicePoint_inner, sliceRadius=self.sliceRadius_inner, sliceInterval=self.sliceInterval_inner)
        coord[:,2] = coord[:,2] + z
        point[:,2] = point[:,2] + z
        self.innerCoord = coord.copy()
        self.innerPoint = point.copy()
        
        print('innerSurface ED completed >_<')
    
    def innerSurface_ES(self,sliceNum=None, slicePoint=None, stretch_radial=None, mode='smooth', radialDir = 'ES', smoothInit = [5,0], ratio=0.9):
        '''
        simplified inner surface motion, use radial stretch to calculate radius (not same direction, especially apical part)
        # mode: 'shape' to maintain the ES inner surface length, 'smooth' to obtain smoother ES inner STL to improve apex shape
        radialDir: radial direction use 'ED' middle surface normal or 'ES' middle surface normal
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum_inner
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint_inner.copy()
        if type(stretch_radial)==type(None):
            stretch_radial = self.shellStretch[:,2].copy()
        
        thickness = self.thickness_inner.copy()
        thickness[:-1] = stretch_radial*self.thickness_inner[:-1]
        thickness[-1] = stretch_radial[-1]*self.thickness_inner[-1]
        self.thickness_inner_ES = thickness.copy()
        
        shellPoint = self.shellPoint_ES.copy()
        
        if radialDir == 'ED':
            radialDirection = self.radialDirection_inner.copy()     # shell ED surface normal
        elif radialDir == 'ES':
            radialDirection = self.radialDirection(shellPoint)     # shell ES surface normal
        else:
            print('Error: radialDIR input "ED" or "ES"')
            sys.exit()
        innerPoint = np.zeros(np.shape(shellPoint))
        innerPoint[:-1,0] = self.thickness_inner_ES[:-1]*radialDirection[:,0] + shellPoint[:-1,0]
        innerPoint[:-1,2] = self.thickness_inner_ES[:-1]*radialDirection[:,1] + shellPoint[:-1,2]
        innerPoint[-1,2] = self.thickness_inner_ES[-1] + shellPoint[-1,2]
        
        if innerPoint[0,0] <= 0:
            flag = 5
            return flag
        
        # innerPoint X smoothing for bad condition
        flag = 0
        flag_X = 1
        while(flag_X!=0):
            sliceBad = []
            temp = innerPoint[smoothInit[0]:-1,0] - innerPoint[(smoothInit[0]+1):,0]
            for i in range(len(temp)):
                if temp[i]<=0:
                    sliceBad.append(smoothInit[0]+1+i)
    #                sliceBad = i+3
                    flag_X = 1
                    flag = 1
            for i in range(len(innerPoint)-1):
                if innerPoint[i,0]<=0:
                    sliceBad = []
                    sliceBad.append(i)
                    flag_X = 2
                    flag = 2
                    break
            if len(sliceBad)==0:
                flag_X = 0
            if flag_X == 1:
                sliceBad = np.flip(np.array(sliceBad))
                for i in range(len(sliceBad)):
                    if sliceBad[i]< len(innerPoint)-2:
                        innerPoint[sliceBad[i],0] = (innerPoint[sliceBad[i]+2,0]+innerPoint[sliceBad[i]-1,0]+innerPoint[sliceBad[i],0]+innerPoint[sliceBad[i]+1,0]+innerPoint[sliceBad[i]+2,0])/5
                    elif sliceBad[i]== len(innerPoint)-2:
                        innerPoint[sliceBad[i],0] = (innerPoint[sliceBad[i]-1,0]+innerPoint[sliceBad[i],0]+innerPoint[sliceBad[i]+1,0])/3
            elif flag_X==2:
                for i in range(sliceBad[0], len(innerPoint)):
                    innerPoint[i,0] = self.innerPoint[i,0] * innerPoint[sliceBad[0]-1,0]/self.innerPoint[sliceBad[0]-1,0]
        '''
        # innerPoint Z smoothing for bad condition
        flag_Z = 1
        while(flag_Z!=0):
            sliceBad = []
            temp = innerPoint[smoothInit[1]:-1,2] - innerPoint[(smoothInit[1]+1):,2]
            for i in range(len(temp)):
                if temp[i]<=0:
                    sliceBad.append(smoothInit[1]+1+i)
                    flag_Z = 1
                    if flag!= 0:
                        flag = 4
                    else:
                        flag = 3
            if len(sliceBad)==0:
                flag_Z = 0
            if flag_Z == 1:
                for i in range(len(sliceBad)):
                    if sliceBad[i]==len(innerPoint)-1:
                        innerPoint[sliceBad[i]-1,2] = innerPoint[sliceBad[i]-2,2]*1/3+innerPoint[sliceBad[i],2]*2/3
                    else:
                        innerPoint[sliceBad[i],2] = innerPoint[sliceBad[i]-1,2]*1/3+innerPoint[sliceBad[i]+1,2]*2/3
        '''
        
        '''
        flag = 0    # condition of innerPoint change, 0 no change, 1 only change X, 2 only change Z, 3 change both X and Z
        
        # innerPoint X ajustment for bad condition
        sliceBad = len(innerPoint)
        
        temp = innerPoint[smoothInit[0]:-1,0] - innerPoint[smoothInit[0]+1:,0]
        for i in range(len(temp)):
            #if temp[i]<=0 and abs(temp[i]/innerPoint[i,0])>0.01:
            if temp[i]<=0:
                sliceBad = np.min([sliceBad, smoothInit[0]+i])
#                sliceBad = i+3
                flag = 1
                break
        
        for i in range(len(innerPoint)-1):
            if innerPoint[i,0]<=0:
                sliceBad = i
                flag = 2
                break

        
        if sliceBad < len(innerPoint):
            for i in range(sliceBad-1, len(innerPoint)):
                innerPoint[i,0] = self.innerPoint[i,0] * innerPoint[sliceBad-1,0]/self.innerPoint[sliceBad-1,0]
#                temp = 0.5
#                temp2 = innerPoint[i,0]
#                while(innerPoint[i,0]>=innerPoint[i-1,0]):
#                #innerPoint[i,0] = self.shellPoint_ES[i,0] * innerPoint[sliceBad-1,0]/self.shellPoint_ES[sliceBad-1,0]
#                    temp = temp*ratio
#                    innerPoint[i,0] = (self.innerPoint[i,0] * innerPoint[sliceBad-1,0]/self.innerPoint[sliceBad-1,0])*(1-temp) + temp2*temp
        '''
        
        # innerPoint Z ajustment for bad condition
        # test
        sliceBad2 = len(innerPoint)
        temp = innerPoint[smoothInit[1]:-1,2] - innerPoint[(smoothInit[1]+1):,2]
        for i in range(len(temp)):
            if temp[i]<=0:
                sliceBad2 = np.min([sliceBad2, smoothInit[1]+i])
                if flag == 0:
                    flag = 3
                else:
                    flag = 4
#                sliceBad = i+3
                break
        if sliceBad2 < len(innerPoint):
            # adustment 1, fixed apex point
            for i in range(sliceBad2-1, len(innerPoint)-1):
                dist = (innerPoint[i-1,2] - innerPoint[-1,2])*(self.innerPoint[i,2] - self.innerPoint[-1,2])/(self.innerPoint[i-1,2] - self.innerPoint[-1,2])
                innerPoint[i,2] = innerPoint[-1,2] + dist
#                temp = 0.5
#                temp2 = innerPoint[i,2] - innerPoint[-1,2]
#                while(innerPoint[i,0]>=innerPoint[i-1,0]):
#                    temp = temp*ratio
#                    #dist = (innerPoint[i-1,2] - innerPoint[-1,2])*(self.shellPoint_ES[i,2] - self.shellPoint_ES[-1,2])/(self.shellPoint_ES[i-1,2] - self.shellPoint_ES[-1,2])
#                    dist = ((innerPoint[i-1,2] - innerPoint[-1,2])*(self.innerPoint[i,2] - self.innerPoint[-1,2])/(self.innerPoint[i-1,2] - self.innerPoint[-1,2]))*(1-temp) + temp2*temp
#                    innerPoint[i,2] = innerPoint[-1,2] + dist
            # adustment 2
#            for i in range(sliceBad2-1, len(innerPoint)):
#                #temp = (innerPoint[i-2,2] - innerPoint[i-1,2])*(self.shellPoint_ES[i-1,2] - self.shellPoint_ES[i,2])/(self.shellPoint_ES[i-2,2] - self.shellPoint_ES[i-1,2])
#                temp = (innerPoint[i-2,2] - innerPoint[i-1,2])*(self.innerPoint[i-1,2] - self.innerPoint[i,2])/(self.innerPoint[i-2,2] - self.innerPoint[i-1,2])
#                innerPoint[i,2] = innerPoint[i-1,2] - temp
        
        
        self.sliceRadius_inner_ES = np.array(innerPoint[:,0]).copy()
        z = innerPoint[-1,2]
        sliceInterval_inner = innerPoint[:-1,2] - innerPoint[1:,2]
        self.sliceInterval_inner_ES = sliceInterval_inner.copy()
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum_inner, slicePoint=self.slicePoint_inner, sliceRadius=self.sliceRadius_inner_ES, sliceInterval=self.sliceInterval_inner_ES)
        coord[:,2] = coord[:,2] + z
        point[:,2] = point[:,2] + z
        self.innerCoord_ES = coord.copy()
        self.innerPoint_ES = point.copy()
        self.flag = flag
        
        return flag    # condition return
    
    def outerSurface(self,sliceNum=None, slicePoint=None, thickness=None, thickRatio=1.):
        '''
        use thickness to obtain outer ED surface
        thickRatio = thick_apex/thick_basal
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint.copy()
        
        if type(thickness) == int or type(thickness) == float and self.sliceNum>1:
#        if len(thickness) == 1 and self.sliceNum>1:
            temp = np.ones(self.sliceNum)*thickness
            self.thickness_outer = temp
        else:
            self.thickness_outer = thickness
        if thickRatio > 1:
            pass
        elif abs(thickRatio-1.)!=0.:
            decrement = (thickRatio-1)/(sliceNum-1)
            ratio = np.arange(1.,thickRatio+decrement,decrement)
            self.thickness_outer = self.thickness_outer*ratio
        
        shellPoint = self.shellPoint.copy()
        radialDirection = self.radialDirection(shellPoint)
        radialDirection = radialDirection * (-1)
        outerPoint = np.zeros(np.shape(shellPoint))
        outerPoint[:-1,0] = self.thickness_outer[:-1]*radialDirection[:,0] + shellPoint[:-1,0]
        outerPoint[:-1,2] = self.thickness_outer[:-1]*radialDirection[:,1] + shellPoint[:-1,2]
        outerPoint[-1,2] = -1*self.thickness_outer[-1] + shellPoint[-1,2]
        
        self.sliceNum_outer = sliceNum
        self.slicePoint_outer = slicePoint.copy()
        self.sliceRadius_outer = np.array(outerPoint[:,0]).copy()
        self.radialDirection_outer = radialDirection.copy()
        
        z = outerPoint[-1,2]
        sliceInterval_outer = outerPoint[:-1,2] - outerPoint[1:,2]
        self.sliceInterval_outer = sliceInterval_outer.copy()
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum_outer, slicePoint=self.slicePoint_outer, sliceRadius=self.sliceRadius_outer, sliceInterval=self.sliceInterval_outer)
        coord[:,2] = coord[:,2] + z
        point[:,2] = point[:,2] + z
        self.outerCoord = coord.copy()
        self.outerPoint = point.copy()
        
        print('outerSurface ED completed >_<')
    
    def outerSurface_ES(self,sliceNum=None, slicePoint=None, stretch_radial=None, mode='smooth', radialDir = 'ES'):
        '''
        simplified outer surface motion, use radial stretch to calculate radius (not same direction, especially apical part)
        mode: 'shape' to maintain the ES outer surface length, 'smooth' to obtain smoother ES outer STL
        to improve apex shape
        '''
        if type(sliceNum)==type(None):
            sliceNum = self.sliceNum_inner
        if type(slicePoint)==type(None):
            slicePoint = self.slicePoint_outer.copy()
        if type(stretch_radial)==type(None):
            stretch_radial = self.shellStretch[:,2].copy()
        
        thickness = self.thickness_outer.copy()
        thickness[:-1] = stretch_radial*self.thickness_outer[:-1]
        thickness[-1] = stretch_radial[-1]*self.thickness_outer[-1]
        self.thickness_outer_ES = thickness.copy()
        
        
        shellPoint = self.shellPoint_ES.copy()
        #radialDirection = self.radialDirection(shellPoint)
        #radialDirection = radialDirection *(-1)
        if radialDir == 'ED':
            radialDirection = self.radialDirection_outer.copy()     # shell ED surface normal
        elif radialDir == 'ES':
            radialDirection = self.radialDirection(shellPoint)     # shell ES surface normal
            radialDirection = radialDirection *(-1)
        else:
            print('Error: radialDIR input "ED" or "ES"')
            sys.exit()
        
        outerPoint = np.zeros(np.shape(shellPoint))
        outerPoint[:-1,0] = self.thickness_outer_ES[:-1]*radialDirection[:,0] + shellPoint[:-1,0]
        outerPoint[:-1,2] = self.thickness_outer_ES[:-1]*radialDirection[:,1] + shellPoint[:-1,2]
        outerPoint[-1,2] = -1*self.thickness_outer_ES[-1] + shellPoint[-1,2]
        
        self.sliceRadius_outer_ES = np.array(outerPoint[:,0]).copy()
        
        z = outerPoint[-1,2]
        sliceInterval_outer = outerPoint[:-1,2] - outerPoint[1:,2]
        self.sliceInterval_outer_ES = sliceInterval_outer.copy()
        
        coord,point = self.surfaceGeneration(sliceNum=self.sliceNum_outer, slicePoint=self.slicePoint_outer, sliceRadius=self.sliceRadius_outer_ES, sliceInterval=self.sliceInterval_outer_ES)
        coord[:,2] = coord[:,2] + z
        point[:,2] = point[:,2] + z
        self.outerCoord_ES = coord.copy()
        self.outerPoint_ES = point.copy()
        
        return 0    # normal return
    
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
            
    def strainCalculation(self,point_ED=None,point_ES=None, strain_type = 'green', strain_border = False):
        '''
        calculate stretch and strain
        strain_type: 'green' strain or 'engineering' strain
        strain_border: longitudinal strain calculated by total length if True
        
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
            if strain_type == 'green':
                strain = 0.5*(stretch**2-1)
            elif strain_type == 'engineering':
                strain = stretch -1
            return stretch, strain
        
        stretch_longit = dist_ES/dist_ED
        stretch_radial = 1/stretch_longit/stretch_circum
        
        if strain_type == 'green':
            strain_longit = 0.5*(stretch_longit**2-1)
            strain_circum = 0.5*(stretch_circum**2-1)
            strain_radial = 0.5*(stretch_radial**2-1)
        elif strain_type == 'engineering':
            strain_longit = stretch_longit-1
            strain_circum = stretch_circum-1
            strain_radial = stretch_radial-1
        else:
            print('Please input green or engineering for strain_type')
        
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
                # apical part may have bad motion
                sliceStretch[i,0] = sliceStretch[i-1,0]
            else:
                sliceStretch[i,0] = (temp2**0.5)/abs(point2[2])
        sliceInterval_new, sliceRadius_new, coord, point = self.simpleStretch(sliceStretch=sliceStretch)
        return sliceInterval_new, sliceRadius_new, coord, point
    
    def borderEndSystoleSolver(self, surface='inner', volume=None, stretch_radial=None, error=0.005, inc=0.1, minDiff=0.0001, minVolume=200, mode='smooth', radialDir = 'ES', stlName=None, sliceNum=None, slicePoint=None, sliceRadius=None, sliceInterval=None):
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
        
        volume_ES = float('-inf')    # avoid accident like 0
        error_ite = (volume_myoc-abs(volume[1]-volume_ES))/volume_myoc
        stretch_radial_backup1 = stretch_radial.copy()
        Incre_max = stretch_radial.copy()*0
        Decre_min = stretch_radial.copy()*100
        
        while abs(error_ite)>error:
            if surface=='inner':
                flag = self.innerSurface_ES(stretch_radial = stretch_radial, mode=mode, radialDir=radialDir)
                if flag == 5:
                    if positive == True:
                        stretch_radial = stretch_radial/(1.0+inc)
                        inc = inc*0.9
                        stretch_radial = stretch_radial*(1.0+inc)
                        continue
#                index = self.innerSurface_ES(stretch_radial = stretch_radial, mode=mode)
#                if index == 1:      # special return, stop iteration
#                    return volume_ES
                coord,connect = self.STLgeneration(coord=self.innerCoord_ES, sliceNum=self.sliceNum_inner, slicePoint=self.slicePoint_inner,stlName = stlName)
                volume_ES = self.STLvolume(coord=coord,connect=connect)
                if volume_ES<=minVolume and error_ite>0:
                #if volume_ES<=minVolume:
                    return volume_ES
            elif surface =='outer':
                self.outerSurface_ES(stretch_radial = stretch_radial, radialDir=radialDir)
                coord,connect = self.STLgeneration(coord=self.outerCoord_ES, sliceNum=self.sliceNum_outer, slicePoint=self.slicePoint_outer,stlName = stlName)
                volume_ES = self.STLvolume(coord=coord,connect=connect)
            else:
                print('input surface error: inner or outer!')
                sys.exit()
            
            error_ite = (volume_myoc-abs(volume[1]-volume_ES))/volume_myoc
            #stretch_radial_backup = stretch_radial
            #print('error_ite: ',error_ite, '; stretch_radial:',stretch_radial[0], '; Incre:',Incre_max[0])
            
            # interation method 1
            if volume_myoc-abs(volume[1]-volume_ES)>0:
                positive = True
                if stretch_radial[0] > Incre_max[0]:
                    Incre_max = stretch_radial.copy()
                try:
                    if positive==positive_backup:
                        #inc = inc*0.8
                        stretch_radial = stretch_radial*(1.0+inc)
                        #pass
                    else:
                        inc = inc*0.9
                        stretch_radial = (Incre_max+Decre_min)/2
                        #stretch_radial = (stretch_radial_backup2+stretch_radial_backup1)/2
                except:
                    stretch_radial = stretch_radial*(1.0+inc)
                    #pass
                
            else:
                positive = False
                if stretch_radial[0] < Decre_min[0]:
                    Decre_min = stretch_radial.copy()
                try:
                    if positive==positive_backup:
                        #inc = inc*0.8
                        stretch_radial = stretch_radial*(1.0-inc*0.1)
                        #pass
                    else:
                        inc = inc*0.9
                        stretch_radial = (Incre_max+Decre_min)/2
                        #stretch_radial = (stretch_radial_backup2+stretch_radial_backup1)/2
                except:
                    stretch_radial = stretch_radial*(1.0-inc)
                    #pass
            if abs(Decre_min[0]-Incre_max[0])<minDiff:      # due to surface smoothing, volume change is not continuous sometime
                return volume_ES
            positive_backup = positive
            #stretch_radial_backup2 = stretch_radial_backup1.copy()
            #stretch_radial_backup1 = stretch_radial.copy()
            
            '''
            # interation method 2
            if volume_myoc-abs(volume[1]-volume_ES)>0:
                positive = True
                try:
                    if positive==positive_backup:
                        #inc = inc*0.8
                        stretch_radial = stretch_radial*(1.0+inc*0.1)
                        #pass
                    else:
                        inc = inc*0.6
                        stretch_radial = stretch_radial*(1.0+inc)
                except:
                    stretch_radial = stretch_radial*(1.0+inc)
                    #pass
                
            else:
                positive = False
                try:
                    if positive==positive_backup:
                        #inc = inc*0.8
                        stretch_radial = stretch_radial*(1.0-inc*0.1)
                        #pass
                    else:
                        inc = inc*0.6
                        stretch_radial = stretch_radial*(1.0-inc)
                except:
                    stretch_radial = stretch_radial*(1.0-inc)
                    #pass
                
            positive_backup = positive
            '''
            
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
    
    









