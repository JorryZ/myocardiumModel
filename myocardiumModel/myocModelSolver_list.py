# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 13:16:10 2019

@author: bieZY (e0348827@u.nus.edu)
Objective: ideal left ventricle myocardium ideal model
Method: 
Modules:
"""


import os
import sys
#sys.path.insert(0, "/home/yaplab/Programming/python3")
#sys.path.insert(0,"C:\\Users\\bieZY\\Documents\\Python3")
sys.path.insert(0,"C:\\Users\\bieZY\\GoogleDrive\\Work - NUS\\2.PhD Project\\2.PyCode")
sys.path.insert(0,"C:\\zhengyu\\PyCode")
import numpy as np
import shellModel as shellModel

RWT_normal = [0.3169, 0.3390, 0.3610, 0.3831, 0.4051, 0.4272]
thickness_normal = [3.5470, 3.7938, 4.0405, 4.2872, 4.5339, 4.7806]
ratio_normal = 1

RWT_HFpEF = [0.4351, 0.4700, 0.5049, 0.5398, 0.5747, 0.6096]
thickness_HFpEF = [4.8690, 5.2596, 5.6501, 6.0407, 6.4312, 6.8218]
ratio_HFpEF = 1

RWT_HFrEF = [0.3169, 0.3390, 0.3610, 0.3831, 0.4051, 0.4272]
thickness_HFrEF = [4.0220, 4.3017, 4.5815, 4.8613, 5.1410, 5.4208]
ratio_HFrEF = 1.1339

strain_longit = [-0.20, -0.18, -0.16, -0.14, -0.12, -0.10, -0.08, -0.06]
strain_circum = [-0.20, -0.18, -0.16, -0.14, -0.12, -0.10, -0.08, -0.06]

result_endoLongit = []
result_endoCircum = []
result_epiLongit = []
result_epiCircum = []
result_SV = []
result_EF = []

for s in range(len(RWT_normal)):
    RWT = RWT_normal[s]
    thickness = thickness_normal[s]
    ratio = ratio_normal
    
    result_endoLongit.append([])
    result_endoCircum.append([])
    result_epiLongit.append([])
    result_epiCircum.append([])
    result_SV.append([])
    result_EF.append([])
    
    for p in range(len(strain_longit)):
        for q in range(len(strain_circum)):
            strain_shell = [strain_longit[p], strain_circum[q]]
    
            savePath = r'D:\3.HFpEF Pig-Teresa\MRI_HFpEF_HFrEF_sham\shellModel\sham476-D42\normal-{0:.4f}_strain-{1:.2f}-{2:.2f}'.format(RWT, strain_shell[0], strain_shell[1])
            dataPath= r'D:\3.HFpEF Pig-Teresa\MRI_HFpEF_HFrEF_sham\shellModel\sham476-D42\shellModel.txt'
            thickness_inner = thickness   # normal is 2
            thickness_outer = thickness
            #strain_shell=[-0.16,-0.16]    # shell strain value: longit. , circum.
            #ratio = 1
            
            model = shellModel.shellModel()
            temp = model.txtStringRead(dataPath)
            
            sliceNum = len(temp[0])
            slicePoint= np.array(temp[1])
            temp2=[]
            for i in range(len(slicePoint)):
                temp2.append(int(slicePoint[i]))
            slicePoint=temp2.copy()
            slicePoint.reverse()
            sliceRadius=temp[0] * ratio     # mm unit
            sliceRadius.reverse()
            zAxis=temp[2]
            
            zAxis.reverse()
            sliceInterval=[]
            for i in range(len(zAxis)-1):
                sliceInterval.append(zAxis[i]-zAxis[i+1])
            
            os.makedirs(savePath, exist_ok=True)
            
            # basal dilation model
            #basalLength = int(len(sliceRadius)/3)
            #maxChange = 1.2
            #for i in range(basalLength):
            #    k = (maxChange-1)/basalLength*(basalLength-i)+1
            #    sliceRadius[i] = sliceRadius[i]*k
                
                
            #sliceRadius = np.array(sliceRadius)
            #sliceRadius[0:basalLength] = sliceRadius[0:basalLength]*1.2
            
            # from basal to apex
            #sliceNum = 6
            #slicePoint= np.array([8,7,6,4,2,1])*5
            #sliceRadius=[8,7,6,4,2,0]
            #sliceNum = 10
            #slicePoint= np.array([12,12,11,10,9,8,7,6,5,1])
            #sliceRadius=[10,10,9,8,7,6,4,3,2,0]     # mm unit
            #sliceInterval=[2,2,2,2,2,2,1.5,1.0,0.5]
            #thickness_inner = 1
            
            ###################################################
            # shell ED surface
            volume = []
            model.shellSurface(sliceNum=sliceNum,slicePoint=slicePoint,sliceRadius=sliceRadius,sliceInterval=sliceInterval)
            
            stlName = r'{:s}\shell.stl'.format(savePath)
            model.STLgeneration(coord=model.shellCoord, sliceNum=model.sliceNum, slicePoint=model.slicePoint, stlName = stlName)
            volumeTotal = model.STLvolume(coord=model.shellCoord, connect=model.connectivity)
            volume.append(volumeTotal)
            
            # shell ES surface
            interval, radius, coord, point = model.shellEndSystoleSolver(strain=strain_shell)
            model.sliceInterval_ES = interval
            model.sliceRadius_ES = radius
            model.shellCoord_ES = coord      # all the points
            model.shellPoint_ES = point
            stlName = r'{:s}\shell-ES.stl'.format(savePath)
            model.STLgeneration(coord=coord, sliceNum=model.sliceNum, slicePoint=model.slicePoint, stlName = stlName)
            volumeTotal = model.STLvolume(coord=coord, connect=model.connectivity)
            volume.append(volumeTotal)
            
            stretch, strain = model.strainCalculation()
            model.shellStretch = stretch.copy()
            model.shellStrain = strain.copy()
            meanStrain_shell=np.mean(strain[:-1,:],axis=0)
            # inner surface motion
            #temp = np.zeros(len(model.sliceRadius))
            #temp[:len(model.sliceRadius_inner)] = model.sliceRadius_inner
            #temp[len(model.sliceRadius_inner):] = model.sliceRadius[len(model.sliceRadius_inner)-1]
            #thickness_ES = stretch[:,2]*(model.sliceRadius-model.sliceRadius_inner)
            
            # inner ED surface
            model.innerSurface(thickness=thickness_inner)
            stlName = r'{:s}\inner.stl'.format(savePath)
            model.STLgeneration(coord=model.innerCoord, sliceNum=model.sliceNum_inner, slicePoint=model.slicePoint_inner, stlName = stlName)
            volumeTotal = model.STLvolume(coord=model.innerCoord, connect=model.connectivity)
            volume.append(volumeTotal)
            
            np.savetxt((r'{:s}\shellCoord.txt'.format(savePath)),model.shellCoord,fmt='%.4f',delimiter=' ')
            np.savetxt((r'{:s}\innerCoord.txt'.format(savePath)),model.innerCoord,fmt='%.4f',delimiter=' ')
            
            #sliceStretch = np.array([[0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8], [0.9,0.8]])
            #sliceStretch = np.ones((28,2))
            #sliceStretch[:,0] = sliceStretch[:,0]*0.8   # interval
            #sliceStretch[:,1] = sliceStretch[:,1]*0.8   # radius
            #interval, radius, coord, point = model.simpleStretch(sliceStretch=sliceStretch)
            
            '''
            #sliceStretch_inner = sliceStretch.copy()
            #temp = model.sliceRadius_ES-stretch[:,2]*(model.sliceRadius-model.sliceRadius_inner)
            #model.sliceRadius_inner_ES = temp
            #sliceStretch_inner[:,1] = temp/model.sliceRadius_inner
            #interval, radius, coord, point = model.simpleStretch(slicePoint=model.slicePoint_inner, sliceRadius=model.sliceRadius_inner, sliceInterval=model.sliceInterval_inner, sliceStretch=sliceStretch)
            #model.sliceRadius_inner_ES = radius
            #model.sliceInterval_inner_ES = interval
            #model.innerCoord_ES = coord
            #model.innerPoint_ES = point
            '''
            
            # inner ES surface
            #model.innerSurface_ES()
            stlName = r'{:s}\inner-ES.stl'.format(savePath)
            volumeTotal = model.borderEndSystoleSolver(surface='inner', volume=volume[0:3], stlName=stlName)
            volume.append(volumeTotal)
            #model.STLgeneration(coord=model.innerCoord_ES, sliceNum=model.sliceNum_inner, slicePoint=model.slicePoint_inner, stlName = stlName)
            stretch, strain = model.strainCalculation(point_ED= model.innerPoint, point_ES = model.innerPoint_ES)
            #meanStrain_inner=np.mean(strain[:-1,:],axis=0)
            temp = model.sliceRadius_inner[:-2].copy()
            temp2 = temp/np.sum(temp)
            
            meanStrain_inner=[]
            for i in range(3):
                meanStrain_inner.append(np.sum(temp2*strain[:-1,i]))
            
            # outer ED surface
            model.outerSurface(thickness=thickness_outer)
            stlName = r'{:s}\outer.stl'.format(savePath)
            model.STLgeneration(coord=model.outerCoord, sliceNum=model.sliceNum_outer, slicePoint=model.slicePoint_outer, stlName = stlName)
            volumeTotal = model.STLvolume(coord=model.outerCoord, connect=model.connectivity)
            volume.append(volumeTotal)
            np.savetxt((r'{:s}\outerCoord.txt'.format(savePath)),model.innerCoord,fmt='%.4f',delimiter=' ')
            
            # outer ES surface
            temp = volume[0:2].copy()
            temp.append(volume[4])
            #model.outerSurface_ES()
            stlName = r'{:s}\outer-ES.stl'.format(savePath)
            volumeTotal = model.borderEndSystoleSolver(surface='outer', volume=temp, stlName=stlName)
            volume.append(volumeTotal)
            #model.STLgeneration(coord=model.outerCoord_ES, sliceNum=model.sliceNum_outer, slicePoint=model.slicePoint_outer, stlName = stlName)
            stretch, strain = model.strainCalculation(point_ED= model.outerPoint, point_ES = model.outerPoint_ES)
            #meanStrain_outer=np.mean(strain[:-1,:],axis=0)
            temp = model.sliceRadius_outer[:-3].copy()
            temp2 = temp/np.sum(temp)
            
            meanStrain_outer=[]
            for i in range(3):
                meanStrain_outer.append(np.sum(temp2*strain[:-2,i]))
            
            
            volume.append(volume[0]-volume[2])      # myocardium volume (inner part)
            volume.append(volume[4]-volume[0])      # myocardium volume (outer part)
            SV = volume[2] - volume[3]
            EF = SV/volume[2]
            volume.append(SV)
            volume.append(EF)
            strainGlobal = np.concatenate((meanStrain_shell,meanStrain_inner,meanStrain_outer)).reshape((-1,3))
            np.savetxt((r'{:s}\volume.txt'.format(savePath)),volume,fmt='%.4f',delimiter=' ')
            np.savetxt((r'{:s}\strainGlobal.txt'.format(savePath)),strainGlobal,fmt='%.4f',delimiter=' ')
            
            result_endoLongit[-1].append(strainGlobal[1,0])
            result_endoCircum[-1].append(strainGlobal[1,1])
            result_epiLongit[-1].append(strainGlobal[2,0])
            result_epiCircum[-1].append(strainGlobal[2,1])
            result_SV[-1].append(SV)
            result_EF[-1].append(EF)
            
            print('Finished, have a good day ^_^')


