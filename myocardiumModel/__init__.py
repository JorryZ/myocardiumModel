# -*- coding: utf-8 -*-
"""
Created on Wed July 29 13:16:10 2020
@author: JorryZ
File: __init__.py
Description: load shellModel function
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: jorry.zhengyu@gmail.com         29July2020             -V1.0.0 test version
                                                        -myocardiumModel version 1.0.0
                                                        
Requirements:
    numpy
    scipy
    trimesh
    meshplex
All rights reserved.
"""
_version='1.0.0'
print('myocardiumModel version',_version)

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import myocardiumModel.shellModel as shellModel                                                        
