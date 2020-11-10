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
                                                        -shellModel version 1.0.0
  Author: jorry.zhengyu@gmail.com         26AUGU2020             -V1.0.1 release version
                                                        -shellModel version 1.0.1 
  Author: jorry.zhengyu@gmail.com         29AUGU2020             -V1.0.2 release version
                                                        -shellModel version 1.0.2  
  Author: jorry.zhengyu@gmail.com         05Sept2020             -V2.0.0 test version
                                                        -shellModel version 2.0.0        
  Author: jorry.zhengyu@gmail.com         11Sept2020             -V2.0.1 test version
                                                        -shellModel version 2.0.1      
  Author: jorry.zhengyu@gmail.com         16Sept2020             -V2.0.2 release version
                                                        -shellModel version 2.0.2
  Author: jorry.zhengyu@gmail.com         10Nov2020              -V2.0.3 release version
                                                        -shellModel version 2.0.3
         
Requirements:
    numpy
    scipy
    trimesh
    meshplex
All rights reserved.
"""
_version='2.0.3'
print('myocardiumModel version',_version)

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import myocardiumModel.shellModel as shellModel                                                        
