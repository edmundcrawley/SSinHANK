"""
 Creates all tables and figures for Chapter 3 of my dissertation
"""
import sys 
import os

base_dir = os.getcwd()
dolo_dir = os.getcwd() + '/Python/dolocode/'
HANK_dir = os.getcwd() + '/Python/'

print('First run the TANK models in Dolo')
sys.path.append(dolo_dir)
os.chdir(dolo_dir)
import TANKdolo_MAIN


print('Next run HANK model')
sys.path.append(HANK_dir)
os.chdir(HANK_dir)
import Main


    
    
