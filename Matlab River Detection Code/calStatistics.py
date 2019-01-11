import arcpy
import os
from arcpy import env
from arcpy.sa import *
import csv

env.overwriteOutput=True
env.workspace=r'D:\1. Codes and Programs\1. river detection comparison methods\00_my_Matlab_codes_for_single_image\sentinel2_T22WEV'

rasters=arcpy.ListRasters("*","tif")
for raster in rasters:
    print (raster)
    arcpy.CalculateStatistics_management(raster)
    arcpy.BuildPyramids_management(raster)
    


            


