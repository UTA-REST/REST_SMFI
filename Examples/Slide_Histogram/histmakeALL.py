import numpy as np
import glob
import os

from PIL import Image

# custom functions
from DefDefDefinitions import PE_Vals
from DefDefDefinitions import Image_Converter
from DefDefDefinitions import FrequencyFilterFunction
from DefDefDefinitions import FFT_Filter
from DefDefDefinitions import Unique_Circle
from DefDefDefinitions import Fit_2D_Gaussian
from DefDefDefinitions import gaussian_2d
from DefDefDefinitions import moving_average

from C_DefDef import *

from time import time
import matplotlib.pyplot as plt

import re
from glob import glob
flatten = lambda l: [item for sublist in l for item in sublist]

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

# save the data here
DataFolder = "HistData/"
t1 = time()
DIR_FULL  = []
NAME_FULL = []

Head_Dir    = "/Users/austinmcdonald/Data/022520/PVA/"
# finding the sample folders under head
Sample_Dirs = glob(Head_Dir+"*/")
Slide_Dirs  = []
Spot_Dirs   = []
# finding the slide folders under sample
for Sample in Sample_Dirs:
    Slide_Dirs.append(glob(Sample+"*/"))
Slide_Dirs = flatten(Slide_Dirs)
# finding the spots folders under slide
for Spot in Slide_Dirs:
    Spot_Dirs.append(glob(Spot+"*/"))
Spot_Dirs = flatten(Spot_Dirs)
Spot_Dirs = natural_sort(Spot_Dirs)

# making names for the saved files
Data_Name = []
for data in Spot_Dirs:
    ll = data.split(Head_Dir)[1].split('/')
    name = ll[0]+"_"
    for l in ll[1::]:
        name +=l
    Data_Name.append(name)
print("************************************ ")
print("Found the following directories ")
print("************************************ ")

for spot in Spot_Dirs:
    print(spot)
print("************************************ ")
print("Starting the loop ")
print("************************************ ")
print( )


for q in range(0,len(Data_Name)):
    NAME = Data_Name[q]
    DIR  = Spot_Dirs[q]

    datamin  = 0
    datamax  = 3000
    numbins  = 500
    bins     = np.linspace(datamin, datamax, numbins)
    HIST     = np.zeros(numbins-1, dtype='int32')

    print("looking at "+DIR)
    end = 300
    DataFiles = glob(DIR+'*.tif')
    DataFiles = natural_sort(DataFiles)
    DataFiles = DataFiles[0:end]
    print("Number of files "+str(len(DataFiles)))

    Row    = int(250)
    Col    = Row
    Yindex = int(256)
    Xindex = int(256)

    eOffset, eCoeff = PE_Vals(DataFiles[0])
    Test = Image_Converter(DataFiles[0], eOffset, eCoeff, Xindex, Yindex, Row, Col)
    Shape = Test.shape[0]
    radius = 160
    Xoff = 260
    Yoff = 250
    X_circle, Y_circle = Unique_Circle(radius, Xoff, Yoff)
    keep = np.where(Y_circle<Shape)
    Y_circle = Y_circle[keep]
    X_circle = X_circle[keep]
    keep = np.where(X_circle<Shape)
    Y_circle = Y_circle[keep]
    X_circle = X_circle[keep]

    #FreqCut=0.03
    #FreqCutWidth=0.04
    #FilterArray = FrequencyFilterFunction(Shape,FreqCut,FreqCutWidth)


    for img in DataFiles:
        HIST_VAL = []
        eOffset, eCoeff = PE_Vals(img)
        ReducedImage = Image_Converter(img, eOffset, eCoeff, Xindex, Yindex, Row, Col)
        #ReducedImage = FFT_Filter(ReducedImage, FilterArray)
        r_vals = ReducedImage[X_circle,Y_circle]

        r_cut  = ((X_circle-Xoff)**2+(Y_circle-Yoff)**2).max()
        HIST_VAL = Hist_maker(r_cut, Yoff, Xoff, ReducedImage.shape[0], ReducedImage)
        htemp, jnk = np.histogram(HIST_VAL, bins)
        HIST += htemp

    np.save(DataFolder + NAME,HIST)
    print("finished "+NAME)
t2 = time()
print(t2-t1)

