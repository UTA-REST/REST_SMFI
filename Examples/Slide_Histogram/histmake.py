import numpy as np
import glob
import os

from PIL import Image

# custom functions
from Image_Formatting import natural_sort
from Image_Formatting import PE_Vals
from Image_Formatting import Image_Converter
from Image_Formatting import FrequencyFilterFunction
from Image_Formatting import FFT_Filter

from Spot_Identification import Unique_Circle
from Trajectory_Maker import moving_average

from Cython_Def import *

from time import time
import matplotlib.pyplot as plt

import re
from glob import glob
flatten = lambda l: [item for sublist in l for item in sublist]


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

CWD = os.getcwd()
CWD = CWD+'/data/'


t1 = time()

""" DIR = [ '/Users/austinmcdonald/Data/012420/A/1/1/',
        '/Users/austinmcdonald/Data/012420/A/1/2/',
        '/Users/austinmcdonald/Data/012420/A/1/3/',
        '/Users/austinmcdonald/Data/012420/A/2/1/',
        '/Users/austinmcdonald/Data/012420/A/2/2/',
        '/Users/austinmcdonald/Data/012420/A/2/3/']
NAME = "A_FFT" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/A-glass/']
NAME = "A_FFT_glass" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/C/1/1/',
        '/Users/austinmcdonald/Data/012420/C/1/2/',
        '/Users/austinmcdonald/Data/012420/C/1/3/',
        '/Users/austinmcdonald/Data/012420/C/2/1/',
        '/Users/austinmcdonald/Data/012420/C/2/2/',
        '/Users/austinmcdonald/Data/012420/C/2/3/']
NAME = "C_FFT" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/C-glass/']
NAME = "C_FFT_glass" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/E/1/1/',
        '/Users/austinmcdonald/Data/012420/E/1/2/',
        '/Users/austinmcdonald/Data/012420/E/1/3/',
        '/Users/austinmcdonald/Data/012420/E/2/1/',
        '/Users/austinmcdonald/Data/012420/E/2/2/',
        '/Users/austinmcdonald/Data/012420/E/2/3/']
NAME = "E_FFT" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/E-glass/']
NAME = "E_FFT_glass" """



DIR = [ '/Users/austinmcdonald/Data/012420/A/1/1/',
        '/Users/austinmcdonald/Data/012420/A/1/2/',
        '/Users/austinmcdonald/Data/012420/A/1/3/',
        '/Users/austinmcdonald/Data/012420/A/2/1/',
        '/Users/austinmcdonald/Data/012420/A/2/2/',
        '/Users/austinmcdonald/Data/012420/A/2/3/']
NAME = "A_"

""" DIR = [ '/Users/austinmcdonald/Data/012420/A-glass/']
NAME = "A_glass" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/C/1/1/',
        '/Users/austinmcdonald/Data/012420/C/1/2/',
        '/Users/austinmcdonald/Data/012420/C/1/3/',
        '/Users/austinmcdonald/Data/012420/C/2/1/',
        '/Users/austinmcdonald/Data/012420/C/2/2/',
        '/Users/austinmcdonald/Data/012420/C/2/3/']
NAME = "C_" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/C-glass/']
NAME = "C_glass" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/E/1/1/',
        '/Users/austinmcdonald/Data/012420/E/1/2/',
        '/Users/austinmcdonald/Data/012420/E/1/3/',
        '/Users/austinmcdonald/Data/012420/E/2/1/',
        '/Users/austinmcdonald/Data/012420/E/2/2/',
        '/Users/austinmcdonald/Data/012420/E/2/3/']
NAME = "E_" """

""" DIR = [ '/Users/austinmcdonald/Data/012420/E-glass/']
NAME = "E_glass" """


datamin  = 0
datamax  = 1000
numbins  = 1000
bins     = np.linspace(datamin, datamax, numbins)
HIST     = np.zeros(numbins-1, dtype='int32')

for x in range(0,len(DIR)):
    
    #Dir  = CWD+DIR[x]
    Dir = DIR[x]
    print("looking at "+Dir)
    end = 350
    DataFiles = glob.glob(Dir+'/*.tif')
    DataFiles = natural_sort(DataFiles)
    DataFiles = DataFiles[0:end]

    Row    = int(250)
    Col    = Row
    Yindex = int(256)
    Xindex = int(256)

    eOffset, eCoeff = PE_Vals(DataFiles[0])
    Test = Image_Converter(DataFiles[0], eOffset, eCoeff, Xindex, Yindex, Row, Col)
    Shape = Test.shape[0]
    radius = 130
    Xoff = 250
    Yoff = 220
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

np.save(NAME,HIST)
print("finished "+NAME)
t2 = time()
print(t2-t1)

