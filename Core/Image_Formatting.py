# -----------------------------------------------------------------------------
#  REST_SMFI | Image_Formatting.py
#
#  Various functions for data analysis
#   * Author: Austin McDonald
#   * Creation date: Oct 2019
# -----------------------------------------------------------------------------

from PIL import Image
import numpy as np
import re

# Nedded for background fit
#import numpy.polynomial.polynomial as poly
# needed for spotfinder
#import scipy.ndimage as ndimage
#import scipy.ndimage.filters as filters

#from scipy import optimize


#######################################
# Sorts the data files even with text in the name. 
#######################################
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)



#######################################
# The PE vals looks in a raw tiff file, it readin all of the symbols and
# splits where the info about the camara settings are
# it returns 2 values that are used in converting the adc counts to PE
#######################################
def PE_Vals(FILE):
    with open(FILE, encoding="utf8", errors='ignore') as f:
        contents = f.read()
    contents = contents.split("Created by Hamamatsu Inc.")[1].split('\n')
    for line in contents:
        if 'eOffset1' in line:
            eOffset = float(line.split("eOffset1 = ")[1])
        elif 'eCoeff1' in line:
            eCoeff = float(line.split("eCoeff1 = ")[1])
    return eOffset, eCoeff


#######################################
# This takes the raw tiff file and crops it to the proper shape and converts
# the adc counts into PE
#######################################
def Image_Converter(FILE,eOffset,eCoeff,Xindex,Yindex,Row,Col):
    TestImage = Image.open(FILE)
    Testspot = np.array(TestImage)[Yindex-Row:Yindex+Row+1,Xindex-Col:Xindex+Col+1]
    Shape = Testspot.shape[0]
    eCoeffM = eCoeff*np.ones(Shape**2).reshape((Shape,Shape))
    eOffsetM = eOffset*np.ones(Shape**2).reshape((Shape,Shape))
    Testspot = eCoeffM*(Testspot - eOffsetM)
    TestImage.close()
    return Testspot



#######################################
# Used to filter the high frequency part of the images.
#######################################
def FrequencyFilterFunction(Shape,FreqCut,FreqCutWidth):
    # the ben classic tanh step function to give
    #  a smooth transition in k space
    Filter=lambda x: (0.5+np.tanh((x-FreqCut)/FreqCutWidth)/2.)
    Freqs = np.fft.fftfreq(Shape)
    FilterArray=np.zeros([Shape,Shape])
    for i in range(0,Shape):
        for j in range(0,Shape):
            Freq2D=(Freqs[i]**2+Freqs[j]**2)**0.5
            FilterArray[i,j]=Filter(Freq2D)
    return FilterArray

def FFT_Filter(ImageArray, FilterArray):
    FFTed = np.fft.fft2(ImageArray)
    FilteredSlideFFT = FFTed*FilterArray
    return np.abs(np.fft.ifft2(FilteredSlideFFT))







