# -----------------------------------------------------------------------------
#  REST_SMFI | Definitions_Dev.py
#
#  Various functions for data analysis
#   * Author: Austin McDonald
#   * Creation date: Oct 2019
# -----------------------------------------------------------------------------

from PIL import Image
import numpy as np
import re

# Nedded for background fit
import numpy.polynomial.polynomial as poly
# needed for spotfinder
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

from scipy import optimize


#######################################
# Sorts the data files even with text in the name. 
#######################################
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


#############
# The PE vals looks in a raw tiff file, it readin all of the symbols and
# splits where the info about the camara settings are
# it returns 2 values that are used in converting the adc counts to PE
#############
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

#############
# This takes the raw tiff file and crops it to the proper shape and converts
# the adc counts into PE
#############
def Image_Converter(FILE,eOffset,eCoeff,Xindex,Yindex,Row,Col):
    TestImage = Image.open(FILE)
    Testspot = np.array(TestImage)[Yindex-Row:Yindex+Row+1,Xindex-Col:Xindex+Col+1]
    Shape = Testspot.shape[0]
    eCoeffM = eCoeff*np.ones(Shape**2).reshape((Shape,Shape))
    eOffsetM = eOffset*np.ones(Shape**2).reshape((Shape,Shape))
    Testspot = eCoeffM*(Testspot - eOffsetM)
    TestImage.close()
    return Testspot



#################################################################
# datafiles is the list of images in the stack ( should be in order)
# This function scans the data and finds maxima in a region that is the size of
# neighborhood_size. It looks for points in the neighborhood_size that are above
# the threshold. The threshold is defined as the mean of the data set + threshold_sigma
# deavations above the mean. the lower the threshold_sigma the more points you will get.
# also lowering the neighborhood_size will yeild more points
#
# Shape and EdgeCut are included so you can cut points out that are near the edge of the image
#################################################################
def Spot_finder(SummedData, neighborhood_size, threshold_sigma, Shape, EdgeCut):

    threshold = np.mean(SummedData)+threshold_sigma*np.std(SummedData)
    data_max = filters.maximum_filter(SummedData, neighborhood_size)
    maxima = (SummedData == data_max)
    data_min = filters.minimum_filter(SummedData, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    AllPairs = []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        y_center = (dy.start + dy.stop - 1)/2
        if x_center>EdgeCut and x_center<Shape-EdgeCut and y_center>EdgeCut and y_center<Shape-EdgeCut:
            AllPairs.append([x_center,y_center])
    return np.array(AllPairs)



#################################################################
# Takes the datafiles and the IDed points and computes the area for each point
# through the whole stack
#################################################################
def Spot_Area(Dict, KEYS, spots, Rbkg, Rsig):
    
    SpotInfo = []
    for keys in KEYS:
        for we in range(0,len(spots)):
            Xindex = int(spots[we][0])
            Yindex = int(spots[we][1])

            spotBKG = np.array(Dict[keys][Yindex-Rbkg:Yindex+Rbkg+1,Xindex-Rbkg:Xindex+Rbkg+1])
            spotSIG = np.array(Dict[keys][Yindex-Rsig:Yindex+Rsig+1,Xindex-Rsig:Xindex+Rsig+1])
            Signal     = spotSIG.sum()/(2*Rsig+1)**2
            Background = (spotBKG.sum()-spotSIG.sum())/((2*Rbkg+1)**2-(2*Rsig+1)**2)
            SpotInfo.append([Xindex,Yindex,Signal,Background])
    return SpotInfo




def moments(data):
    """Returns (height, x, y, width_x, width_y, offset)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y, np.mean(data)


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


def gaussian_2d(xy_mesh, amp, xc, yc, sigma_x, sigma_y, offset):
    # unpack 1D list into 2D x and y coords
    (x, y) = xy_mesh
    # make the 2D Gaussian matrix
    gauss = amp * np.exp(-((x - xc) ** 2 / (2 * sigma_x ** 2) + (y - yc) ** 2 / (2 * sigma_y ** 2))) / (
                2 * np.pi * sigma_x * sigma_y)
    gauss = gauss + offset
    # flatten the 2D Gaussian down to 1D
    return np.ravel(gauss)


def Fit_2D_Gaussian(array2D):
    #https://kippvs.com/2018/06/non-linear-fitting-with-python/
    Shape = array2D.shape[0]
    x = np.arange(Shape)
    y = np.arange(Shape)
    xy_mesh = np.meshgrid(x, y)

    guess_vals = moments(array2D)
    
    fit_params, cov_mat = optimize.curve_fit(gaussian_2d, xy_mesh, np.ravel(array2D), p0=guess_vals)

    fit_errors = np.sqrt(np.diag(cov_mat))
    expect = gaussian_2d(xy_mesh, *fit_params).reshape(np.outer(x, y).shape)
    fit_residual = (array2D - expect)/np.sqrt(expect)
    chi_squared = np.sum(fit_residual**2)/((Shape-1)*(Shape-1))
    #fit_Rsquared = 1 - np.var(fit_residual) / np.var(array2D)

    return fit_params, fit_errors, chi_squared



def Unique_Circle(r, Xoff, Yoff, Shape):
    # found the 1000000 is where it tops out
    theta = np.linspace(0,2*np.pi,1000000)
    Xcirc = np.rint(r*np.cos(theta)+Xoff)
    Ycirc = np.rint(r*np.sin(theta)+Yoff)
    # this finds the unique ones
    base=Xcirc.max()+1
    c=Xcirc+base*Ycirc
    val,ind=np.unique(c,return_index=True)
    
    X_circle = Xcirc[ind].astype(int)
    Y_circle = Ycirc[ind].astype(int)
    keep = np.where(Y_circle<Shape)
    Y_circle = Y_circle[keep]
    X_circle = X_circle[keep]
    keep = np.where(X_circle<Shape)
    Y_circle = Y_circle[keep]
    X_circle = X_circle[keep]
    
    return X_circle, Y_circle


class Trajectory():
    def __init__(self):
        name           = ""
        Mean_diff      = 0 
        Chi_Step       = 0
        Chi_Line       = 0
        Step_Time      = 0
        Times          = np.array([])
        Signal         = np.array([])
        Signal_er      = np.array([])
        Background     = np.array([])
        Background_er  = np.array([])
        
def Step_fit(Time, Signal, Sigma):
    ch = []
    for x in range(1,len(Time)-1):
        DataLeft  = Signal[:x]
        DataRight = Signal[x:]
        sigmasLeft  = Sigma[:x]
        sigmasRight = Sigma[x:]
        mean  = np.mean(Signal)
        meanI  = np.mean(DataLeft)
        meanE  = np.mean(DataRight)

        chiHIGH = np.sum((DataLeft  - meanI)**2/abs(sigmasLeft )**2)/len(DataLeft )
        chiLOW  = np.sum((DataRight - meanE)**2/abs(sigmasRight)**2)/len(DataRight)
        cc = chiHIGH+chiLOW
        ch.append(cc)
    chi = min(ch)
    loc = np.where(ch==chi)[0][0]+1 #+1 for the loop offset
    MeanDiff = (np.mean(Signal[:loc])) - ((np.mean(Signal[loc:])))
    return MeanDiff, chi, loc


def Line_fit(Time, Signal, Sigma):
    coefs=np.polyfit(Time,Signal,1)    
    ffit  = coefs[0]*Time+coefs[1]
    chiFit  = np.sum((ffit - Signal)**2/abs(Sigma)**2)/748
    return chiFit



def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n



## Fit 2d surfaces
def Surface_Fit(Data,Power):
    
    ypos, xpos  = np.indices(Data.shape) 
    #Znew = np.zeros(xpos.shape)
    X = xpos.flatten()
    Y = ypos.flatten()
    Z = Data.flatten()
    
    tmp_A = []
    tmp_b = []
    for i in range(len(X)):
        LIST = []
        for q in reversed(range(1,Power+1)):
            LIST.append(X[i]**q)
            LIST.append(Y[i]**q)
        LIST.append(1)
        tmp_A.append(LIST)
        tmp_b.append(Z[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    fit = (A.T * A).I * A.T * b
    return fit

def Fit_Surface(Data,Power,fit):
    ypos, xpos  = np.indices(Data.shape) 
    Znew = np.zeros(xpos.shape)

    for r in range(xpos.shape[0]):
        for c in range(xpos.shape[1]):
            LIST = []
            for q in reversed(range(1,Power+1)):
                LIST.append(xpos[r,c]**q)
                LIST.append(ypos[r,c]**q)
            LIST.append(1)

            Znew[r,c] = sum([x*y for x,y in zip(fit,LIST)])
    return Znew