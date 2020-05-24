# -----------------------------------------------------------------------------
#  REST_SMFI | Spot_Identification.py
#
#  Various functions for data analysis
#   * Author: Austin McDonald
#   * Creation date: Oct 2019
# -----------------------------------------------------------------------------

from PIL import Image
import numpy as np

# Nedded for background fit
import numpy.polynomial.polynomial as poly
# needed for spotfinder
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

from scipy import optimize



#################################################################
# Make a circle for cutting the image spot out.
#################################################################
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



