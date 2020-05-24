# -----------------------------------------------------------------------------
#  REST_SMFI | Trajectory_Maker.py
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

