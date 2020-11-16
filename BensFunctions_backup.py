import math
import numpy as np
import numpy.ma as ma
import matplotlib
from astropy.io import fits
from scipy.stats import mode
import matplotlib.pyplot as plt
import os, sys
from scipy.optimize import curve_fit
import matplotlib.animation as ani
get_ipython().magic(u'matplotlib inline')
import fnmatch as fnm
import scipy.special as sse
from scipy.stats import chisquare
import Johnston_MakeCube as mc
import BensFunctions as bf
import heapq
import itertools


def getCube(files):
    """
    This function makes a data cube based on a given list of fits files. data cube
    
    --------------------------------------------------------------------
    Inputs: files -> list of filenames, including full directory
    
    Example file name: 'C:\\Users\\benja\\Documents\\Johnston_MAVEN\\Sept_SEP_Event\\mvn_iuv_l1a_incorona-orbit05700- \
                fuvdark_20170907T130400_v13_r01.fits.gz'
    """
    
    First = files[0] #Open the first file in the inputed list
    Firstimg = fits.open(First) #Open the image
    
    Firstprim = Firstimg[0] #this gives the header
    Firstprim_data = Firstimg[0].data #Get data from the header
    Shape = Firstprim_data.shape #Find the shape of the first file
    print(Shape)

    Cube = Firstprim_data #renaming

    for File in files[1:]: #loop through remaining files 
        Image = fits.open(File)
        Image_data = fits.getdata(File)

        Prim = Image[0]
        Prim_data = Image[0].data
        
        Cube = np.dstack([Cube, Prim_data]) #Stack each image, dstack stacks in sequence
        
    if np.array(Shape).shape == (3,):
        Cube = np.reshape(Cube, (Shape[0], Shape[1], len(files), Shape[2]))
    if np.array(Shape).shape == (2,):
        Cube = np.reshape(Cube, (Shape[0], len(files), Shape[1]))
    print(Cube.shape)
    return Cube


def STD(data):
    """
    Finds the standard deviation of a dataset which doesn't have to be gaussian,  n is the "standard deviation"
    
    inputs: data -> histogram values, i.e. data = vals when vals, bins, patches = plt.hist(image.flatten(), ...)
    
    """
    maxLoc = np.argmax(data) #Find the maximum location in the dataset
    
    #Now count from the peak to get where the data is 68.2% of the peak
    
    n = 0 #Initialize the counter
    tot = np.sum(data)
    
    while n < maxLoc: #while the std is less than 0.682
        n = n+1 #update
        area = np.sum(data[maxLoc-n:maxLoc+n]) #calculate the area from the top down until 68% is filled in
        per = area/tot  #calculate percent
        if per >=0.682: #if the percentage is over 68%, break and n is the std, as it is the number of steps out from the peak
            break
    return n