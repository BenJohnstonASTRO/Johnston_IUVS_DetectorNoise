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



def find_all(pattern, path):
    """ FROM KYLE 
    Find all files with a specified pattern.

    Args:
        pattern: a Unix-style string to search for. Ex '*.pdf'
        path: a Unix-style string of the path to search for the name. Ex. '/Users/kyco2464/'

    Returns:
        a list of the complete paths containing the pattern
    """
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnm.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result



def findSTD(data):
    """
    Finds the standard deviation of a dataset which doesn't have to be gaussian
    
    n is the "standard deviation"
    
    """
    maxLoc = np.argmax(data) #Find the maximum location in the dataset
    
    #Now count from the peak to get where the data is 68.2% of the peak
    
    n = 0 #Initialize the counter
    tot = np.sum(data)
    
    while n < maxLoc:
        n = n+1
        area = np.sum(data[maxLoc-n:maxLoc+n]) #find the area under the curve from the peak down
        per = area/tot #calculate the 
        if per >=0.682:
            break
    return n



def getFiles(shape, a, orbitnums): 
    """
    This will give you a list of the orbits with a specifi data structure (shape, i.e. (2, 7, 256))
    shape -> tuple containing the shape of the data to make the cube in
    """

    #Initiaite the loop to load all images
    first = a[0] #first file name
    firstimg = fits.open(first) #open the first image
    firstimg_data = fits.getdata(first) #open the first image's data

    firstprim = firstimg[0] #index zero is where the data is stored
    firstprim_data = firstimg[0].data #.data gets the data
    cube = firstprim_data #renaming
    
    counter = 0 #count each time through the look, initialize
    countOrb = [] #empty lists to keep track of which orbit we are on later
    files = []
    orbitvals=[]
    for file in a[1:]: #loop through remaining files
        counter +=1 #step the counter
        image = fits.open(file) #open each image and data, same as before
        image_data = fits.getdata(file)

        prim = image[0]
        prim_data = image[0].data

        if prim_data.shape == shape: #if the current image opened is of the desired shape
            countOrb.append(counter) #append the count value
            files.append(a[counter]) #append the filename of this current image
            orbitvals.append(orbitnums[counter]) #and give the orbit number corresponding to this orbit
            #print(countOrb, '\n')
    
    return files, orbitvals #return the files and the orbit numbers



def getCube(files):
    First = files[0]
    Firstimg = fits.open(First)
    Firstimg_data = fits.getdata(First)

    #Ask Nick
    Firstprim = Firstimg[0]
    Firstprim_data = Firstimg[0].data
    Shape = Firstprim_data.shape
    print(Shape)

    Cube = Firstprim_data

    for File in files[1:]: #loop through remaining files 
        Image = fits.open(File)
        Image_data = fits.getdata(File)

        Prim = Image[0]
        Prim_data = Image[0].data
        
        Cube = np.dstack([Cube, Prim_data]) #Stack each image, dstack stacks in sequenc
    if np.array(Shape).shape == (3,):
        Cube = np.reshape(Cube, (Shape[0], Shape[1], len(files), Shape[2]))
    if np.array(Shape).shape == (2,):
        Cube = np.reshape(Cube, (Shape[0], len(files), Shape[1]))
    print(Cube.shape)
    return Cube




def findOutliers(CUBE, threshold, std):
    MoDe = mode(CUBE, axis=None)[0][0]
    maskedCube = ma.masked_less(CUBE, MoDe+(std*threshold))
    return maskedCube
    
    

    
def getOutlierPlot(cube, name, filenames, maskedCube, v_min, v_max):
    wri = ani.FFMpegWriter(fps=10)

    fig = plt.figure(figsize=(8,8), dpi=100)

    with wri.saving(fig, name, 100):

        data = np.log10(cube[:,:,:,:])

        for i in range(len(cube[0,0,:,0])):
            for j in range(len(cube[:,0,0,0])):

                plt.clf()

                b = filenames[i].split('_')
                plt.title("{}\nimg {}".format(b[6],j))
                plt.xlabel("Spectral Number")
                plt.ylabel("Spatial Number")

                width = len(cube[0,0,0,:])
                height = len(cube[0,:,0,0])

                plt.imshow(np.log10(maskedCube[j,:,i,:]), aspect=(width/height), origin="lower", cmap="CMRmap", vmin=v_min, vmax=v_max)

                cbar = plt.colorbar()
                cbar.set_label("Log Brightness (DN)")

                wri.grab_frame()




def getCubePlot(name, CUBE, a, v_min, v_max):
    wri = ani.FFMpegWriter(fps=10)

    #Define figure size
    fig = plt.figure(figsize=(8,8), dpi=150)

    #Tells where to save the file
    with wri.saving(fig, r'C:\Users\benja\Documents\Johnston_MAVEN\FurtherSepAnalysis\{}'.format(name), 100):

        data = np.log10(CUBE)

        for i in range(len(CUBE[0,0,:,0])): #Loop over all orbits
            for j in range(len(CUBE[:,0,0,0])): #Loop over the dark images taken per orbit

                plt.clf() #Clear all previous figures before executing  - "clean slate"

                b = a[i].split('_')
                plt.title("file: {}\nimg {}".format(b[6],j))
                plt.xlabel("Spectral Number")
                plt.ylabel("Spatial Number")

                width = len(CUBE[0,0,0,:]) 
                height = len(CUBE[0,:,0,0])

                plt.imshow(np.log10(CUBE[j,:,i,:]), aspect=(width/height), \
                           origin="lower", cmap="CMRmap", vmin=v_min, vmax=v_max)

                cbar = plt.colorbar() #Set a colorbar
                cbar.set_label("Log Brightness (DN)") #Set label

                wri.grab_frame()
                
                

def getCounts(CUBE, threshold):
    start = np.amin(CUBE)
    end = np.amax(CUBE)
    N_Arr = CUBE.flatten()
    MODE = mode(CUBE, axis=None)[0][0]

    plt.figure(figsize=(6, 4), dpi=150)
    hist_vals, Bins, Patches = plt.hist(N_Arr, histtype="step", range=(start, end), bins=int(end-start), label='all dark')
    
    STD = findSTD(hist_vals)
    
    plt.xscale("log")
    plt.title('Histogram of flattened arrays')
    plt.xlabel('Log DN')
    plt.ylabel('# Occurences')
    plt.axvspan(xmin=(MODE-STD),xmax=(MODE+STD), color="orange", alpha=.4, label="one sigma")
    plt.legend()

    MaskedCube = ma.masked_less(CUBE, MODE+(STD*threshold))
    #How many bins lie above a threshold, n*sigma
    
    #ThisCube = np.mean(ThisCube, axis=0)
    Counts = []
    for i in range(len(CUBE[0,0,:,0])): #Look over each orbit
        mask = np.mean(MaskedCube[:, :, i, :], axis=0)
        Num = ma.count(mask) #Count the non-masked elements, 
        Counts = np.append(Counts, Num)
    Counts = np.array(Counts)
    return Counts, STD




def outlierHist(inImg, outImg, bins = 100):
    inFlat = inImg.flatten()
    outFlat = outImg.flatten()
    #.flatten() stacks the values from each dimension

    plt.figure(figsize=(6, 4), dpi=80)
    plt.hist(inFlat, histtype="step", bins=bins, label='FUV Dark possible Event')
    plt.hist(outFlat, histtype="step", bins=bins, label='FUV Dark Typical')
    plt.xscale("log")
    plt.title('Histogram of flattened arrays')
    plt.xlabel('Log DN')
    plt.ylabel('# Occurences')
    plt.legend()
    
    
    
def doAThing(orbits, orbMax):
    return (np.where(np.array(orbits) <= orbMax))
    
    

    
def identifyEvent(orb, dark, orb2, dark2, data, outliers, fileList, aspect):
    #ANALYZING ORBIT 2182

    plt.figure(figsize=(10,4), dpi=300)
    data = np.log10(data)
    MIN = np.amin(data)
    MAX = np.amax(data)

    typical = (data[dark, :, orb, :])
    plt.subplot(1,2,1)
    plt.imshow(typical, aspect=(256/7.0), origin="lower", cmap="CMRmap", vmin=MIN, vmax=MAX)
    plt.title('orbit {}'.format(fileList[orb].split('_')[6], dark))
    cbar = plt.colorbar()
    plt.xlabel('Spatial Bin')
    plt.ylabel('Spectral Bin')

    possible = (data[dark2, :, orb2, :])
    plt.subplot(1,2,2) #LOG IMAGE
    plt.imshow(possible, aspect=(256/7.0), origin="lower", cmap="CMRmap", vmin=MIN, vmax=MAX)
    cbar = plt.colorbar()
    plt.title('orbit {}, dark {}'.format(fileList[orb2].split('_')[6], dark2))
    cbar.set_label("Log Brightness (DN)")
    plt.xlabel('Spatial Bin')
    plt.ylabel('Spectral Bin')
    plt.tight_layout()


    plt.figure(figsize=(10,4), dpi=300)
    dataOutlier = np.log10(outliers)

    typicalOutlier = (dataOutlier[dark, :, orb, :])
    plt.subplot(1,2,1)
    plt.imshow(typicalOutlier, aspect=(aspect), origin="lower", cmap="CMRmap", vmin=MIN, vmax=MAX)
    plt.title('Outliers, orbit {}'.format(fileList[orb].split('_')[6], dark))
    plt.text(220, 6.1, "{} bins".format(ma.count(typicalOutlier)), ha="center", va="center", bbox=dict(boxstyle="round", fc=(1,1,1)))
    cbar = plt.colorbar()
    plt.xlabel('Spatial Bin')
    plt.ylabel('Spectral Bin')

    possibleOutlier = (dataOutlier[dark2, :, orb2, :])
    plt.subplot(1,2,2) #LOG IMAGE
    plt.imshow(possibleOutlier, aspect=(aspect), origin="lower", cmap="CMRmap", vmin=MIN, vmax=MAX)
    plt.text(220, 6.1, "{} bins".format(ma.count(possibleOutlier)), ha="center", va="center", bbox=dict(boxstyle="round", fc=(1,1,1)))
    cbar = plt.colorbar()
    plt.title('Outliers, orbit {}, dark {}'.format(fileList[orb2].split('_')[6], dark2))
    cbar.set_label("Log Brightness (DN)")
    plt.xlabel('Spatial Bin')
    plt.ylabel('Spectral Bin')
    plt.tight_layout()


    outlierHist(possible, typical, 100)