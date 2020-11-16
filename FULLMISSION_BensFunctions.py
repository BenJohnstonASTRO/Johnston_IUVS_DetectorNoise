#Import necessary modules

import numpy as np
import numpy.ma as ma
import matplotlib
from astropy.io import fits
from scipy.stats import mode
import matplotlib.pyplot as plt
import os, sys
import math
import fnmatch as fnm
plt.rcParams.update({'figure.max_open_warning': 0})



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



def getFiles(Type, directory):
    """
    
    This will make a list of .fits image files of a specific observing type in a given directory
    
    ---Inputs---
    1. Type - the type of observation, e.g. apoapse, periapse, inlimb, outlimb, outdisk, indisk, etc
            - format: string, e.g. 'apoapse'
    2. directory - the folder in which to look for the .fits files
                 - formart r'{directory name here}'
        
    ---Outputs---
    1. s - a list of files with the given name in a directory 
    
    ---Example---
    s = getFiles('apoapse', directory)
    print('Number of files: {}'.format(len(s)))
    
    """
    
    plt.clf() #Clear any unwanted plots for memory's sake
    
    s = find_all("*{}*.fits*".format(Type), directory) #Use the function from Kyle to find the files
    
    s.sort() #re-arrange the files in order of ascending orbit #, if necessary
    
    return s #return the list of files



def getOrbitRange(files):
    """
    
    Returns the orbit range for the list of files provided
    
    ---Inputs---
    1. files - list of files, i.e. the output from the getFiles() function
    
    ---Outputs---
    1. allOrbits - a list of the orbits associated with the files
    
    ---Example---
    orbits = getOrbitRange(s) #Get the orbits
    print('Number of orbits: {}'.format(len(orbits))) #get the number of orbits/files
    
    """
    
    allOrbits = [] #Define an empty list to add orbit number to in the next steps
    
    for i in range(len(files)):  #Loop over all the files, len(files) gives the number of files
        orbit = files[i][:-40][-5:]  #this indexes the orbit value, as a string, out of each filename
        orbit = int(orbit) #Convert the string to an integer   
        allOrbits.append(orbit) #add the orbit to the list
    
    allOrbits = np.array(allOrbits) #Convert the list of orbits to an array, to make calculations easier
    
    return allOrbits #return the orbit array



def getSpecifics(files, setShape, allOrbits, orbitRange, Type):
    """
    
    This function gets the shape of each file, the number of spectral bins, and the associated orbits.
    There is an optional argument, 'orbitRange', which will "zoom in" on a specified orbit range
    The function will return the orbits, shapes, and files in a specified orbit range of a specified type and shape
    
    ---Inputs---
    1. files - a list of files, i.e. the output from 'getFiles'
    2. setShape - a tuple, i.e. (2, 7, 256), to represent the desired shape
    3. allOrbits - complete orbit list, i.e. from the outputof 'getOrbitRange'
    4. orbitRange - enter a tuple, i.e. (5700, 5790), to select an orbit range
    5. Type - the obersvational mode, i.e. apoapse, periapse, as a string
    
    ---Outputs---
    1. Bins - a list of the spectral bins corresponding to the desired shape
    2. Shapes - a list of the shapes corresponding to the desired shape
    3. Orbits - all orbits corresponding to the desired shape
    4. Files - a complete list of files corresponding to the desired shape
    5. Headers - a list of every header from the files
    
    ---Example---
    bins, files, shapes, orbits = getSpecifics(s, (2, 10, 256), orbits, (0, 12000))
    
    """

    shapes = [] #create empty lists for the shape, # bins, orbits, and headers
    Bins = []
    sNew = []
    Orbits = []
    headers = []


    for i in range(len(files)): #Loop over all files
        file = files[i] #load each file
        image = fits.open(file) #open each image and data
        image_data = fits.getdata(file) #get the data
        prim = image[0] #Save the data to prim_data
        prim_data = image[0].data
        header = prim.header #Find the header
        binVal = header[3] #index the spectral bin value out of the header

        shape = prim_data.shape #Find the shape of the file
        
        if (allOrbits[i] >= orbitRange[0]) and (allOrbits[i] <= orbitRange[1]): 
            #If the orbit is within the specified range
            if shape == setShape: #If the shape of the file matches the desired shape, do the following
                headers.append(header) #add the header, shape, bin, and orbit Number to the empty list
                sNew.append(file)
                shapes.append(shape)
                Bins.append(binVal)
                Orbits.append(allOrbits[i])

    Bins = np.array(Bins) #Rename the original variables
    Files = np.array(sNew)
    Shapes = np.array(shapes)
    Orbits = np.array(Orbits)
    
    if len(Orbits) == 0:
        print('No orbits of shape {} in orbit range {}-{} of type {}'.format(setShape, orbitRange[0], orbitRange[1], Type))
    
    return (Bins, Files, Shapes, Orbits, headers) #return the bins, files, shapes, and orbits



def getSTD(data, width=1):
    """
    
    Finds the standard deviation of a dataset which doesn't have to be gaussian,  n is the "standard deviation"
    
    ---Inputs---
    1. data - histogram values, i.e. data = vals when vals, bins, patches = plt.hist(cube.flatten(), ...)
       the histogram should be the histogram of the entire data cube
    
    ---Outputs---
    1. n*width - the pseudo standard deviation
    
    ---Example---
    v, b, p = plt.hist(cube.flatten(), histtype="step", range=(vmi, vma), bins=int(vma-vmi))

    std = STD(v)
    
    """
    
    maxLoc = np.argmax(data) #Find the maximum location in the dataset
    
    #Now count from the peak to get where the data is 68.2% of the peak
    
    n = 0 #Initialize the counter
    tot = np.sum(data) #find the total value of the data combined
    
    while n < maxLoc: #while the std is less than 0.682
        n = n+1 #update
        area = np.sum(data[maxLoc-n:maxLoc+n]) #calculate the area from the top down until 68% is filled in
        per = area/tot  #calculate percent
        if per >=0.682: #if the percentage is over 68.2%, stop the code;  n is the std, as it is the number of steps out from the peak
            break
    return n*width



def makeCube(files, setShape):
    """
    Makes a data cube, i.e. a "stack" of spatial-spectral images. The dimensions a re as follows:
    -if setshape has 3 elements, i.e. (2, 7, 256), and the number of files is 100, the shape of the cube will 
     be (2, 7, 100, 256) - 2 is the #darks per file, 7 is the # spatial bins, 100 is the #files, 256 is # spectral bins
    -if setShape has 2 elements, i.e. (7, 165), there is only 1 dark, and the cube will be (7, 100, 165)
    
    ---Inputs---
    1. files - list of files, from getFiles output
    2. setShape - desired shape of files for cube
    
    ---Outputs---
    1. cube - a data cube containing spatial-spectral-time images
    
    ---Example---
    cube = makeCube(s, (2, 10, 256))
    
    """

    first = files[0]
    firstimg = fits.open(first)
    firstimg_data = fits.getdata(first)

    firstprim = firstimg[0]
    firstprim_data = firstimg[0].data

    FourDCube = firstprim_data

    for file in files[1:]: #loop through remaining files 
        image = fits.open(file)
        image_data = fits.getdata(file)

        prim = image[0]
        prim_data = image[0].data
        shape = prim_data.shape

        FourDCube = np.dstack([FourDCube, prim_data]) #Stack each image, dstack stacks in sequence

    if len(setShape) == 3:
        cube = np.reshape(FourDCube, (setShape[0], setShape[1], len(files), setShape[2]))
    if len(setShape) == 2:
        cube = np.reshape(FourDCube, (setShape[1], len(files), setShape[2]))
    print('\n \nThe shape of the cube is {}'.format(cube.shape))
    return cube



def getHist(cube, scale, threshold=6, limit=False):
    """
    Makes a histogram of a data cube
    
    ---Inputs---
    1. cube - the data cube for a specific shape, i.e. the output from makeCube()
    2. scale - log or lin, sets x-scale
    3. threshold - threshold to look for outliers, default 6
    4. limit - enter a 2x1 tuple(i.e. (0, 10000)) to limit the x-span of the histogram plot
    
    ---Outputs---
    1. vmi - minimum value of all the images
    2. vma - maximum value of all the histograms
    3. cubeSTD - the standard deviation of the histogram, calcualted using getSTD
    4. histogram plot - a plot of the histogram
    
    ---Example--
    minVal, maxVal, std, mode = getHist(cube, 'log')
    
    """
    
    vmi = np.amin((cube)) #Calculate the minimum and maximum of the cube
    vma = np.amax((cube))
    if math.isnan(float(vmi)) == True:
        vmi = 1000
        print('minimum value is nan, using default minimum of 1000')
    if math.isnan(float(vma)) == True:
        vma = 30000
        print('maximum value is nan, using default maximum of 30000')
    
    
    print('min = {},  max= {}'.format(vmi, vma))
    
    #Calculate the mode of the cube, i.e. the value that appears most often in the cube
    cubeMode = mode(cube, axis=None)[0][0] 

    plt.figure(figsize=(8, 5), dpi=100) #make a new matplotlib figure
    vals, bins, patches = plt.hist(cube.flatten(), histtype="step", range=(vmi,vma), bins=int(vma-vmi))
    
    cubeSTD = getSTD(vals) #Calculate the standard deviation of the histogram
    
    plt.axvline(cubeMode, label='mode', color='r') #plot a line at the mode and the threshold on the histogram and the 
    plt.axvline(cubeMode + cubeSTD*threshold, color='k', label='{}\u03C3 threshold'.format(threshold))
    #Calculate the histogram
    plt.xlabel('Log # DN (data number)') #set x, y, and title labels
    plt.ylabel('frequency')
    plt.title('cube histogram')
    plt.xscale(scale) #set scale, lin or log
    plt.legend() #show the legend
    if limit == True: #plot limits  if specified by user
        plt.xlim(limit[0], limit[1])
    plt.show() #show plot
    
    return (vmi, vma, cubeSTD, cubeMode)



def plotPostage(cube, orbits, setShape, dark):
    """
    This makes "postage stamp" size plots of the fuvdark images using CMRmap
    
    ---Inputs---
    1. cube - the data cube, i.e. the output from makeCube()
    2. orbits - the complete list of lrbits, i.e. the orbit output from getSpecifics()
    3. setShape - the shape of the cube
    4. dark - the desired dark frame to use, either 0 (for the first) or 1 (for the second)
    
    ---Outputs---
    1. matplotlib figures - postage stamp images
    
    ---Example---
    plotPostage(cube, allOrbits, (2, 7, 256), 0)
    
    """
    
    vmi = np.amin((cube)) #Calculate the minimum and maximum of the cube
    vma = np.amax((cube))
    if math.isnan(float(vmi)) == True:
        vmi = 1000
        print('minimum value is nan, using default minimum of 1000')
    if math.isnan(float(vma)) == True:
        vma = 30000
        print('maximum value is nan, using default maximum of 30000')
    
    stop = 1 #This is a counter for determining the number of stamps to make in a row

    plt.figure(figsize=(20, 8), dpi=100) #Make a figure
    plt.tight_layout(w_pad=6) #Make each subplot compacts
        
    for i in range(len(orbits)): #loop over all orbits
        
        if len(setShape) == 3: #Determione how to index cube based on the shape of the inputted cube
            data = cube[dark, :, i, :]
            
        if len(setShape) == 2: #When the shappe is 2, there is implicitly only 1 dark
            data = cube[:, i, :]
    
        plt.subplot(1, 8, stop) #make a new subplot

        plt.title("orbit {}".format(orbits[i])) #Set title
        
        #Show the image, vmin and vmax set the colormap bounds to span from lowest values to possible aurora
        plt.imshow(np.log10(data), aspect=setShape[-1]/setShape[-2], cmap="CMRmap", vmin=np.log10(vmi), vmax=np.log10(vma))
        
        stop +=1 #Update the counter

        if stop == 9: #When the counter reaches 9, make a new row for the subplots, i.e. only 8 subplots per row
            plt.figure(figsize=(20, 8)) #Make a new figure
            plt.tight_layout(w_pad=6)
            stop = 1 #reset the counter
            
            
            
def plotPostageHist(cube, orbits, setShape, dark, cubeMode, std, binwidth=100, threshold=6):
    """
    This makes "postage stamp" size phistogram plots of the fuvdark images comparing each to a typical FUVDarl histogram
    
    ---Inputs---
    1. cube - the data cube, i.e. the output from makeCube()
    2. orbits - the complete list of lrbits, i.e. the orbit output from getSpecifics()
    3. setShape - the shape of the cube
    4. dark - the desired dark frame to use, either 0 (for the first) or 1 (for the second)
    5. cubeMode - the mode  of the cube, i.e. the "mode" output from getHist()
    6. std - the std of the cube, ie. the "std" output from getHist()
    7. binwidth - the width of the histogram bins, default to 100 based on experimenting.
    8. threshold - the threshold to look for outliers, default is 6
    
    ---Outputs---
    1. matplotlib figures - postage stamp images
    
    ---Example---
    plotPostageHist(cube, allOrbits, (2, 7, 256), 0, Mode, std)
    
    """
    
    vmi = np.amin(cube) #Find the min and max of the cube
    vma = np.amax(cube)
    if math.isnan(float(vmi)) == True:
        vmi = 600
        print('minimum value is nan, using default minimum of 600')
    if math.isnan(float(vma)) == True:
        vma = 35000
        print('maximum value is nan, using default maximum of 35000')
    
    stop = 1 #This is a counter for determining the number of stamps to make in a row

    plt.figure(figsize=(20, 2), dpi=100) #Make a figure
    plt.tight_layout(w_pad=6) #Make each subplot compacts
        
    for i in range(len(orbits)): #loop over all orbits
        
        if len(setShape) == 3: #Determione how to index cube based on the shape of the inputted cube
            flatData = cube[dark, :, i, :].flatten() #Flatten the data to make a histogram
            typicalData = cube[dark, :, -1, :].flatten() #Make a typical file too
            
        if len(setShape) == 2: #When the shappe is 2, there is implicitly only 1 dark
            flatData = cube[:, i, :].flatten()
            typicalDatadata = cube[:, i, :].flatten()
    

        plt.subplot(1, 7, stop)

        b = orbits[i] #get the orbit
        plt.title("orbit {}".format(b)) #set title and labels
        plt.xlabel("log DN")
        plt.ylabel("frequency")
        plt.xscale("log")
        plt.xlim(left=vmi, right=vma)  
        
        #make the histogram of each image
        vals, bins, patches = plt.hist(flatData, histtype="step", bins=range(int(vmi),int(vma+binwidth),int(binwidth)), label='E')

        #make a "typical" histogram for comparison
        vals, bins, patches = plt.hist(typicalData, histtype="step", bins=range(int(vmi),int(vma+binwidth),int(binwidth)), label='T')

        plt.axvline(x=(cubeMode+(std*threshold)), color="k", label="""{}\u03C3""".format(threshold)) #make a line at the 6 sigma threshold
        
        if i == 0: #only plot the legend on the very first image
            plt.legend()
        stop +=1

        if stop == 8:
            plt.figure(figsize=(20, 2), dpi=100)
            plt.tight_layout()
            stop = 1
            
            
            
def findOutliers(orbits, cube, cubeMode, std, dark, setShape, threshold=6):
    """
    This function finds the number of "outliers" in an image, i.e. the number of pixels that are above a threshold
    I.e. the number of bins that are above 6 sigma of a typical FUVdark image, by making a "mask"
    A mask deletes values less than a specific values
    
    ---Inputs---
    1. orbits - the list of orbits
    2. cube - the data cube
    3. mode - the mode of the cube
    4. std - the standard deviation of the cube
    5. dark - the dark to use, either 0 or 1
    6. setShape - desired shape, must match other functions
    7. threshold - default to 6 for 6-sigma anlaysis, can change to other value 
    
    ---Outputs---
    1. counts - list of the number of outliers, linear
    2. percents - list of the percent outliers
    3. medians - list of the medians of each image
    
    ---Example---
    counts, percents, medians = findOurliers(orbits, cube, cubeMode, cubeStd, 0)

    """
    
    #mask all values less than the mode+std*threshold, i.e. mask values below 6 sigma
    mask = ma.masked_less(cube, cubeMode+(std*threshold)) #make the mask

    counts = [] #make empty lists for the counts, median values, and percent values
    medians = []
    percents = []

    for i in range(len(orbits)):#Loop over all the orbits
        if len(setShape) == 3: #Determione how to index cube based on the shape of the inputted cube
            data = cube[dark, :, i, :]
            maskedData = mask[dark, :, i, :]
            
        if len(setShape) == 2: #When the shappe is 2, there is implicitly only 1 dark
            data = cube[:, i, :]
            maskedData = mask[:, i, :]
            
        count = ma.count(maskedData) #Count the number of values in the mask
        med = np.median(data) #find the median of the deata

        counts.append(count) #append the number of counts
        medians.append(med) #append the median value of the data
        perCount = 100*count/ma.count(data) #find the percent outliers by dividing by the total pixels
        percents.append(perCount) #append the percent value
    
    counts = np.array(counts) #make the lists arrays for easier computation later
    percents = np.array(percents)
    medians = np.array(medians)
    
    return (counts, percents, medians)



def plotOutliers(orbits, counts, percents, medians, dark, THRESHOLD=6):
    """
    This makes4 plots: orbit vs counts, orbit vs log10(counts), and orbit vs percent counts above 6 sigma, and orbit vs median of each image
    
    ---Inputs---
    1. orbits - list of orbits
    2. counts - list of counts, i.e. the output from findOutliers()
    3. percents - list of percents, i.e. the output from findOutliers()
    4. THRESHOLD - default to 6, threshold to look for outliers. MUST MATCH THE INPUT FOR findOutliers()
    
    ---Outputs---
    1. plots of orbits vs counts, log10counts, and percent counts
    
    ---Example---
    plotOutliers(orbits,  counts, percents)
    
    
    """
    
    plt.figure(figsize=(10,4), dpi=100)
    plt.tight_layout()
    
    #PLOT LINEAR OUTLIERS
    plt.subplot(1, 2, 1)
    plt.scatter(orbits, counts, s=5, c='k')
    plt.ylabel('# outliers above {}\u03C3'.format(THRESHOLD))
    plt.xlabel('Orbit number')
    plt.title('{}\u03C3 Outliers vs Orbit Number ({} dark only)'.format(THRESHOLD, dark))

    #PLOT LOG OUTLIERS
    plt.subplot(1, 2, 2)
    plt.scatter(orbits, np.log10(counts), s=5, c='r')
    plt.ylabel('Log # outliers above {}\u03C3'.format(THRESHOLD))
    plt.xlabel('Orbit number')
    plt.title('Log {}\u03C3 Outliers vs Orbit Number ({} dark only)'.format(THRESHOLD, dark))
    plt.show()
    
    #PLOT PERCENT OUTLIERS
    plt.figure(figsize=(10,4), dpi=100)
    plt.tight_layout()

    plt.subplot(1, 2, 1)
    plt.scatter(orbits, percents, s=5, c='b')
    plt.ylabel('# outliers above {}\u03C3 / total pixel count'.format(THRESHOLD))
    plt.xlabel('Orbit number')
    plt.title('Percent {}\u03C3 Outliers vs Orbit Number ({} dark only)'.format(THRESHOLD, dark))
    plt.tight_layout()
    
    #PLOT MEDIANS
    plt.subplot(1, 2, 2)
    plt.scatter(orbits, medians, s=5, c='b')
    plt.ylabel('Median value of each image'.format(THRESHOLD))
    plt.xlabel('Orbit number')
    plt.title('Median vs Orbit Number ({} dark only)'.format(dark))
    plt.tight_layout()
    plt.show()

    
    
def plotPrettyLogCounts(orbits, counts, dates, c, m, s, saveName, legendLoc, ylim, Type, step, dark, frameOff=False, thres=6):
    """
    This makes a custom plot of orbits vs log(counts) for comparison to other figure (i.e. f rom papers, etc)
    The figure has 2 x-axes: the lower is orbit # and the upper is the corresponding date, the y axis is on the right side
    
    ---Inputs---
    1. orbits - the list of orbits, i.e. the output of getSpecifics()
    2. counts - the counts, i.e. the output of findOutliers()
    3. dates - a list of dates corresponding to the orbit list, get the dates from the headers and manually enter a list
    4. c - the coolor of the plot, i.e. 'red', see matplotlib documentation, string
    5. m - the marker style, see matplotlib documentation "markersStyle", string
    6. s - the marker size for the scatter plot, number
    7. saveName - string ending in .png, i.e. 'test.png', the name of the file to save to
    8. legendLoc - 2-number tuple, determines where to place the legend on the plot
    9. ylim - this specifies the ylimits, enter a 2-number tuple
    10. Type - the observational mode, must match other functions, enter a string
    11. step - how many orbits to step between x-ticks, enter a number
    12. dark - number, which dark to use if more than 1 are available
    13. frameOff - this is default to False, where the axes will be shown.  If set to True, the axes will be hidden
    14. thres - threshold to look for outliers, default is 6. must match other functions
    
    ---Outputs---
    1. A plot of orbits vs log counts with dual axes, orbit # and corresponding date
    
    ---Example---
    dates = ['Sep 07', 'sep 09', 'Sep 11', 'sep 13', 'Sep 14', 'sep 16', 'Sep 18', 'sep 20', 'Sep 22', 'sep 24', 'Sep 26']
    plotPrettyLogCounts(orbits, cubeCounts, dates, 'r', '*', 5, 'test.png', (0.85, 0.94), (0, 2), 'apoapse', 10)
    
    """  
    
    f = plt.figure(figsize=(10,4), dpi=300) #make the figure and axes
    ax = f.add_subplot(111) #make the bottom x-axis
    ax2 = ax.twiny() #make the top x-axis

    ax.yaxis.tick_right() #set the y-axis positions
    ax.yaxis.set_label_position("right")
    ax.yaxis.label.set_color(c) #set the y-axis color and position
    ax.set_ylabel('Log outliers above {}\u03C3'.format(thres))

    if frameOff == False: #only do the following if the frame axes are shown
        plt.ylabel('Log Outliers above {}\u03C3'.format(thres)) #set the x and y labels, titles, and ylimits
        plt.xlabel('Orbit date')
        plt.title('Log {}\u03C3 Outliers vs Orbit Number ({} dark only)'.format(thres, dark))
        ax.set_xlabel('Orbit #') #set label
    plt.ylim(ylim[0], ylim[1])
    
    plt.scatter(orbits, np.log10(counts), label=Type, marker=m, c=c, s=s) #plot the orbits vs log10(counts)

    ax2.set_xticks(np.arange(min(orbits), max(orbits), step=step)) #set the top axis with the entered string of dates
    ax2.set_xticklabels(dates, fontsize=8) #set labels
    ax2.set_xlim(left=min(orbits), right=max(orbits)) #set limits

    ax.set_xticks(np.arange(min(orbits), max(orbits), step=step)) #set the bottom axis with orbit numbers
    ax.set_xlim(left=min(orbits), right=max(orbits)) #set labels
    
    if frameOff == True: #this determines whether or not to turn the axes off
        #if true, hide all axes labels and tick marks
        ax.set_xticks([]) #turn off xtick marks and labels
        ax.set_xticklabels([])
        ax.set_xlabel('')
        ax2.set_xticks([])
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax.set_yticks([])
        ax.set_ylabel('')
        ax.set_yticklabels('')

    plt.legend(loc=legendLoc, frameon=False); #show legend
    plt.savefig(str(saveName), bbox_inches='tight', transparent=True) #save the figure and show the plot
    
    
    
    
def plotPrettyPercents(orbits, percents, dates, c, m, s, saveName, legendLoc, ylim, Type, step, dark, frameOff=False, thres=6):
    """
    This makes a custom plot of orbits vs percent counts (#Outliers/total pixel count) for comparison to other 
    figure (i.e. from papers, etc)
    The figure has 2 x-axes: the lower is orbit # and the upper is the corresponding date, the y axis is on the right side
    
    ---Inputs---
    1. orbits - the list of orbits, i.e. the output of getSpecifics()
    2. percentss - the percent counts, i.e. the output of findOutliers()
    3. dates - a list of dates corresponding to the orbit list, get the dates from the headers and manually enter a list
    4. c - the color of the plot, i.e. 'red', see matplotlib documentation, string
    5. m - the marker style, see matplotlib documentation "markersStyle", string
    6. s - the marker size for the scatter plot, number
    7. saveName - string ending in .png, i.e. 'test.png', the name of the file to save to
    8. legendLoc - 2-number tuple, determines where to place the legend on the plot
    9. ylim - this specifies the ylimits, enter a 2-number tuple
    10. Type - the observational mode, must match other functions, enter a string
    11. step - how many orbits to step between x-ticks, enter a number
    12. dark - number, which dark to use if more than 1 are available
    13. frameOff - this is default to False, where the axes will be shown.  If set to True, the axes will be hidden
    14. thres - threshold to look for outliers, default is 6. must match other functions
    
    ---Outputs---
    1. A plot of orbits vs log counts with dual axes, orbit # and corresponding date
    
    ---Example---
    dates = ['Sep 07', 'sep 09', 'Sep 11', 'sep 13', 'Sep 14', 'sep 16', 'Sep 18', 'sep 20', 'Sep 22', 'sep 24', 'Sep 26']
    plotPrettyLogCounts(orbits, cubeCounts, dates, 'r', '*', 5, 'test.png', (0.85, 0.94), (0, 2), 'apoapse', 10)
    
    """  
    
    f = plt.figure(figsize=(10,4), dpi=300) #make the figure and axes
    ax = f.add_subplot(111) #make the bottom x-axis
    ax2 = ax.twiny() #make the top x-axis

    ax.yaxis.tick_right() #set the y-axis positions
    ax.yaxis.set_label_position("right")
    ax.yaxis.label.set_color(c) #set the y-axis color and position
    ax.set_ylabel('Outliers above {}\u03C3 / total pixel count'.format(thres))

    if frameOff == False: #only do the following if the frame axes are shown
        plt.ylabel('Outliers above {}\u03C3 / total pixel count'.format(thres)) #show labels and title if the axes are shown
        plt.xlabel('Orbit date')
        plt.title('Percent {}\u03C3 Outliers vs Orbit Number ({} dark only)'.format(thres, dark))
        ax.set_xlabel('Orbit #') #set label
    plt.ylim(ylim[0], ylim[1])
    
    plt.scatter(orbits, percents, label=Type, marker=m, c=c, s=s) #scatter plot of orbits vs percents

    ax2.set_xticks(np.arange(min(orbits), max(orbits), step=step)) #set the top axis with the entered string of dates
    ax2.set_xticklabels(dates, fontsize=8) #set labels
    ax2.set_xlim(left=min(orbits), right=max(orbits)) #set limits

    ax.set_xticks(np.arange(min(orbits), max(orbits), step=step)) #set the bottom axis with orbit numbers
    ax.set_xlim(left=min(orbits), right=max(orbits)) #set labels
    
    if frameOff == True: #this determines whether or not to turn the axes off
        #if true, hide all axes labels and tick marks
        ax.set_xticks([]) #turn off xtick marks and labels
        ax.set_xticklabels([])
        ax2.set_xticks([])
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax.set_yticks([])
        ax.set_ylabel('')
        ax.set_yticklabels('')
        ax.set_xlabel('')

    plt.legend(loc=legendLoc, frameon=False); #show legend
    plt.savefig(str(saveName), bbox_inches='tight', transparent=True) #save the figure and show the plot
    
    
    

def plotTemps(orbits, files):
    """
    """

    caseT = []
    detT = []
    for i in range(len(orbits)):
        file = files[i] #load each file
        image = fits.open(file)

        case = image[1].data[0][-2]
        det = image[1].data[0][-1]
        caseT.append(case)
        detT.append(det)

    caseT = np.array(caseT)
    detT = np.array(detT)

    fig = plt.figure(figsize=(8,5), dpi=100) #make the figure and axes
    ax1 = fig.add_subplot(111)

    ax1.set_xlabel('Orbit #')
    ax1.set_ylabel('Case Temp (DN)', color='k')
    ax1.scatter(orbits, caseT, c='k', label='Case T', s=8, alpha=0.8)
    ax1.tick_params(axis='y', labelcolor='k')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Detector Temp (DN)', color='r')  # we already handled the x-label with ax1
    ax2.scatter(orbits, detT, c='r', label='Det T', s=8, alpha=0.5)
    ax2.tick_params(axis='y', labelcolor='r')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.title('Case & Detector Temperature orbits {}-{}'.format(min(orbits), max(orbits)))
    plt.show();
    
    
    
    
def runAll(directory, Type, dark, shape, orbitRange, scale, dates, name0, name1, loc, lim0, lim1, c, m, s, step, F=False):
    """
    This executes all functions, in sequence
    """

    if F == False:
        files = getFiles(Type, directory) #Get a list of all files
    else:
        files = F
    orbits = getOrbitRange(files) #Get a list of all orbits

    #print('TOTAL Number of orbits of type {}: {}'.format(Type, len(orbits)))
    #print('TOTAL Number of files of type {}: {}'.format(Type, len(files)))

     #Get files, shapes, bins, and orbits of a specific shape and update files and orbits
    bins, files, shapes, orbits, headers = getSpecifics(files, shape, orbits, orbitRange,  Type)

    #print('\n \nNumber of orbits with shape {} and type {}: {}'.format(shape, Type, len(orbits)))
    print('Number of files with shape {} and type {}: {}'.format(shape, Type, len(files)))
    
    if len(orbits) != 0:
        Cube = makeCube(files, shape) #Make a cube
        Min, Max, cubeStd, cubeMode = getHist(Cube, scale) #get the min and max of the cube, std, and mode

        cubeCounts, cubePercents, cubeMedians = findOutliers(orbits, Cube, cubeMode, cubeStd, dark, shape) #get the counts
        
        plotPostage(Cube, orbits, shape, dark) #Plot the images

        #plotOutliers(orbits, cubeCounts, cubePercents, cubeMedians, dark) #plot the orbits vs outliers
        
        #plotTemps(orbits, files) #make orbit vs case temp and orbit vs detector temp plots

        return (Cube, cubeCounts, cubePercents, cubeMedians, orbits, files, headers, cubeStd, cubeMode)
   
    
    
