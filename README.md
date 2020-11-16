This work continues previous work done by David Gomez and Majd Mayyasi. Goal: characterize the particle sensitivity of the MAVEN IUVS instrument throughout the mission, 
identify (using IUVS data) orbits with single-particle or particle-barrage interactions with the IUVS detector, and correlate these detections with known SEP/ICME events.

----------------------------------------- Mission Long Analysis ----------------------------------------- 

main function file: "FULLMISSION_Bens_Functions_allDark"

python notebooks: all "FULLMISSION_fuvdark_{}" (mission long dark frame analysis for all observing modes), and "ShapesOfFiles" (tells' the shape of every .fits file 
for all observing modes)

- These files will take a long time to run and require a suffieicnt processor (min 3.0 GHz) and ram (min 16 GB)

- the main function files executes functions in a specific order, displays postage stamp images, histograms, and outlier plots.

- Run the python notebooks with the proper directory to the folder containing the fuvdark images.

- There should be no need to re-run the "ShapesOfFiles" notebooks as I have already classified each observing mode and each file using a batch method, which is 
  visible in the "FULLMISSION" notebooks.  The "FULLMISSION_{}" notebooks are distinguished by observing mode, and loop through all specified shapes for all
  orbits (between 0-12000), searching every file from the MAVEN mission with the correct matching orbit and data shape.  A data cube is created for each shape in every 
  observing mode. The images are all plotted using the same colormap.  The 6-sigma outliers foor each image are calculated using the mode and standard deviation of the 
  Data Number histogram, and the outlier plot is shown in the notebooks.
  
 - Due to GitHub size contraints, the output from each cell has been cleared.  
 - Please contact me for a detailed "walkthrough" of the function files and python notebooks, as I will be able to give a better verbal presentation of these files.



----------------------------------------- December 2014 SEP Event Analysis -----------------------------------------

main function file: "BensFunctions, BensFunction_allLight"

python notebooks: all "25dec2014_fuvdark_{}" (dark frame analysis), "25dec2014_fuv" (for light frame analysis), and "{}_ShapesOfFiles" 



----------------------------------------- July 2017 SEP Event Analysis -----------------------------------------

main function file: "BensFunctions, BensFunction_allLight"

python notebooks: all "23july2017_fuvdark_{}" (dark frame analysis), "23july2017_fuv" (for light frame analysis), and "{}_ShapesOfFiles" 



----------------------------------------- September 2017 SEP Event Analysis -----------------------------------------

main function file: "BensFunctions, BensFunction_allLight"

python notebooks: all "13sept2017_fuvdark_{}" (dark frame analysis), "13sept2017_fuv" (for light frame analysis), and "{}_ShapesOfFiles" 



------------------------------------------------------------------------------------------------------------------------------------

Notebooks showing localized impact analysis not uploaded, and there are many more notebooks for each SEP event for the remaining observing modes.  Contact me with any issues, 
questions, or requests for more notebooks at Ben.Johnston@lasp.colorado.edu
