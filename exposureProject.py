import numpy as np
import PgmPpmFormatter as ppf
import math
import sys
import os, sys

######################################
## CODE ALGORITHM BY KIIBATI ADEOJO
## BUG FIX RESEARCH BY TRIET NGUYEN
######################################

def main():
    filename = ""
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print("Please provide the name of a file in the form image.pgm.")
        return

    data = ppf.readImageFile(filename)
    theImageArray = data[0]
    imageWidth, imageHeight = data[1]
    maxpixel = data[2]
    stringRC = ""+ str(imageHeight) + " " + str(imageWidth) + "\n"
    stringMax = str(maxpixel) + "\n"
    differenceThreshold = 4
    newList = []
    bigList = []

    binSpace = maxpixel/5 # 5 because we have 5 different levels of histogram exposure we are considering
    histogram = []

    hist,bins = np.histogram(theImageArray, bins = [0,51,102,153,204,255]) 
    overExposed = True
    if ( (hist[4] > hist[0]) and ( (hist[4]/hist[0]) >= differenceThreshold) ):
        overExposed = True
        print ("This image is Overexposed")
    elif ( (hist[4] < hist[0]) and ( (hist[0]/hist[4]) >= differenceThreshold) ):
        overExposed = False
        print ("This image is Underexposed")

    for i in range(0, len(theImageArray)):

        if (overExposed):
            if( theImageArray[i] >= 204):
                newList.append(255)
            else:
                newList.append(0)
        else:
            if( theImageArray[i] <= 51):
                newList.append(0)
            else:
                newList.append(255)
        
        bigList.append(newList)
        newList = []
    npArray = np.array(bigList) 

    arrangedArray = np.reshape( npArray, (imageWidth, imageHeight * 3))

    invertedFileName = filename.replace('.ppm', '') + "_Inverted.ppm"
    ppf.createFile(invertedFileName, "P3\n", stringRC, stringMax, arrangedArray)
    
################################################################
## CONVERTING A PPM FILE TO PGM FILE
################################################################

    pgmVersion = invertedFileName.replace('.ppm', '.pgm')
    newPgmArray = []
    for i in range(0, len(npArray)):
        if ( (i+1) % 3 == 0 ):
            greyscaleCalc =  (0.1140*npArray[i]) + (0.5870*npArray[i-1]) + (0.29890*npArray[i-2]) 
            # print greyscaleCalc
            value = int(greyscaleCalc)
            newPgmArray.append(int(value))

    fileArray = np.array(newPgmArray)
    ppf.createFile(pgmVersion, "P2\n", stringRC, stringMax, fileArray)

# ##############################################################
# #Smooth out image twice
# #Take any number between 51 and 204, and change it to white
# #If pixel in pgm is white then that pixel is black in the ppm image

    theImageArrayEdit = ppf.smoothImage(1.4, invertedFileName, "P2", stringRC, stringMax, fileArray, (imageWidth, imageHeight) )
    #Take any number between 51 and 204, and change it to white
    # OR
    #Take any number > 0 and < 51 , and change it to black 
    for i in range(0, len(theImageArrayEdit)):
        if(overExposed):
            if( ( int(theImageArrayEdit[i]) < 204 and int(theImageArrayEdit[i]) > 51) ):
                theImageArrayEdit[i] = 255
        else:
            if( ( int(theImageArrayEdit[i]) < 51 and int(theImageArrayEdit[i]) > 0) ):
                theImageArrayEdit[i] = 0            

    # Turn pgm to ppm
    # Turn overexposed or underexposed section of the picture into a uniform color
    # Turn other pixels into their original colored pixel
    for i in range(0, len(theImageArrayEdit)):
        if(overExposed):
            if( theImageArrayEdit[i] >= 204 ):
                theImageArray[(i+1) *3 - 1] = 0
                theImageArray[(i+1) *3 - 2] = 0
                theImageArray[(i+1) *3 - 3] = 0
        else:
            if( theImageArrayEdit[i] <= 51 ):
                theImageArray[(i+1) *3 - 1] = 255
                theImageArray[(i+1) *3 - 2] = 255
                theImageArray[(i+1) *3 - 3] = 255
 
    ppf.createFile(filename.replace('.ppm', '')+"_finalExposure.ppm", "P3\n", stringRC, stringMax, theImageArray)

if __name__== "__main__":
    main()