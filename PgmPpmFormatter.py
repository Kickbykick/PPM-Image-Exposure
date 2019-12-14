import numpy as np

######################################
## PGM/PPM FORMATTER BY KIIBATI ADEOJO
######################################
GAUSS5x5 = np.array([[1, 3, 6, 3, 1], 
                     [3, 15, 25, 15, 3], 
                     [6, 25, 41, 25, 6], 
                     [3, 15, 25, 15, 3], 
                     [1, 3, 6, 3, 1]])/253
KERNEL1 = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]])
KERNEL2 = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]])
DarkRange = 50

def readImageFile(fileName):
    with open(fileName) as pgmFile:
        allLines = pgmFile.readlines()
    
    # assert (allLines[0].strip() == 'P3' || allLines[0].strip() == 'P2')

    #Ignoring Comments
    for str_line in list(allLines): #Turn Lines into a list and check for any string that starts with # and remove it
        if str_line[0] == '#':
            allLines.remove(str_line)

    #list
    data_list = []
    for line in allLines[1:]:
        data_list.extend([int(i) for i in line.split()])

    resultMulti = (np.array(data_list[3:]),(data_list[1],data_list[0]),data_list[2])
    # newdata = np.reshape(resultMulti[0],resultMulti[1]), (data_list[0],data_list[1]), data_list[2]

    return resultMulti

def smoothImage(sigma, basename, fileType, str_fileRowsCols, fileMaxPixel, fileArray, rowsColsTuple):
    maxpixel = int(fileMaxPixel)
    array = np.reshape( fileArray, (rowsColsTuple[0], rowsColsTuple[1]) )
    numrows, numcols = array.shape
    
    size = int(5) // 2
    x, y = np.mgrid[-size:size+1, -size:size+1]
    normal = 1 / (2.0 * np.pi * sigma**2)
    g =  np.exp(-((x**2 + y**2) / (2.0*sigma**2))) * normal

    smooth_array2d = convolve2D(array, g)
    smooth_name2d = basename.replace('.ppm', '') + '-2Dsmooth.pgm'
    writeimage(smooth_name2d, fileType, fileMaxPixel, smooth_array2d)

    return np.array( smooth_array2d ).flatten()


def createFile(fileName, fileType, fileRowsCols, fileMaxPixel, fileArray):
    theFile = open(fileName,"w")
    theFile.write(fileType)
    theFile.write(fileRowsCols)
    theFile.write(fileMaxPixel)
    np.savetxt(theFile, fileArray, fmt="%s")  

#######################################################################################
#### Codes from Professor J, Vaughan
#######################################################################################
def convolve2D(array, kernel):
    """
    This function convolves an array of pixel values with a 2D kernel.
    The kernel should have dimensions (2k+1) x (2k+1).
    
    Position (i,j) in the array is always matched with the center of the kernel,
    position (k,k).
    
    If i and j are at least a distance k from the boundaries of the array, then
    the entire kernel can be matched to a (2k+1) x (2k+1) slice of the array
    centered on position (i,j).  We take this slice, multiply it component-wise
    by the kernel, and sum the results.  This number becomes entry (i,j) in
    the convolved image.
        
    If position (i,j) is such that a (2k+1) x (2k+1) slice of the array
    centered on this position would exceed the boundaries of the array, then we
    adjust the slice to stop at the boundaries, and we truncate the kernel to 
    match.
    
    For example, if (i,j) = (0,0), then the slice of the array would be
    array[0:k+1, 0:k+1], and the corresponding slice of the kernel would be
    kernel[k:2k+1, k:2k+1].
    
    If all of the entries in the kernel are positive, then we assume that the
    entries in the kernel sum to 1.  In this case, if the kernel is truncated,
    then we sum the entries in the truncated kernel and rescale by the
    reciprocal of this number.  This way, the entries in the truncated kernel
    also sum to 1.
    
    If at least one entry in the kernel is not positive, then we simply 
    truncate the kernel and do not rescale it.
    """

    # Decide whether to rescale the kernel after truncating
    if np.amax(kernel) > 0:
        renormalize = True
    else:
        renormalize = False

    # Assumes that the dimensions of the kernel have the form (2k+1) x (2k+1)
    kernel_size = kernel.shape[0]
    margin = int((kernel_size-1)/2)
    
    # Get the dimensions of the pixel array
    numrows, numcols = array.shape

    # The result has the same dimensions as the pixel array
    conv_array = np.zeros((numrows, numcols))

    # Iterate over the pixel array to perform the convolution
    for ai in range(numrows):
        # First indices for the slice of the pixel array
        pu=max(ai-margin,0)
        pd=min(ai+margin+1,numrows)
        # First indices for the slice of the kernel
        ku=max(0,margin-ai)
        kd=kernel_size-max(ai+margin+1-numrows,0)

        for aj in range(numcols):            
            # Second indices for the slice of the pixel array
            pl=max(aj-margin,0)
            pr=min(aj+margin+1,numcols)            
            # Second indices for the slice of the kernel
            kl=max(0,margin-aj)
            kr=kernel_size-max(aj+margin+1-numcols,0)

            if kl == 0 and ku == 0 and kr == kernel_size and kd == kernel_size:
                truncated = False
            else:
                truncated = True

            # Find the sum of the truncated kernel, if neeed
            if not truncated or not renormalize:
                kernel_sum = 1
            else:
                kernel_sum = np.sum(kernel[ku:kd,kl:kr])
            
            # Multiply the slice of the array by the slice of the kernel
            # component-wise, take the sum, renormalize if appropriate
            conv_array[ai,aj] =  np.sum(array[pu:pd,pl:pr]*kernel[ku:kd,kl:kr])/kernel_sum
            
    return conv_array

def readimage(longname):
    """
    The header of a PGM file has the form

    P2
    20 30
    255

    where P2 indicates a PGM file with ASCII encoding,
    and the next three numbers are the number of columns in the image,
    the number of rows in the image, and the maximum pixel value.
    
    Comments are indicated by # and could occur anywhere in the file.
    
    The remaining entries are pixel values.

    This function returns the file type, the maximum pixel value, 
    and the pixel data in the form of a NumPy array whose dimensions 
    are the number of rows and number of columns in the image.
    """

    # dummy values that are returned if the file cannot be opened
    filetype = 'none'
    maxpixel = 0
    array = np.array([0])

    try:
        image_file = open(longname)
    except:
        print("Failed to open the file named " + longname + ".")

   
    """
    The first entry in a pgm file that is not a comment is the file type,
    which is saved as a string.  The remaining entries are all saved as integers.
    """
        
    # list to hold all integer entries in the file
    longlist = []

    firstword = True
    for line in image_file:
        words = line.split()
        wi = 0
        comment = False
        while not comment and wi < len(words):
            word = words[wi]
            wi += 1
            if not word.startswith('#') and not firstword:
                # an entry that is not part of a comment and is not 
                # the first word is an integer
                longlist.append(int(word))
            elif word.startswith('#'):
                # this is a comment
                # drop the rest of the line
                comment = True
            elif firstword:
                # the first word that is not a comment is the file type
                filetype = word
                firstword = False

    image_file.close()
    
    if filetype == 'P2':
        """
        After the file type, the first three integers in a pgm file are the
        number of columns in the image, the number of rows in the image, and
        the maximum pixel value.
        """
            
        numcols = longlist[0] # number of columns in the image
        numrows = longlist[1] # number of rows in the image
        maxpixel = longlist[2] # maximum pixel value
        
        """
        The remaining integers are pixel values.
        Arrange them in a NumPy ndarray using the dimensions of the image that 
        were read from the header
        """
        array = np.reshape(np.array(longlist[3:]), (numrows, numcols))
        
    elif filetype != 'none':
        # for the moment, only P2 files are supported
        print(filetype + " is not a recognized file type.")
    
    return filetype, maxpixel, array

def writeimage(longname, filetype, maxpixel, array):
    """
    This function receives a file name, the file type, the maximum 
    pixel value, and an array of pixel data.
    It writes a PGM file with the given name, header and pixel data.
    The pixel data must be in a NumPy ndarray whose dimensions are the 
    number of rows and number of columns in the image.  Entries in this 
    array are rounded and cast to integers before they are written to the file.
    """

    try:
        image_file = open(longname,'w')
    except:
        print("Failed to write to the file named " + longname + ".")
        return
        
    if filetype == "P2":
        # obtain the dimensions of the image from the shape of the array
        numrows, numcols = array.shape

        """
        The header of a PGM file has the form
    
        P2
        20 30
        255
    
        where P2 indicates a PGM file with ASCII encoding,
        and the next three numbers are the number of columns in the image,
        the number of rows in the image, and the maximum pixel value.
        
        All subsequent numbers are pixel values.
        """
    
        # write the file header
        image_file.write(filetype+'\n')
        image_file.write("{0} {1}\n".format(numcols,numrows))
        image_file.write("{}\n".format(int(maxpixel)))

        """
        Values in the array of pixel data are rounded and cast to integers before
        they are printed to the file.
        """
        for i in range(numrows):
            for j in range(numcols):
                image_file.write("{} ".format(int(round(array[i,j]))))
            image_file.write('\n')
                    
    else:
        # for the moment, only P2 files are supported
        print("This was not a P2 file.  No result was printed.")

    image_file.close()
    return

if __name__== "__main__":
    main()
