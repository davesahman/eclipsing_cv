import glob
import numpy
import matplotlib.pyplot

filenames = glob.glob('data/*.csv')
filenames = filenames[0:4] 
for f in filenames:
    print (f)
    
    data = numpy.loadtxt(fname=f, delimiter=',') # load in the data
    # calculate the average brightness over all stars (rows)
    ave_brightness = data.mean(axis=0) 
    # divide by the average brightness
    processed_data = data/ave_brightness
    
    image  = matplotlib.pyplot.imshow(processed_data) 
    matplotlib.pyplot.show(image) 
    