def process(filename):
    data = numpy.loadtxt(fname=filename, delimiter=',') 
    ave_brightness = data.mean(axis=0)
    processed_data = data / ave_brightness
    return processed_data
def plot(processed_data):
    image  = matplotlib.pyplot.imshow(processed_data) 
    matplotlib.pyplot.show(image) 
def detect_variables(processed_data):
    # variation of each star
    deviations = processed_data.std(axis=1)
    # calculate the 'mean' of the variations 
    mean_deviation = deviations.mean()
    std_deviations = deviations.std()
    for star_deviation in deviations:
        if (star_deviation - mean_deviation > 3.0*std_deviations):
            print('variable star in this data')
def detect_variables(processed_data):
    # variation of each star
    deviations = processed_data.std(axis=1)
    # calculate the 'mean' of the variations 
    mean_deviation = deviations.mean()
    std_deviations = deviations.std()
    for star_deviation in deviations:
        if (star_deviation - mean_deviation > 3.0*std_deviations):
            print('variable star in this data')
import numpy
import matplotlib.pyplot
import glob

filenames = glob.glob('data/*.csv')
for f in filenames:
    print(f)
    processed_data = process(f)
    plot(processed_data)
    detect_variables(processed_data)