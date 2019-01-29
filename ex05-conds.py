import numpy

#load data
data = numpy.loadtxt(fname='data/star_data_04.csv', delimiter=',')

# calculate the average brightness at each time, over all stars (rows)
ave_brightness = data.mean(axis=0)

# divide by the average brightness
processed_data = data/ave_brightness

# find the standard deviation of each star
deviations = processed_data.std(axis=1)
print(deviations)
mean_deviation = deviations.mean()
std_deviations = deviations.std()
for star_deviation in deviations:
    if (star_deviation - mean_deviation > 3.0*std_deviations):
        print('variable star in this data')
import glob

filenames = glob.glob('data/*.csv')
for f in filenames:
    print (f)
    
    data = numpy.loadtxt(fname=f, delimiter=',') # load in the data

    # calculate the average brightness over all stars (rows)
    ave_brightness = data.mean(axis=0) 

    # divide by the average brightness
    processed_data = data/ave_brightness

    # standard deviation of each star
    deviations = processed_data.std(axis=1)
    
    mean_deviation = deviations.mean()
    std_deviations = deviations.std()
    for star_deviation in deviations:
        if (star_deviation - mean_deviation > 3.0*std_deviations):
            print('variable star in this data')
