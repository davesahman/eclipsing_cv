import numpy
import matplotlib
import matplotlib.pyplot
data = numpy.loadtxt(fname='data/star_data_01.csv', delimiter=',')
rows = data.shape[0]
cols = data.shape[1]
a = numpy.zeros(rows)
for i in numpy.arange(0,rows):
    a[i] = (data[i:].std())
print(a[:])    
matplotlib.pyplot.plot(a[:])
matplotlib.pyplot.show()