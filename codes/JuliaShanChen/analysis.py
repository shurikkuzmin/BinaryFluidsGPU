import pylab
import numpy
import os

print os.getcwd()
profile=numpy.loadtxt("height00220.dat")

pylab.plot(profile[128, :])
pylab.show();
