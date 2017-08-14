import os
import sys
import re
from matplotlib import pyplot
from pylab import genfromtxt

path = "./"
files = [f for f in os.listdir(path) if f.endswith('.dat')]
files.sort()
for file in files:
    mat = genfromtxt(file)
    pyplot.plot(mat[:,0],mat[:,1],label=file)
    pyplot.grid()
    pyplot.legend()
    pyplot.title("Solution Heat equation")
    pyplot.xlabel("x")
    pyplot.ylabel("u")
    pre, ext = os.path.splitext(file)
    
pyplot.savefig("solution.pdf")
