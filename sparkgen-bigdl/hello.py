import numpy as np
from numpy import genfromtxt
from scipy.stats import f, t
import math

from LaTeXPrinter import create_table
from dataCollector import collector

raw_data = genfromtxt('data/24-10-2019_1324.csv', delimiter=',')
titles, names, parameters, DATA = collector(raw_data)

DATA = DATA[:][:10][:]
print DATA