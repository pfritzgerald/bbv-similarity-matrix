#!/usr/bin/env python

#Created by Spencer Hance and Trevor Gale on January 18th 2015
#Northeastern University Computer Architecture Research Group
#Licensed under MIT License

import sys
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import re
import random
from scipy.misc import comb
import argparse
import warnings


def parseBBV(input_filename):
    """Parses a Basic Block Vector and converts data
    into a Numpy array
    """
    with open(input_filename, 'r') as f:
        input_list = []
        # Opens file into a list
        for line in f.readlines():
            # Ignores BBV comments, which are any line that starts with a "#"
            if not line.strip().startswith('#'):  
                input_list.append(line.split())   
    
    # Removes empty list elements
    input_list = filter(None, input_list)
    num_intervals = len(input_list)
    
    # Determines the total number of basic blocks
    max_list = []
    for line in input_list:
        for j in range(0, len(line)):
            m = re.search(":(\d+):(\d+)", line[j])
            max_list.append(int(m.groups()[0]))
    num_bb = max(max_list)
    
    # Initializes array and adds basic block data
    bbv_array = np.zeros((num_intervals, num_bb))
    for i in range(0, num_intervals):
        for j in range(0, len(input_list[i])):
            m = re.search(":(\d+):(\d+)", input_list[i][j])
            bbv_array[i, int(m.groups()[0])-1] = int(m.groups()[1])
    
    # Update user on current progress
    print 'Parsing Completed\n'
    return bbv_array 


def reduceArray(bbv_array):
    """Takes in numpy array of bbv vectors and reduces dimensions to 15.
    Returns the reduced array
    """
    # Initializes an array with the same number of rows 
    # as the BBV numpy array and 15 columns
    random_array = np.zeros((bbv_array.shape[1], 15))
    
    # Fills the array with a random float between -1 and 1
    for i in range(0, random_array.shape[0]):
        for j in range(0, random_array.shape[1]):
            random_array[i, j] = random.uniform(-1,1)

    # Takes the dot product of the two arrays to reduce
    # the total dimensions to 15
    reduced_array = np.dot(bbv_array, random_array)
    return reduced_array


def mDistCompute(a, b):
    """Takes in two 1D arrays and computes sum of
    manhattan distances.  This function is an inner function of mDist()
    """
    # Initialize the sum value
    sum_dist = 0
    
    # Both arrays must be of of the same length
    length = len(a) 
    
    # Compute sum of differences
    for i in range(0, length):
         sum_dist += abs(a[i]- b[i])   
    return sum_dist


def mDist(bbv_array):
    """Takes in bbv array and calls mDistCompute to compute
    manhattan distance between the vectors. Returns an 
    array with differences.
    """
    # Determines the size of the array
    mDist_length = bbv_array.shape[0]
    
    # Initializes a new array to store distance values
    mDist_array = np.zeros((mDist_length, mDist_length))
    
    # Determines total number of steps for progress bar
    total_steps = float(comb(mDist_length, 2, exact=True))
    
    # Initializes step counter for progress bar
    step = 0

    # Compute distances by using mDistCompute() for each comparison
    print 'Computing Manhattan Distances'
    for i in range(0, mDist_length):
        for j in range(1+i, mDist_length):
            sum_dist = mDistCompute(bbv_array[i], bbv_array[j])
            mDist_array[i, j] = sum_dist
        
        # Calculations for progress counter    
        step += len(range(1+i, mDist_length))   
        sys.stdout.write('\r')
        sys.stdout.write('Completion: ' + \
                str(int(round((step/total_steps)*100))) + '%')
        sys.stdout.flush()
    print '\n'
    return mDist_array


def normMatrix(mDist_values):
    """Takes in array of manhattan distance values and
    returns the array normalized to the maximum value
    """
    #Renames input to norm_array
    norm_array = mDist_values

    #Determines the largest distance to normalize to
    max_val = max(max(l) for l in norm_array)
    
    # Update user on current progress
    print 'Normalizing Matrix\n'
    
    #Replaces every value with the new normalized value
    for i in range(0, norm_array.shape[0]):
        for j in range(0, norm_array.shape[1]):
            norm_array[i, j] /= max_val
    return norm_array


def plotNormData(norm_values, show=True):
    """Takes in  normalized values and plots 
    the data
    """
    # Initialize lists for plt.scatter
    x, y, colors = [], [], []

    # Determines the height of the array for the graph's Y-Value
    yval = norm_values.shape[0]
    
    # The size of each point
    # Dividing by 4.5 usually provides enough granularity, however this should 
    # be adjusted if a different resolution requirement is needed
    SIZE = yval/4.5
    
    # Update user on current progress
    print 'Plotting Norm Data\n'
    
    #Adds data to x, y, and colors lists
    for i in range(0, yval):
        for j in range(i, yval):
            x.append(j)
            y.append(i)
            colors.append(norm_values[i,j])    
    
    #Plots data with gray colormap and aligns both axes to 0
    plt.scatter(x, y, c = colors, cmap=cm.gray, s = SIZE)
    plt.xlim(0)
    plt.ylim(0)

    #Inverts y axis to show similarity accurately
    plt.gca().invert_yaxis()
    if show == True:
         plt.show()
    

def commandParser(): 
    """Uses argparse module to parse command line options
    """
    parser = argparse.ArgumentParser(description='Similarity Matrix Generator \
            for Basic Block Vectors')
    parser.add_argument('-i',dest='filename', required=True, help='input BBV file',
            metavar='file')
    parser.add_argument('-s','--simmatrix', help='Create and display a similarity matrix' ,
            action='store_true')
    parser.add_argument('-dr','--do-not-reduce',
    help='Do not reduce input matrix for similarity matrix', action='store_true')
    args = parser.parse_args()
    
    if not args.filename:
        print 'Error: Not enough input arguments'
    if args.do_not_reduce:
        print 'Starting Similarity Matrix Process (with unreduced array)\n'
        plotNormData(normMatrix(mDist(parseBBV(args.filename))))    
    else:
        print 'Starting Similarity Matrix Process\n'
        plotNormData(normMatrix(mDist(reduceArray(parseBBV(args.filename)))))
    

def main():
    """Main Function"""
    commandParser()


if __name__ == '__main__':
    main()

