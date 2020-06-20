README

Input file = "2002QF15moginputs.txt"

mog function takes a three-line file (as a numpy array) as a parameter and calculates the position, velocity, and range vectors, as well as the orbital elements.

mogpermute function takes a four-line file and slices it into two three-line files, each representing a permutation of the original four-line file.
mogpermute function calculates the orbital elements and vectors of each of the three-line files by calling the mog function

No user inputs required while running

The precess time is hard-coded and allocated to memory as the variable t_0

KimPOD.py is the main file, KimPfunctions.py contains most of the functions used in the main, babyod2.y contains the code which converts orbital elements into orbital vectors.
Only KimPOD.py needs to be opened because it imports and calls functions from the other two files as necessary.

Code runs for all possible real, positive roots of lagrange which return positive values for the range vectors.

Code corrects for light travel time.

All angles in degrees, distances in AU, and velocities in AU/day in the output.

Monte Carlo included. 

montecarlo.py runs the monte carlo simulation and writes all of the data to a csv file

montetest.py passes in the csv file created by montecarlo.py and plots it in matplotlib. Only first orbital element distribution done so far.

Disregard mogfunction. It is only used in the montecarlo file.