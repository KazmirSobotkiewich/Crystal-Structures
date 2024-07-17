# This code plots the electrostatic potential map for VO2 in the Monoclinic M1 phase.
# The number of miller indices can be increased by changing "n" in "miller = f.Miller(n)".
# This will increase the precision of the plots, but can make certain structural features disappear when "n" is too large.

""" Import libraries """

import numpy as np
import matplotlib.pyplot as plt
import Functions as f

""" Get the miller indices """

miller_indices = f.miller(6)

""" Set up a grid of points to be inputted in the electrostatic potential map function """

x = np.arange(0.0, 12.5, 0.1) # First and second argument set the range over which to plot
y = np.arange(0.0, 12.5, 0.1) # Third argument sets the width of the grid, decrease this to increase the # of points
X, Y = np.meshgrid(x, y)

""" Calculate the electrostatic potential for all points in the grid """

Z = f.electrostatic_potential_M1(X, Y, miller_indices)

""" Plot the Electrostatic Potential Map for VO2 in the Monoclinic M1 phase """

im = plt.imshow(Z, cmap=plt.cm.RdBu_r, extent=(0.0, 12.5, 0.0, 12.5))  # This labels the range of the plot correctly
plt.rcParams.update({'font.size': 16})
plt.colorbar(im).set_label('Electrostatic Potential (A.U.)')  
plt.title('Cross Section of $VO_2$ in $M_1$ Phase')
plt.xlabel('Along $a_R + b_R\:(\AA)$')
plt.ylabel('Along $c_R\:(\AA)$')
#plt.savefig('Cross-Section-VO2-M1.png') # If we want to automatically save the figure
plt.figure()
plt.show()