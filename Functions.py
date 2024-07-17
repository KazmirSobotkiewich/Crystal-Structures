import numpy as np
from itertools import product
import matplotlib.pyplot as plt

### Define all constants ###

## Lattice parameters (in Angstroms)

# Monoclinic M1
M1a = 5.743
M1b = 4.517
M1c = 5.375

# Rutile R
Ra = 4.517
Rb = 4.517
Rc = 2.872

## Lattice vectors

# Monoclinic M1
M1a1 = np.array([0, 0, -5.742]) # a1 = a
M1a2 = np.array([-4.517, 0, 0]) # a2 = b
M1a3 = np.array([0, 4.526, -2.900]) # a3 = c

# Rutile R
Ra1 = np.array([4.56, 0, 0]) # a1 = a
Ra2 = np.array([0,  4.56, 0]) # a2 = b
Ra3 = np.array([0, 0,  2.86]) # a3 = c

## Reciprocal lattice vectors Function

def reciprocal_lattice_vectors(a1, a2, a3):
    """
    Input the 3 primitive lattice vectors as arrays and this will
    return the 3 primitive reciprocal lattice vectors as arrays.
    """
    b1 = (2*np.pi*np.cross(a2, a3))/(np.dot(a1, np.cross(a2, a3)))
    b2 = (2*np.pi*np.cross(a3, a1))/(np.dot(a2, np.cross(a3, a1)))
    b3 = (2*np.pi*np.cross(a1, a2))/(np.dot(a3, np.cross(a1, a2)))
    return b1, b2, b3

## Reciprocal lattice vectors

M1b1, M1b2, M1b3 = reciprocal_lattice_vectors(M1a1, M1a2, M1a3) # b1, b2 and b3
Rb1, Rb2, Rb3 = reciprocal_lattice_vectors(Ra1, Ra2, Ra3) # b1, b2 and b3

## Location of atoms for M1 (12 atoms)

# Calculated as weighted ratios of lattice vectors
M1r1 = 0.61*M1a1 + 0.19*M1a2 + 0.21*M1a3 # O1
M1r2 = 0.61*M1a1 + 0.31*M1a2 + 0.71*M1a3 # O1
M1r3 = 0.10*M1a1 + 0.21*M1a2 + 0.20*M1a3 # O1
M1r4 = 0.10*M1a1 + 0.29*M1a2 + 0.70*M1a3 # O1
M1r5 = 0.90*M1a1 + 0.71*M1a2 + 0.30*M1a3 # O2
M1r6 = 0.90*M1a1 + 0.79*M1a2 + 0.80*M1a3 # O2
M1r7 = 0.39*M1a1 + 0.69*M1a2 + 0.29*M1a3 # O2
M1r8 = 0.39*M1a1 + 0.81*M1a2 + 0.79*M1a3 # O2
M1r9 = 0.76*M1a1 + 0.03*M1a2 + 0.97*M1a3 # V
M1r10 = 0.24*M1a1 + 0.53*M1a2 + 0.53*M1a3 # V
M1r11 = 0.76*M1a1 + 0.48*M1a2 + 0.47*M1a3 # V
M1r12 = 0.24*M1a1 + 0.97*M1a2 + 0.03*M1a3 # V
# Put these into an array
M1r = np.array([M1r1, M1r2, M1r3, M1r4, M1r5, M1r6, M1r7, M1r8, M1r9, M1r10, M1r11, M1r12])

## Location of atoms for R (6 atoms)

# Calculated as weighted ratios of lattice vectors
Rr1 = 0.3001*Ra1 + 0.3001*Ra2 + 0.00*Ra3 # O
Rr2 = -0.3001*Ra1 - 0.3001*Ra2 + 0.00*Ra3 # O
Rr3 = 0.8001*Ra1 + 0.8001*Ra2 + 0.50*Ra3 # O
Rr4 = 0.1999*Ra1 + 0.1999*Ra2 + 0.50*Ra3 # O
Rr5 = 0.00*Ra1 + 0.00*Ra2 + 0.00*Ra3 # V (0,0,0)
Rr6 = 0.50*Ra1 + 0.50*Ra2 + 0.50*Ra3 # V (1/2,1/2,1/2)
# Put these into an array
Rr = np.array([Rr1, Rr2, Rr3, Rr4, Rr5, Rr6])

## Elastic atomic scattering factors of electrons for neutral atoms and s up to 2.0 Angstrom^-1

# Vanadium

Va = np.array([0.3876, 1.2750, 1.9109, 2.8314, 1.8979])
Vb = np.array([0.2967, 2.3780, 8.7981, 35.9528, 101.7201])

# Oxygen

Oa = np.array([0.0974, 0.2921, 0.6910, 0.6990, 0.2039])
Ob = np.array([0.2067, 1.3815, 4.6943, 12.7105, 32.4726])

### Define all functions ###

def atomic_form_factor_R(h, k, l):
    """
    Calculates the atomic form factors of O and V in the VO2 Rutile R phase,
    given the miller indices "h, k, l". Returns a 2D array containing the atomic
    form factor for Oxygen and then Vanadium.
    """
    G = h*Rb1 + k*Rb2 + l*Rb3
    s = ((G[0]**2 + G[1]**2 + G[2]**2)**(1/2))/(4*np.pi)
    n1, fO, n2, fV = 0, 0, 0, 0
    for i in Oa:
        fO += Oa[n1]*np.exp(-Ob[n1]*s**2)
        n1 += 1
    for i in Va:
        fV += Va[n2]*np.exp(-Vb[n2]*s**2)
        n2 += 1
    return(np.array([fO, fV]))
    
def atomic_form_factor_M1(h, k, l):
    """
    Calculates the atomic form factors of O and V in the VO2 Monoclinic M1 phase,
    given the miller indices "h, k, l". Returns a 2D array containing the atomic
    form factor for Oxygen and then Vanadium.
    """
    G = h*M1b1 + k*M1b2 + l*M1b3
    s = ((G[0]**2 + G[1]**2 + G[2]**2)**(1/2))/(4*np.pi)
    n1, fO, n2, fV = 0, 0, 0, 0
    for i in Oa:
        fO += Oa[n1]*np.exp(-Ob[n1]*s**2)
        n1 += 1
    for i in Va:
        fV += Va[n2]*np.exp(-Vb[n2]*s**2)
        n2 += 1
    return(np.array([fO, fV]))

def S_R(h, k, l):
    """
    Calculates the structure factor for VO2 in the Rutile R phase, given
    the miller indices "h, k, l". Returns a 2D array containing the real and 
    imaginary parts of the structure factor.
    """
    S_re, S_im, n = 0, 0, 0
    G = h*Rb1 + k*Rb2 + l*Rb3
    s = ((G[0]**2 + G[1]**2 + G[2]**2)**(1/2))/(4*np.pi)
    n1, fO, n2, fV = 0, 0, 0, 0
    for i in Oa:
        fO += Oa[n1]*np.exp(-Ob[n1]*s**2)
        n1 += 1
    for i in Va:
        fV += Va[n2]*np.exp(-Vb[n2]*s**2)
        n2 += 1
    f_R = [fO, fO, fO, fO, fV, fV]
    for i in Rr:
        x = -1*np.dot(G, Rr[n])
        S_re += np.cos(x)*f_R[n]
        S_im += np.sin(x)*f_R[n]
        n += 1
    return np.array([S_re, S_im])

def S_M1(h, k, l):
    """
    Calculates the structure factor for VO2 in the Monoclinic M1 phase, given
    the miller indices  "h, k, l". Returns a 2D array containing the real and 
    imaginary parts of the structure factor.
    """
    S_re, S_im, n = 0, 0, 0
    G = h*M1b1 + k*M1b2 + l*M1b3
    s = ((G[0]**2 + G[1]**2 + G[2]**2)**(1/2))/(4*np.pi)
    n1, fO, n2, fV = 0, 0, 0, 0
    for i in Oa:
        fO += Oa[n1]*np.e**(-Ob[n1]*s**2)
        n1 += 1
    for i in Va:
        fV += Va[n2]*np.e**(-Vb[n2]*s**2)
        n2 += 1
    f_M1 = [fO, fO, fO, fO, fO, fO, fO, fO, fV, fV, fV, fV]
    for i in M1r:
        x = -1*np.dot(G, M1r[n])
        S_re += np.cos(x)*f_M1[n]
        S_im += np.sin(x)*f_M1[n]
        n += 1
    return np.array([S_re, S_im])

def phase_M1(h,k,l):
    """
    Calculates the complex angle phase in the Monoclinic M1 transition phase 
    given the miller indices "h, k, l".
    """
    S_im = S_M1(h,k,l)[1] #Imaginary part of the Structure Factor.
    S_re = S_M1(h,k,l)[0] #Real part of the Structure Factor.
    return np.arctan2(S_im, S_re) #Return the complex angle

def phase_R(h,k,l):
    """
    Calculates the complex angle phase in the Rutile R transition phase 
    given the miller indices "h, k, l".
    """
    S_im = S_R(h,k,l)[1] #Imaginary part of the Structure Factor.
    S_re = S_R(h,k,l)[0] #Real part of the Structure Factor.
    return np.arctan2(S_im, S_re) #Return the complex angle

                
def miller(n):
    """
    Calculate all reflections from [-n,-n,-n] to [n,n,n]. Increasing "n", will
    increase the precision of the precision of any calculation using miller indices.
    """
    rdm = list(range(-n, n+1))
    return list(product(rdm,repeat = 3))

def electrostatic_potential_R(x,y, miller, a = (1,1,0), b = (0,0,1)):
    """
    Calculate the Electrostatic Potential for the Rutile R phase of VO2.
    The axes along which the cross section is sliced is given by "a" and "b".
    By default, we choose: x-axis = a + b and y-axis = c, because of the symmetry of VO2.
    """
    n, epot = 0, 0
    for i in miller:
        G = miller[n][0]*Rb1 + miller[n][1]*Rb2 + miller[n][2]*Rb3
        I = (S_R(*miller[n])[0])**2 + (S_R(*miller[n])[1])**2
        epot += (I**(1/2))*np.cos(phase_R(*miller[n]))*np.cos(x*np.dot(a,G) + y*np.dot(b,G))
        n += 1
    return epot

def electrostatic_potential_M1(x,y, miller, a = (1,1,0), b = (0,0,1)):
    """
    Calculate the Electrostatic Potential for the Monoclinic M1 phase of VO2.
    The axes along which the cross section is sliced is given by "a" and "b".
    By default, we choose: x-axis = a + b and y-axis = c, because of the symmetry of VO2.
    """
    n, epot = 0, 0
    for i in miller:
        G = miller[n][0]*M1b1 + miller[n][1]*M1b2 + miller[n][2]*M1b3
        I = (S_M1(*miller[n])[0])**2 + (S_M1(*miller[n])[1])**2
        epot += (I**(1/2))*np.cos(phase_M1(*miller[n]))*np.cos(x*np.dot(a,G) + y*np.dot(b,G))
        n += 1
    return epot

def Mag_G_M1(miller):
    """
    Returns a 1D list of the magnitudes of G for all miller indices for VO2 in the M1 phase.
    """
    magM1 = []
    for i in miller:
        x = M1b1*i[0]+M1b2*i[1]+M1b3*i[2] # G for M1
        y = ((x[0]**2 + x[1]**2 + x[2]**2)**(1/2)) # |G| for M1
        magM1.append(y)
    return magM1
        
def Mag_G_R(miller):
    """
    Returns a 1D list of the magnitudes of G for all miller indices for VO2 in the R phase.
    """
    magR = []
    for i in miller:
        x = Rb1*i[0]+Rb2*i[1]+Rb3*i[2] # G for M1
        y = ((x[0]**2 + x[1]**2 + x[2]**2)**(1/2)) # |G| for M1
        magR.append(y)
    return magR

def I_M1(h,k,l):
    """
    Returns the intensity for given miller indices in M1 phase of VO2.
    """
    x = S_M1(h,k,l)
    return x[0]**2 + x[1]**2

def I_R(h,k,l):
    x = S_R(h,k,l)
    I = x[0]**2 + x[1]**2
    return I

### Test code ###
# =============================================================================
# if __name__ == "__main__":
#     """Calculate the miller indices"""
#     miller = Miller(4)
# 
#     """Plot the Electrostatic Potential Maps"""
#     x = np.arange(0.0, 10, 0.1) # First and second argument set the range over which to plot
#     y = np.arange(0.0, 10, 0.1) # Third argument sets the width of the grid, therefore make this small to make the # of points large
#     X, Y = np.meshgrid(x, y)
# 
#     # Uncomment the function you wish to plot, i.e. Z
#     Z = Epot_M1(X, Y, miller)
# 
#     im = plt.imshow(Z, cmap=plt.cm.RdBu_r, extent=(0.0, 10, 0.0, 10))  # This labels the range of the plot correctly
#     plt.colorbar(im).set_label('Electrostatic Potential')  
#     plt.title('Cross Section of Electrostatic Potential for R Phase')
#     plt.xlabel('$a_r + b_r [\AA]$')
#     plt.ylabel('$c_r [\AA]$')
#     plt.rcParams.update({'font.size': 22})
#     plt.figure()
#     plt.show()
# =============================================================================
    