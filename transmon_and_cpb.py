# -*- coding: utf-8 -*-
"""
Python version: 3.8.5
IDE: Spyder 4.2.0

Originally crated on July-November 2018 by Joaquín Fernández Rossier.

Modified and extended by Arturo Mena López on 2021

H = 4Ec*(-i d/dphi  - ng)^2 - Ej*cos(phi)
"""


import numpy as np
import scipy.linalg as spy
import matplotlib.pyplot as plt



def josephson(phi, Ej):

	"""Josephson cosine potential."""
	
	return - Ej * np.cos(phi)



def ham(phimax, N, potfun, param, ng, Ej, Ec):
    
    '''
    Description:
        Discretize the hamiltonian in the line (-xmax, xmax) 
        
    INPUTS
        * phimax is a real number. It goes in radians.
        * N  (integer) is the number of points in the grid.
          Therefore, there are N-1 segments.
        * potfun is a function f(x,param), to be provided externally, in this case it will be the Josephson potential.
        * param is a set of parameters for the potential function (Ej in out case). 
        * ng is the effective offset charge of the transmon.
        * Ej is the Josephson energy.
        * Ec is the charge energy.
        
    OUTPUT: 
        A dimension N numpy array (matrix) with the discretized Hamiltonian.
    '''
    
    dphi=2.*phimax/float(N-1) #size of step
        
    epsilon_a = -1/dphi**2 + 1j*ng/dphi
    epsilon_b = -1/dphi**2 - 1j*ng/dphi
    
    # Create a square matrix of size (N,N) filled with 0.0j complex values.
    mat = np.full((N,N), 0.0j)
       
    # We now fill the diagonal
    phi=-phimax
    for i in np.arange(N):
        # potfun -> -Ej * np.cos(psi) # Josephson
        mat[i,i] = 2/dphi**2 + ng**2 + potfun(phi,param)/(4*Ec)
        phi=phi+dphi
        
    # We now fill the positive co-diagonal, the (i,i+1) terms.
    for i in range(N-1):
        mat[i,i+1]= epsilon_a

    # We now fill the negative co-diagonal, the (i+1,i) terms.
    for i in range(N-1):
        mat[i+1,i]= epsilon_b
        
    # Finally we fill the two corners outside the diagonal.
    mat[0,N-1] = epsilon_b # Upper right corner.
    mat[N-1,0] = epsilon_a # Lower left corner.
    
         
    return mat 



def plotstates(mat, kval, phimax, N, potfun, param, Ej, Ec):
    '''
    INPUTS:
    * mat is a  hamiltonian matrix, obtained with HAM.
    * kval is the number of states that are going to be plotted
    * xmax and N define the grid (same as in HAM to define mat)
    * potfun: pofun(x,param) is a used-defined function that defines the potential

    OUTPUT
    a figure showing V(x), the energy levels and the wave functions.
    '''

    dphi = 2.*phimax/float(N-1) #size of step
    
    # Create a list corresponding to the discreete phase space (philist) and 
    # a list with the values of the potential (potlist).
    
    potlist = []
    philist = []
    phi = -phimax
    
    for k in range(N):
        philist.append(phi)
        potlist.append(potfun(phi,param))
        phi = phi+dphi
    
    #Create the figure and its axis.
    fig = plt.figure(figsize=(7,5), constrained_layout=True)
    ax = fig.add_subplot(111)
    # plt.rc('text', usetex=True) # Use LaTeX font for text.
    
    # Plot the potential divided by Ej.
    ax.plot(philist, np.array(potlist)/Ej, color="black", lw=3, alpha=0.6, zorder=50)
    
    # Get the diagonalized hamiltonian matrix (ener) its eigenfunctions (wave).
    # Remember that mat is the tridiagonal hamiltonian matrix.
    # ener has the form E/(4Ec), with E being the energy of the state and 
    # Ec the charge energy.
    ener, wave = spy.eigh(mat) 
    
    # We modify ener so we get E/Ej instead of E/(4Ec)
    ener = ener * 4*Ec / Ej
    
    # For each i'th energy state, plot the wave function and its corresponding 
    # energy.
    
    enerlist=[]
    wavelist=[]
    
    for i in range(kval):
        wavelist.append(np.abs(wave[:,i])) # Add the i'th wave function modulus.
        # Plot the i'th state. We multiply wavelist[i] by 3 for better visualization.
        wave_plot_i, = ax.plot(philist, ener[i]+3*wavelist[i]) 
        color_i = wave_plot_i.get_color() # Get the color of the line plotted in wave_plot_i.
        
        # Create a list (subener) with the repeated energy value of the i'th 
        # state.
        # It will be used to plot a flat line indicating the energy value in 
        # the plot with its height.
        subener=[] 
        # Append the value for each phi in the space.
        for phi in philist:
            subener.append(ener[i]) 
            
        enerlist.append(subener)  # Add the flat energy line to a set of energies.
        ax.plot(philist, enerlist[i], color=color_i, linestyle="--", label='_nolegend_') # Plot the flat line.
        
        ax.fill_between(philist, ener[i]+3*wavelist[i], enerlist[i], 
                        alpha=0.3, 
                        color=color_i, 
                        label="State $|{0}\\rangle$".format(i))
    
    # The lower limit of the plot will be the minimum value of the potential.
    # ymin = 1.1*np.min(potlist)/Ej
    # ymax = 1.1*np.max(potlist)/Ej
    
    # Set the plot labels and limits
    ax.set_title("$E_J/E_C = {}$".format(round(Ej/Ec)), fontsize=16)
    ax.set_xlabel('$\\varphi$', fontsize=16)
    ax.set_ylabel('$E/E_{J}$', fontsize=16)
    ax.legend(loc="best", fontsize=12, framealpha=0.8)
    ax.set_xticks([-2*np.pi, -3*np.pi/2, -np.pi, -np.pi/2, 0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax.set_xticklabels(["$-2\\pi$",
                        "$-\\frac{3\\pi}{2}$",
                        "$-\\pi$",
                        "$-\\frac{\\pi}{2}$",
                        "$0$", 
                        "$\\frac{\\pi}{2}$", 
                        "$\\pi$", 
                        "$\\frac{3\\pi}{2}$", 
                        "$2\\pi$"], 
                         fontsize=14)
    # ax.set_ylim([ymin,ener[kval]+0.5])
    # ax.set_ylim([ymin,ymax])
    ax.set_xlim([-phimax, phimax])
    plt.show()

    return



def plot_Em_ng(ng_max, kval, Ej, Ec, num_ngs=50, N=150, phimax=1.5*np.pi):
    
    """
    Plot the relative eigenenergies Em/E01 in function of the 
    effective offset charge ng.
    
    INPUT:
        * ng_max: Upper bound of the ng-space. ng is the charge offset.
        * kval: Number of energy states that to be plotted.
        * Ej: Josephson energy.
        * Ec: Charging energy.
        * num_ngs: Number of points in the ng-space.
        * N: Parameter for ham(). Integer. Number of points on the phase grid.
        * phimax: For ham(). Upper limit for the phase space.
    
    OUTPUT:
        * Plot of Em/E01 in function of ng.     
    """
    
    ngs = np.linspace(-ng_max, ng_max, num_ngs)
    
    # Get the energy difference E01 between the ground state and the 1st 
    # excited state at the sweet spot ng=1/2.
    mat = ham(phimax, N, potfun=josephson, param=Ej, ng=1/2, Ej=Ej, Ec=Ec)
    ener = spy.eigh(mat, eigvals_only=True)  
    E0 = ener[0] * 4*Ec
    E1 = ener[1] * 4*Ec
    E01 = E1 - E0 
    
    enerlist = [] # rows will correspond to the different ng values
                  # columns will correspond to the different energy states.
    
    for ng in ngs:
        
        mat = ham(phimax, N, potfun=josephson, param=Ej, ng=ng, Ej=Ej, Ec=Ec)
        ener = spy.eigh(mat, eigvals_only=True) 
        # We modify ener so we get E instead of E/(4Ec)
        ener = ener * 4*Ec
        
        EmE01 = ener / E01 # Calculate the relative eigenenergies Em/E01 for this ng.
        
        enerlist.append(EmE01) # enerlist now is a list of numpy arrays.
    
    enerlist = np.array(enerlist) # enerlist is now a 2d numpy array.
                                  # Its rows correspond to variation with ng.
                                  # Its columns correspond to the Em/E01 values, 
                                  # where m indicates the energy level.
    
    EmE01_list = [] # The rows will correspond to each energy level. 
                    # The columns will contain the different Em/E01 in function of ng. 
    for i in range(kval):
        EmE01_list.append(enerlist[:,i])
    
    #Create the figure and its axis.
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    # plt.rc('text', usetex=True)

    for i in range(kval):
        ax.plot(ngs, EmE01_list[i], lw=3, label="State $|{}\\rangle$".format(i))
    
    ymin = np.min(EmE01_list)-1
    ymax = np.max(EmE01_list)+1
    
    for ng_i in np.arange(-ng_max+1/2, ng_max, 1/2):
        ax.plot([ng_i, ng_i], [ymin, ymax], c="k", ls="--", lw=2, zorder=0, alpha=0.7)
    
    ax.set_xlabel('$n_g$', fontsize=16)
    ax.set_ylabel('$E_m/E_{01}$', fontsize=16)
    
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([-ng_max, ng_max])
    ax.set_xticklabels(["-1", "-3/4", "-1/2", "-1/4", "0", "1/4", "1/2", "3/4", "1"])
    ax.tick_params(axis='both', which='major', direction="in", labelsize=11)
    
    ax.set_title("$E_J/E_C = {}$".format(round(Ej/Ec, 3)), fontsize=16)
    ax.legend(loc="upper right", fontsize=12)
    
    plt.tight_layout()
    plt.show()
    
    return enerlist, EmE01_list
    


def plot_Em_ng_subplots(kval, 
                        ng_max = 2, 
                        Ej_list = [1, 5, 10, 50], 
                        Ec_list = [1, 1, 1, 1], 
                        num_ngs = 50, 
                        N = 150, 
                        phimax = 1.5*np.pi):
    
    """
    Plot the relative eigenenergies Em/E01 in function of the 
    effective offset charge ng.
    
    INPUTS:
        * kval is the number of states that are going to be plotted.
        * ng_max and num_ngs define the ng grid.
        * N is the same as in the ham() function.
        * Ej_list are the Josephson energies to be used for the plots.
        * Ec_list are the charge energies to be used for the plots.
        * phi_max defines the limits of the plot in the phase space.
    
    OUTPUT
        * A figure showing the relative eigenenergies Em/E01 in function of the 
          effective offset charge ng.
    """

    ngs = np.linspace(-ng_max, ng_max, num_ngs)
    
    # Create the figure and its 4 axis for the subplots.
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,6), constrained_layout=0)
    
    # Reshape the axis matrix to access it with only a single for loop.
    axs = np.reshape(axs, (1,4))[0]
    
    # List for storing the values of E01 in each subplot.
    E01_list = []
    
    for i in range(4):
        
        axi = axs[i]
        Ej  = Ej_list[i]
        Ec  = Ec_list[i]
        
        # Get the energy difference E01 between the ground state and the 1st 
        # excited state at the sweet spot ng=1/2.
        mat = ham(phimax, N, potfun=josephson, param=Ej, ng=1/2, Ej=Ej, Ec=Ec)
        ener = spy.eigh(mat, eigvals_only=True)  
        E0 = ener[0] * 4*Ec
        E1 = ener[1] * 4*Ec
        E01 = E1 - E0 
        E01_list.append(E01) # Store the value of E01.
        
        enerlist = [] # Rows will correspond to the different ng values,
                      # columns will correspond to the different energy states.
        
        for ng in ngs:
            
            mat = ham(phimax, N, potfun=josephson, param=Ej, ng=ng, Ej=Ej, Ec=Ec)
            ener = spy.eigh(mat, eigvals_only=True) 
            ener = ener * 4*Ec # We modify ener so we get E instead of E/(4Ec)
            
            EmE01 = ener / E01 # Calculate the relative eigenenergies Em/E01 for this ng.
            
            enerlist.append(EmE01) # enerlist now is a list of numpy arrays.
        
        enerlist = np.array(enerlist) # enerlist is now a 2d numpy array.
                                      # Its rows correspond to variation with ng.
                                      # Its columns correspond to the Em/E01 values, 
                                      # where m indicates the energy level.
        
        EmE01_list = [] # The rows will correspond to each energy level. 
                        # The columns will contain the different Em/E01 in function of ng. 
        for j in range(kval):
            EmE01_list.append(enerlist[:,j])
        
        # Plot the results.
        for j in range(kval):
            axi.plot(ngs, EmE01_list[j], lw=2)
    
        # Subplot settings
        ymin = np.min(EmE01_list)-1
        ymax = np.max(EmE01_list)+1
        
        # Set x and y label only to some axis so they can share them for better visualization.
        if (i == 2) or (i == 3): # Bottom left and bottom right subplots.
            axi.set_xlabel('$n_g$', fontsize=15)
        if (i == 0) or (i == 2): # Top left and bottom left subplots.
            axi.set_ylabel('$E_m/E_{01}$', fontsize=15)
        
        axi.set_ylim([ymin,ymax])
        axi.set_xlim([-ng_max, ng_max])
        
        axi.set_title("$E_J/E_C = {}$".format(round(Ej/Ec)), fontsize=16) # round() Ej/Ec to an integer
    
        # Only for ng_max=2. Adjust if ng_max has ohter value.
        axi.set_xticks([-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2])
        axi.set_xticklabels(["$-2$", "$-3/2$", "$-1$", "$-1/2$", "$0$", "$1/2$", "$1$", "$3/2$", "$2$"])
    
    # Common legend for all subplots.
    fig.legend(labels = ["State $|0\\rangle$", "State $|1\\rangle$", "State $|2\\rangle$"], 
               ncol=3, 
               loc='lower center',
               fontsize=15)
    
    # Figure geometry, adjust if necessary.
    fig.subplots_adjust(left=0.085, right=0.985, top=0.945, bottom=0.175, hspace=0.28, wspace=0.16)

    plt.show()
    
    return np.asarray(E01_list)



def getener(mat,kval):
    
    ''' 
    Returns the first kval eigenvalues of mat.
    '''
    
    ener=spy.eigvalsh(mat)
    ener0=[]
    for i in range(kval):
        ener0.append(ener[i])
        
    return ener0
