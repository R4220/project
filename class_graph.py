import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

parameters = {'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 14}
plt.rcParams.update(parameters)
colors = ['#425840', '#9ACD32', '#b32323', '#191971', '#006400', '#b0edef', '#470303'] # 0DarkGrey, 1YellowGreen, 2FireBrick, 3MidnightBlue, 4Darkgreen, 5PaleTurquoise, 6Darkred


class graph: # da modificare con i nuovi grafici
    """
    Radial Distribution Function (RDF) Calculator.

    Parameters:
    -----------
    filename : str
        The filename for saving the plot.
    Rmax : float
        Maximum distance for RDF calculation.
    atoms : list of str
        Atom types for which RDF is calculated.
    N : int
        Number of bins for histogram.

    Attributes:
    -----------
    filename : str
        The filename for saving the plot.
    Rmax : float
        Maximum distance for RDF calculation.
    type : list of str
        Atom types for which RDF is calculated.
    N : int
        Number of bins for histogram.
    count : ndarray
        Histogram counts for RDF.
    R : ndarray
        Radial distance array.
    dR : float
        Bin size for histogram.
    norm : ndarray
        Normalization array for RDF calculation.
    condition : bool
        True if RDF is calculated for the same type of atoms, False otherwise.
    at1 : ndarray
        Array to store positions of atoms of type 1.
    N1 : int
        Number of atoms of type 1.
    at2 : ndarray
        Array to store positions of atoms of type 2.
    N2 : int
        Number of atoms of type 2.
    matrix : ndarray or None
        Transformation matrix for position conversion.

    Methods:
    --------
    add_position1(x, y, z):
        Add position of an atom of type 1.
    add_position2(x, y, z):
        Add position of an atom of type 2.
    RDF():
        Calculate the Radial Distribution Function.
    normalization(iteration_obj):
        Normalize RDF counts based on the system volume and atom counts.
    plot():
        Plot the Radial Distribution Function.

    Examples:
    ---------
    rdf_calculator = RDF("example", 10.0, ["A", "A"], 100)
    rdf_calculator.add_position1(0.0, 0.0, 0.0)
    rdf_calculator.add_position2(1.0, 1.0, 1.0)
    rdf_calculator.RDF()
    rdf_calculator.normalization(iteration_obj)
    rdf_calculator.plot()
    """


    def __init__(self, filename, Rmax, atoms, N_bin): # da modificare con i nuovi grafici
        """
        Initialize an instance of the 'Class_Name' class.

        Parameters
        ----------
        filename : str
            Name of the file.
        Rmax : float
            Maximum value for R.
        atoms : list
            List of atoms.
        N : int
            Number of elements.

        Attributes
        ----------
        filename : str
            Name of the file.
        Rmax : float
            Maximum value for R.
        type : list
            List of atoms.
        N : int
            Number of bins.
        count : numpy.ndarray
            Array of zeros with size N.
        R : numpy.ndarray
            Array containing values from 0 to Rmax with N elements.
        dR : float
            Value of the second element in R.
        norm : numpy.ndarray
            Normalized array based on R values.
        condition : bool
            Condition based on the equality of the two elements in atoms.
        at1 : numpy.ndarray
            Bidimensional array for atoms type 1.
        N1 : int
            Number of elements in at1.
        at2 : numpy.ndarray
            Bidimensional array for atoms type 2.
        N2 : int
            Number of elements in at2.
        matrix : NoneType
            Placeholder for a matrix attribute.

        Notes
        -----
        This class is designed to represent a certain type of object, providing various attributes for configuration and calculations.

        Examples
        --------
        >>> instance = ClassName("example.txt", 10.0, ["Oxygen", "Oxygen"], 100)
        >>> print(instance.filename)
        'example.txt'
        >>> print(instance.Rmax)
        10.0
        >>> print(instance.type)
        ['Oxygen', 'Oxygen']
        >>> print(instance.N)
        100
        """
        self.filename = filename

        self.Rmax = Rmax
        self.type = atoms
        self.N_bin = N_bin
        self.count = np.zeros(N_bin)
        self.R = np.linspace(0, Rmax, N_bin)
        self.dR = self.R[1]
        self.norm = np.multiply([((i + self.dR)**3 - i**3) for i in self.R[:N_bin]], np.pi * 4 /3)
        self.condition = (atoms[0] == atoms[1])
        self.at1 = np.array([], dtype=float).reshape(0, 3)
        self.N1 = 0
        self.at2 = np.array([], dtype=float).reshape(0, 3)
        self.N2 = 0

        self.Ek = []
        self.Up = []
        self.T = []
        self.F = np.array([], dtype=float).reshape(0, 3)

        self.time = []
   

    def add_position1(self, x, y, z):
        """
        Add the coordinates of the considered atom to the list 'at1'.

        Parameters
        ----------
        x : float
            X-coordinate of the atom.
        y : float
            Y-coordinate of the atom.
        z : float
            Z-coordinate of the atom.

        Notes
        -----
        This method appends the provided coordinates (x, y, z) to the list 'at1' representing the positions of the atom.

        Examples
        --------
        >>> instance = ClassName("example.txt", 10.0, ["Oxygen", "Oxygen"], 100)
        >>> instance.add_position1(1.0, 2.0, 3.0)
        >>> print(instance.at1)
        array([[1.0, 2.0, 3.0]])
        """
        self.at1 = np.vstack([self.at1, np.array([x, y, z], dtype=float)])

    
    def add_position2(self, x, y, z):
        """
        Add the coordinates of the considered atom to the list 'at2'.

        Parameters
        ----------
        x : float
            X-coordinate of the atom.
        y : float
            Y-coordinate of the atom.
        z : float
            Z-coordinate of the atom.

        Notes
        -----
        This method appends the provided coordinates (x, y, z) to the list 'at2' representing the positions of the atom.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> group_instance.atoms[0].add_position2(1.0, 2.0, 3.0)
        >>> print(group_instance.atoms[0].at2)
        array([[1.0, 2.0, 3.0]])
        """
        self.at2 = np.vstack([self.at2, np.array([x, y, z], dtype=float)])


    def RDF(self):
        """
        Calculate the radial distribution function (RDF) for the group.

        Notes
        -----
        This method calculates the distances between the two defined atomic species and adds them to the counting list 'count'.
        The periodic conditions are used to account for the minimum image criterion. In order to do so, reducted coordinates are also used.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> group_instance.Add_atom("atom_2", 15.0)
        >>> group_instance.atoms[0].add_position1(1.0, 2.0, 3.0)
        >>> group_instance.atoms[1].add_position2(2.0, 3.0, 4.0)
        >>> rdf_values = group_instance.RDF()
        >>> print(rdf_values)
        array([0. , 0.5, 1. , 0. , 0. ])
        """
        dist = []

        # Equal atomic species
        if self.condition:
            n = len(self.at1)
            rpos = self.at1
            
            # Generate the reducted coordinates
            for i, _pos in enumerate(self.at1):
                    rpos[i] = np.dot(np.linalg.inv(self.matrix), _pos)

            # Calculate the distances
            for k in range(n - 1):
                rdiff = rpos[k+1:] - rpos[k]
                int_pos = np.round(rdiff)
                rdiff = rdiff - int_pos
                diff = rdiff
                for i, _rpos in enumerate(rdiff):
                        diff[i] = np.dot(self.matrix, _rpos)
                r = np.linalg.norm(diff)
                dist.extend(r[r < self.Rmax])

        # Different atomic species
        else:
            rpos1 = self.at1
            rpos2 = self.at2
            
            # Generate the reducted coordinates
            for i, _pos in enumerate(self.at1):
                rpos1[i] = np.dot(np.linalg.inv(self.matrix), _pos)
            for i, _pos in enumerate(self.at2):
                rpos2[i] = np.dot(np.linalg.inv(self.matrix), _pos)

            # Calculate the distances
            for i in rpos1:
                rdiff = i - rpos2
                rdiff = rdiff - np.round(rdiff)
                diff = rdiff
                for i, _rpos in enumerate(rdiff):
                    diff[i] = np.dot(self.matrix, _rpos)
                r = np.linalg.norm(diff)
                dist.extend(r[r < self.Rmax])

        self.N1 = len(self.at1)
        self.N2 = len(self.at2)

        self.count += np.histogram(dist, bins=self.N_bin, range=(0, self.Rmax))[0]

        # Reset the arrays for the next time step
        self.at1 = np.array([], dtype=float).reshape(0, 3)
        self.at2 = np.array([], dtype=float).reshape(0, 3)


    def normalization(self, iteration_obj):
        """
        Normalize the radial distribution function (RDF).

        Parameters
        ----------
        iteration_obj : object
            Object representing the iteration and containing necessary information.

        Notes
        -----
        This method calculates the volume 'V' and uses it to normalize the radial distribution function 'count' by dividing it by the appropriate density ('rho_1').
        
        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> group_instance.Add_atom("atom_2", 15.0)
        >>> iteration_obj = Iteration()  # Assuming there is a class named 'Iteration'
        >>> group_instance.normalization(iteration_obj)
        """
        V = np.dot(iteration_obj.az, np.cross(iteration_obj.ax, iteration_obj.ay)) ** 2

        if self.condition:
            rho_1 = V / (self.N1 * (self.N1 - 1) * iteration_obj.N_iteration)
        else:
            rho_1 = V / (self.N1 * self.N2 * iteration_obj.N_iteration)

        self.norm = rho_1 * self.norm
        self.count = np.divide(self.count, self.norm)


    def extracting_values(self, iteration_obj):
        """
        Extract and store relevant values from the iteration object.

        Parameters
        ----------
        iteration_obj : iteration
            An instance of the `iteration` class containing simulation data.

        Notes
        -----
        This method extracts and stores forces (`Ftot`), kinetic energy (`Ek`), potential energy (`Up`),
        total force (`F`), and simulation time (`time`) from the provided `iteration_obj` for each group in the system.

        Examples
        --------
        >>> group_instance = group("example_type", 1)
        >>> iteration_instance = iteration([group_instance])
        >>> group_instance.Ek = np.array([10.0])  # Initial kinetic energy
        >>> group_instance.Ftot = np.array([[1.0, 2.0, 3.0]])  # Initial total force
        >>> iteration_instance.U_pot = 50.0  # Initial potential energy
        >>> iteration_instance.N_iteration = 1  # Initial number of iterations
        >>> iteration_instance.dt = 0.001  # Initial time interval
        >>> iteration_instance.extracting_values(iteration_instance)
        >>> print(group_instance.Ek)
        [10.0 20.0]  # Updated kinetic energy
        >>> print(group_instance.Up)
        [50.0]  # Updated potential energy
        >>> print(group_instance.F)
        [[1.0 2.0 3.0]
         [0.0 0.0 0.0]]  # Updated total force
        >>> print(group_instance.time)
        [0.0 0.001]  # Updated simulation time
        """
        F = np.array([], dtype=float).reshape(0, 3)
        self.Ek = np.append(self.Ek, 0)
        for gr in iteration_obj.groups:
            F = np.append(F, gr.Ftot.reshape(1,3), axis=0)
            #self.T[-1] += gr.T
            self.Ek[-1] += gr.Ek
                        
        self.Up = np.append(self.Up, iteration_obj.U_pot)
        self.F = np.append(self.F, np.sum(F, axis=0).reshape(1, 3), axis=0)
        #self.T[-1] = self.T[-1]/len(iteration_obj.groups)
        self.time = np.append(self.time, iteration_obj.dt * iteration_obj.N_iteration)
        

    def plot_RDF(self):
        """
        Plot the radial distribution function.

        Generates a plot of the radial distribution function (RDF) and saves it as an image file named 'RDF_{filename}.png'.

        Notes
        -----
        This method uses Matplotlib to create a plot of the RDF based on the calculated radial distances and counts. The plot
        is saved as an image file in PNG format.

        Examples
        --------
        To generate and save the RDF plot:

        >>> instance = graph("example", 10.0, ["Oxygen", "Oxygen"], 100)
        >>> instance.RDF()  # Assuming RDF values have been calculated
        >>> instance.plot_RDF()

        The RDF plot will be saved as 'RDF_example.png'.
        """
        fig = plt.figure(figsize=(10, 6.18033988769))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.R, self.count, color = colors[3], label='Power spectrum')
        ax.set_xlabel('r ($A$)')
        ax.set_ylabel('g(r)')
        ax.grid()
        plt.savefig(f'RDF_{self.filename}.png')


    def plot_energy(self):
        """
        Plot the graph of kinetic energy and potential energy in function of time.

        Generates a plot of the kinetic energy and potential energy as a function of time and saves it as an image file named 'E_{filename}.png'.

        Notes
        -----
        This method uses Matplotlib to create a plot of the kinetic energy and potential energy over time. The plot is saved as an image file in PNG format.

        Examples
        --------
        To generate and save the energy plot:

        >>> instance = graph("example", 10.0, ["Oxygen", "Oxygen"], 100)
        >>> instance.extracting_values(iteration_obj)  # Assuming values have been extracted
        >>> instance.plot_energy()

        The energy plot will be saved as 'E_example.png'.
        """
        fig = plt.figure(figsize=(10, 6.18033988769))
        axE = fig.add_subplot(1, 1, 1)
        axU = axE.twinx()
        
        axU.plot(self.time, self.Up * 0.001, color = colors[3], label='Potential energy')
        axE.plot(self.time, self.Ek * 0.001, color = colors[6], linewidth = 2,  label='Kinetic energy')
        axE.plot(self.time, self.Ek * 0.001, color = colors[6], linewidth = 2,  label='Kinetic energy')
        axE.set_xlabel('T (ps)')
        axE.set_ylabel('E$_k$ (keV)')
        axE.legend()
        axE.set_xticklabels([str(int(tick)) for tick in axE.get_xticks()])
        axE.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        axE.set_yticklabels([str(float(tick)) for tick in axE.get_yticks()])
        axE.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        axU.set_ylabel('U (keV)')
        axU.legend()
        axU.set_yticklabels([str(float(tick)) for tick in axU.get_yticks()])
        axU.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
        #axE.grid(color=colors[0], linestyle='-')
        #axU.grid(color=colors[0], linestyle='--')

        plt.savefig(f'E_{self.filename}.png')


    def plot_forces(self, iteration_obj):
        """
        Plot the forces acting on the system in function of time.

        Notes
        -----
        This method generates a plot of the forces (components along x, y, and z acting on the system over time
        and saves it as an image file named 'F_{filename}.png'.

        Examples
        --------
        instance = ClassName("example.txt", 10.0, ["Oxygen", "Oxygen"], 100)
        instance.add_position1(1.0, 2.0, 3.0)
        instance.add_position2(2.0, 3.0, 4.0)
        instance.RDF()
        instance.normalization(iteration_obj)  # Assuming there is an 'iteration_obj'
        instance.plot_forces()
        """
        fig = plt.figure(figsize=(10, 3*6.18033988769 +2))
        
        #ax = fig.add_subplot(1, 1, 1)
        Fx = []
        Fy = []
        Fz = []
        for f in self.F:
                Fx = np.append(Fx, f[0])
                Fy = np.append(Fy, f[1])
                Fz = np.append(Fz, f[2])

        axX = fig.add_subplot(3, 1, 1)
        axX.plot(self.time, Fx * 0.001, color = colors[0], label='F$^{tot}_x$')
        
        axZ = fig.add_subplot(3, 1, 3)
        axZ.plot(self.time, Fz * 0.001, color = colors[0],  label='F$^{tot}_z$')

        axY = fig.add_subplot(3, 1, 2)
        axY.plot(self.time, Fy * 0.001, color = colors[0],  label='F$^{tot}_y$')

        for i, gr in enumerate(iteration_obj.groups):
            Fx = []
            Fy = []
            Fz = []
            for f in gr.Ftot_store:
                Fx = np.append(Fx, f[0])
                Fy = np.append(Fy, f[1])
                Fz = np.append(Fz, f[2])
            axX.plot(self.time, Fx * 0.001, color = colors[i+1],  label=f'{gr.id_group}$_x$')
            axY.plot(self.time, Fy * 0.001, color = colors[i+1],  label=f'{gr.id_group}$_y$')
            axZ.plot(self.time, Fz * 0.001, color = colors[i+1],  label=f'{gr.id_group}$_z$')

        axX.set_xlabel('T (ps)')
        axX.set_ylabel('F$_x$ (nN)')
        axX.grid()
        axX.legend()
        axX.set_xticklabels([str(int(tick)) for tick in axX.get_xticks()])
        axX.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        axX.set_yticklabels([str(float(tick)) for tick in axX.get_yticks()])
        axX.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        axY.set_xlabel('T (ps)')
        axY.set_ylabel('F$_y$ (nN)')
        axY.grid()
        axY.legend()
        axY.set_xticklabels([str(int(tick)) for tick in axY.get_xticks()])
        axY.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        axY.set_yticklabels([str(float(tick)) for tick in axY.get_yticks()])
        axY.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        axZ.set_xlabel('T (ps)')
        axZ.set_ylabel('F$_z$ (nN)')
        axZ.grid()
        axZ.legend()
        axZ.set_xticklabels([str(int(tick)) for tick in axZ.get_xticks()])
        axZ.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        axZ.set_yticklabels([str(float(tick)) for tick in axZ.get_yticks()])
        axZ.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        
        plt.savefig(f'F_{self.filename}.png')


    def plot_temperature(self, iteration_obj):
        """
        Plot the temperature of the system in function of time.

        Notes
        -----
        This method generates a plot of the temperature acting on the system over time
        and saves it as an image file named 'T_{filename}.png'.

        Examples
        --------
        instance = ClassName("example.txt", 10.0, ["Oxygen", "Oxygen"], 100)
        instance.add_position1(1.0, 2.0, 3.0)
        instance.add_position2(2.0, 3.0, 4.0)
        instance.RDF()
        instance.normalization(iteration_obj)  # Assuming there is an 'iteration_obj'
        instance.plot_temperature()
        """
        fig = plt.figure(figsize=(10, 6.18033988769))
        
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.T, color = colors[0], label='T$^{tot}$')
        for i, gr in enumerate(iteration_obj.groups):
            ax.plot(self.time, gr.T, color = colors[i+1],  label=f'{gr.id_group}')

        ax.set_xlabel('T (ps)')
        ax.set_ylabel('T (K)')
        ax.grid()
        ax.legend()
        ax.set_xticklabels([str(int(tick)) for tick in ax.get_xticks()])
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.set_yticklabels([str(float(tick)) for tick in ax.get_yticks()])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        
        plt.savefig(f'T_{self.filename}.png')

        