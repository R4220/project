import numpy as np
import matplotlib.pyplot as plt

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


    def __init__(self, filename, Rmax, atoms, N): # da modificare con i nuovi grafici
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
        self.N = N
        self.count = np.zeros(N)
        self.R = np.linspace(0, Rmax, N)
        self.dR = self.R[1]
        self.norm = np.multiply([((i + self.dR)**3 - i**3) for i in self.R[:N]], np.pi * 4 /3)
        self.condition = (atoms[0] == atoms[1])
        self.at1 = np.array([], dtype=float).reshape(0, 3)
        self.N1 = 0
        self.at2 = np.array([], dtype=float).reshape(0, 3)
        self.N2 = 0

        self.Ek = [0]
        self.Up = [0]
        self.T = [0]
        self.F = np.array([0, 0, 0], dtype=float).reshape(1, 3)

        self.time = [0]
   

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

        self.count += np.histogram(dist, bins=self.N, range=(0, self.Rmax))[0]

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


    '''def kinetic_energy(self, iteration_obj):
        for gr in iteration_obj.groups:
            self.Ek[-1] += gr.Ek
        self.Ek = np.append(self.Ek, 0)

    def potential_energy(self, iteration_obj):
        for gr in iteration_obj.groups:
            self.Up[-1] += gr.Ftot
        self.Up = np.append(self.Up, 0)


    def temperature(self, iteration_obj):
        for gr in iteration_obj.groups:
            self.T[-1] += gr.T
        self.T[-1] = self.T[-1]/len(iteration_obj.groups)
        self.Up = np.append(self.Up, 0)


    def force(self, iteration_obj):
        for gr in iteration_obj.groups:
            self.force[-1] += gr.Ftot
        self.force = np.append(self.force, [0, 0, 0])'''


    def extracting_values(self, iteration_obj):
        F = np.array([], dtype=float).reshape(0, 3)
        for gr in iteration_obj.groups:
            F = np.append(F, gr.Ftot.reshape(1,3), axis=0)
            #self.T[-1] += gr.T
            self.Ek[-1] += gr.Ek
        #print('tot', F)
                
        self.Up = np.append(self.Up, iteration_obj.U_pot)
        #print('gr', np.sum(F, axis=0).reshape(1,3))
        #print(self.F)
        self.F = np.append(self.F, np.sum(F, axis=0).reshape(1, 3), axis=0)
        #print(self.F)
        #self.T[-1] = self.T[-1]/len(iteration_obj.groups)
        self.Ek = np.append(self.Ek, 0)
        self.time = np.append(self.time, iteration_obj.dt * iteration_obj.N_iteration)
        

    def plot_RDF(self):
        """
        Plot the radial distribution function.

        Notes
        -----
        This method generates a plot of the radial distribution function (RDF) and saves it as an image file named '{filename}.png'.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> group_instance.Add_atom("atom_2", 15.0)
        >>> group_instance.RDF()
        >>> group_instance.normalization(iteration_obj)  # Assuming there is an 'iteration_obj'
        >>> group_instance.plot()
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.R, self.count, label='Power spectrum')
        ax.set_xlabel('r ($A$)')
        ax.set_ylabel('g(r)')
        ax.grid()
        plt.savefig(f'RDF_{self.filename}.png')
        #plt.show()


    def plot_energy(self):
        #print('Up: ', self.Up, len(self.Up))
        #print('\nF: ', self.F, len(self.F))
        #print('\nEk: ', self.Ek, len(self.Ek))
        #print('\nTime: ', self.time, len(self.time))
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time, self.Ek, label='Kinetic energy')
        ax.plot(self.time, self.Up, label='Potential energy')
        ax.set_xlabel('T (ps)')
        ax.set_ylabel('E (eV)')
        ax.grid()
        ax.legend()
        plt.savefig(f'E_{self.filename}.png')
        #plt.show()

    def plot_forces(self):
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        print('F', self.F)
        Fx = []
        Fy = []
        Fz = []
        for f in self.F:
            Fx = np.append(Fx, f[0])
            Fy = np.append(Fy, f[1])
            Fz = np.append(Fz, f[2])

        print('Fx', Fx)
        ax.plot(self.time, Fx, label='F$_x$')
        ax.plot(self.time, Fy, label='F$_y$')
        ax.plot(self.time, Fz, label='F$_z$')
        ax.set_xlabel('T (ps)')
        ax.set_ylabel('F (pN)')
        ax.grid()
        ax.legend()
        plt.savefig(f'F_{self.filename}.png')
        #plt.show()



        