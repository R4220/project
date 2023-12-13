import numpy as np
from class_atom import atom

class group:
    """
    This class represents a group of atoms used to divide the whole atoms in a molecular dynamics simulation.

    Attributes
    ----------
    type : str
        Atomic type in the group.
    id_group : int
        Identification number of the group.
    atoms : list
        List of atom instances in the group.
    id_tot : numpy.ndarray
        Array containing the identification numbers of the atoms inside the group.
    DOF : float
        Degree of freedom of the group.
    Ek : float
        Kinetic energy of the group.
    Ftot : numpy.ndarray
        Total force acting on the group.
    force : numpy.ndarray
        Bidimensional array representing the forces against the atoms in the group.
    Vtot : numpy.ndarray
        Bidimensional array representing the total velocity of the atoms in the group.
    velocity : numpy.ndarray
        Bidimensional array representing the velocity of the atoms in the group.

    Methods
    -------
    __init__(type, id_group)
        Initialize a Group instance.
    Add_atom(name, mass)
        Add an atom to the group.
    Add_force(x, y, z)
        Add the coordinates of the force acting on the considered atom from the current time step to the force list.
    Kinetic_energy(dt)
        Calculate the kinetic energy of the group.
    Generate(dt)
        Generate output data for the group.

    Notes
    -----
    This class provides functionality to simulate molecular dynamics at the group level, including methods to add atoms, calculate kinetic energy, and generate output data.

    Examples
    --------
    >>> group_instance = group("C O H", 0)
    >>> print(group_instance.type)
    'C O H'
    >>> print(group_instance.id_group)
    0
    """

    def __init__(self, type, id_group):
        """
        Initialize a Group instance.

        Parameters
        ----------
        type : str
            Atomic type in the group.
        id_group : int
            Identification number of the group.

        Attributes
        ----------
        atoms : list
            List of atom instances in the group.
        type : str
            Atomic type in the group.
        id_group : int
            Identification number of the group.
        id_tot : ndarray
            Array containing the identification numbers of the atoms inside the group.
        DOF : float
            Degree of freedom of the group.
        Ek : float
            Kinetic energy of the group.
        Ftot : ndarray
            Total force acting on the group.
        force : ndarray
            Bidimensional array representing the forces against the atoms in the group.
        Vtot : ndarray
            Bidimensional array representing the total velocity of the atoms in the group.
        velocity : ndarray
            Bidimensional array representing the velocity of the atoms in the group.

        Notes
        -----
        This constructor sets up a 'Group' object with the specified attributes, including the group's type and the identification number of the group.
        Additionally, it initializes other attributes like 'atoms', 'id_tot', 'DOF', 'Ek', 'Ftot', 'force', 'Vtot', and 'velocity' with default values.

        Examples
        --------
        >>> group_instance = group("C O H", 0)
        >>> print(group_instance.type)
        'C O H'
        >>> print(group_instance.id_group)
        0
        """
        self.atoms = []
        self.type = type
        self.id_group = id_group
        self.id_tot = np.array([], dtype=int)
        self.DOF = 0.0
        self.Ek = 0.0

        self.Ftot = np.array([], dtype=float).reshape(0, 3)
        self.force = np.array([], dtype=float).reshape(0, 3)

        self.Vtot = np.array([], dtype=float).reshape(0, 3)
        self.velocity = np.array([], dtype=float).reshape(0, 3)


    def Add_atom(self, name, mass):
        """
        Add an atom to the group.

        Parameters
        ----------
        name : str
            Name of the atom.
        mass : float
            Mass of the atom.

        Notes
        -----
        This method adds an 'atom' instance to the 'atoms' list of the group, creating the atom with the specified name, mass, movement status, and identification number of the group.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> print(len(group_instance.atoms))
        1
        >>> print(group_instance.atoms[0].name)
        'atom_1'
        >>> print(group_instance.atoms[0].mass)
        12.0
        """
        self.atoms = np.append(self.atoms, atom(name, mass, self.id_group))

    
    def Add_force(self, x, y, z):
        """
        Add the coordinates of the force acting on the considered atom from the current time step to the force list.

        Parameters
        ----------
        x : float
            X-coordinate of the force.
        y : float
            Y-coordinate of the force.
        z : float
            Z-coordinate of the force.

        Notes
        -----
        This method appends the specified force vector (x, y, z) to the 'force' list of the group, representing the forces acting on the atoms in the group during the current time step.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_force(1.0, 2.0, -0.5)
        >>> print(group_instance.force)
        array([[ 1.0,  2.0, -0.5]])
        >>> group_instance.Add_force(0.5, -1.0, 2.5)
        >>> print(group_instance.force)
        array([[ 1.0,  2.0, -0.5],
                   [ 0.5, -1.0,  2.5]])
        """
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])


    def Kinetic_energy(self, dt):
        """
        Calculate the kinetic energy of the group.

        Parameters
        ----------
        dt : float
            Time interval.

        Notes
        -----
        This method calculates the kinetic energy of the group based on the velocities of its atoms. It iterates over each atom in the group, generates the velocity array for each atom type, and computes the total velocity ('Vtot') of the atoms in the group. Then, it removes the mean velocity and calculates the thermal kinetic energy ('Ek') using the formula:

        .. math::

            E_k = \\frac{1}{2} m_i \\sum_{i} \\left\\| \\mathbf{v}_i - \\mathbf{V}_{\\text{tot}} \\right\\|^2

        where:
        - \( m_i \) is the mass of the atom.
        - \( \\mathbf{v}_i \) is the velocity vector of the atom.
        - \( \\mathbf{V}_{\\text{tot}} \) is the total velocity of the atoms in the group.

        The calculated kinetic energy is stored in the 'Ek' attribute.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> group_instance.Add_atom("atom_2", 15.0)
        >>> group_instance.atoms[0].velocity = np.array([[1.0, 2.0, -0.5]])
        >>> group_instance.atoms[1].velocity = np.array([[0.5, -1.0, 2.5]])
        >>> group_instance.Kinetic_energy(0.1)
        >>> print(group_instance.Ek)
        1.8251677817711523e-05
        """
        # Generate the velocity array for each atom type
        for at in self.atoms:
            at.generate_velocity(dt)
            self.velocity = np.append(self.velocity, at.velocity)

        # Calculate the mean velocity
        self.Vtot = np.sum(self.velocity, axis=0) / len(self.id_tot)

        # Removing the mean velocity to calculate the thermal energy
        for at in self.atoms:
            #self.Ek += 0.5 * float(at.mass) * np.sum(np.linalg.norm(at.velocity, axis=1) ** 2) * 0.0001036426948415943
            self.Ek += 0.5 * float(at.mass) * np.sum(np.linalg.norm(at.velocity - self.Vtot, axis=1) ** 2) * 0.0001036426948415943
        # Reset the arrays for the next time step
        self.velocity = np.array([], dtype=float).reshape(0, 3)
        
        
    def Generate(self, dt):
        """
        Generate output data for the group.

        Parameters
        ----------
        dt : float
            Time interval.

        Returns
        -------
        numpy.ndarray
            Array containing output data for each atom in the group.

        Notes
        -----
        This method calculates the kinetic energy, temperature, and total force of the group. It generates an output line for each atom in the group, including information such as the atom's name, position, velocity, force, and group identification number.
        To calculate the temperature the method use the formula:

        .. math::

            T = \\frac{2 E_k}{DOF * K_B}

        where:
        - \( E_k \) is the kinetic energy.
        - \( DOF \) is th enumber of degree of freedom.
        - \( K_B \) is the Boltzmann constant in eV/K.

        Examples
        --------
        >>> group_instance = group("example_type", 1, True)
        >>> group_instance.Add_atom("atom_1", 12.0)
        >>> group_instance.Add_atom("atom_2", 15.0)
        >>> group_instance.atoms[0].position = np.array([[1.0, 2.0, -0.5]])
        >>> group_instance.atoms[0].velocity = np.array([[1.0, 2.0, -0.5]])
        >>> group_instance.atoms[0].force = np.array([[0.5, -1.0, 2.5]])
        >>> group_instance.Generate(0.1)
        array(['atom_1\\t  1.0\\t  2.0\\t  -0.5\\t  1.0\\t  2.0\\t  -0.5\\t  0.5\\t  -1.0\\t  2.5\\t  1',
               'atom_2\\t  0.0\\t  0.0\\t  0.0\\t  0.0\\t  0.0\\t  0.0\\t  0.0\\t  0.0\\t  0.0\\t  1'],
              dtype='<U78')
        """
        # Calculate the kinetic energy
        self.Kinetic_energy(dt)

        # Calculate temperature
        self.T = (2 * self.Ek) / (self.DOF * 8.617333262145e-5)

        # Calculate total force
        self.Ftot = np.sum(self.force, axis=0)

        # Generate output line
        body = np.array([], dtype=str)
        for at in self.atoms:
            for i in range(at.N):
                    body = np.append(body, f'{at.name}\t  {at.position[i][0]}\t  {at.position[i][1]}\t  {at.position[i][2]}\t  {at.velocity[i][0]}\t  {at.velocity[i][1]}\t  {at.velocity[i][2]}\t  {at.force[i][0]}\t  {at.force[i][1]}\t  {at.force[i][2]}\t  {at.id_group}')

        # Reset the arrays for the next time step
            at.position_past = at.position
            at.position = np.array([], dtype=float).reshape(0, 3)
            at.force = np.array([], dtype=float).reshape(0, 3)
        self.force = np.array([], dtype=float).reshape(0, 3)
        return body
