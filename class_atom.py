import numpy as np

class atom:
    """
    Represents an 'atom' object used in molecular dynamics simulations.

    Parameters
    ----------
    name : str
        Name of the atom.
    mass : float
        Mass of the atom.
    id_group : int
        Identification number of the group to which the atom belongs.

    Attributes
    ----------
    name : str
        Name of the atom.
    mass : float
        Mass of the atom.
    id_group : int
        Identification number of the group to which the atom belongs.
    DOF : int
        Degree of freedom of the atom.
    N : int
        Total count of atom IDs.
    id : ndarray
        Array containing the identification numbers of the atom.
    position_past : ndarray
        Bidimensional array representing the past positions of the atom.
    position : ndarray
        Bidimensional array representing the current positions of the atom.
    velocity : ndarray
        Bidimensional array representing the velocity of the atom.
    force : ndarray
        Bidimensional array representing the force acting on the atom.

    Notes
    -----
    This class sets up an 'atom' object with the specified attributes, including the atom's name, mass, and the identification number of its group.
    Additionally, it initializes other attributes like 'DOF', 'N', 'id', 'position_past', 'position', 'velocity', and 'force' with default values.

    Examples
    --------
    >>> atom_instance = atom("example", 12.0, 1)
    >>> print(atom_instance.name)
    'example'
    >>> print(atom_instance.mass)
    12.0
    >>> print(atom_instance.id_group)
    1
    """


    def __init__(self, name, mass, id_group):
        """
        Initialize an 'atom' object.

        Parameters
        ----------
        name : str
            Name of the atom.
        mass : float
            Mass of the atom.
        id_group : int
            Identification number of the group to which the atom belongs.

        Attributes
        ----------
        name : str
            Name of the atom.
        mass : float
            Mass of the atom.
        id_group : int
            Identification number of the group to which the atom belongs.
        DOF : int
            Degree of freedom of the atom.
        N : int
            Total count of atom IDs.
        id : list
            List of atom IDs.
        position_past : ndarray
            Array representing the past positions of the atom.
        position : ndarray
            Array representing the current positions of the atom.
        velocity : ndarray
            Array representing the velocity of the atom.
        force : ndarray
            Array representing the forces acting on the atom.

        Notes
        -----
        This constructor sets up an 'atom' object with the specified attributes, including the atom's name, mass, movement status, and the identification number of its group. 
        Additionally, it initializes other attributes like 'DOF', 'N', 'id', 'position_past', 'position', 'velocity', and 'force' with default values.

        Examples
        --------
        >>> atom_instance = atom("example", 12.0, True, 1)
        >>> print(atom_instance.name)
        'example'
        >>> print(atom_instance.mass)
        12.0
        >>> print(atom_instance.move)
        True
        >>> print(atom_instance.id_group)
        1
        """
        self.name = name
        self.mass = mass
        self.id_group = id_group
        self.DOF = 0
        self.N = 0
        self.id = []

        self.position_past = np.array([], dtype=float).reshape(0, 3)
        self.position = np.array([], dtype=float).reshape(0, 3)
        self.velocity = np.array([], dtype=float).reshape(0, 3)
        self.force = np.array([], dtype=float).reshape(0, 3)
        

    def add_id(self, i):
        """
        Add the identification number 'i' to the list of atom IDs.

        Parameters
        ----------
        i : int
            The identification number to be added to the list.

        Notes
        -----
        This method appends the provided identification number to the list of atom IDs (`id`) and increments the total count of IDs (`N`).

        Examples
        --------
        >>> atom_instance = atom("example", 12.0, True, 1)
        >>> atom_instance.add_id(42)
        >>> print(atom_instance.id)
        array([42])
        >>> print(atom_instance.N)
        1
        """
        self.id = np.append(self.id, i)
        self.N += 1


    def add_position_past(self, x, y, z):
        """
        Add the coordinates of the considered atom from the starting configuration to the list position_past.

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
        This method appends the provided coordinates (x, y, z) to the list 'position_past' representing the past positions of the atom.

        Examples
        --------
        >>> atom_instance = atom("example", 12.0, True, 1)
        >>> atom_instance.add_position_past(1.0, 2.0, 3.0)
        >>> print(atom_instance.position_past)
        array([[1.0, 2.0, 3.0]])
        """
        self.position_past = np.vstack([self.position_past, np.array([x, y, z], dtype=float)])


    def add_position(self, x, y, z):
        """
        Add the coordinates of the considered atom to the list position.

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
        This method appends the provided coordinates (x, y, z) to the list 'position' representing the current positions of the atom.

        Examples
        --------
        >>> atom_instance = atom("example", 12.0, True, 1)
        >>> atom_instance.add_position(1.0, 2.0, 3.0)
        >>> print(atom_instance.position)
        array([[1.0, 2.0, 3.0]])
        """
        self.position = np.vstack([self.position, np.array([x, y, z], dtype=float)])

        
    def add_force(self, x, y, z):
        """
        Add the coordinates of the force acting on the considered atom to the list force.

        Parameters
        ----------
        x : float
            X-component of the force.
        y : float
            Y-component of the force.
        z : float
            Z-component of the force.

        Notes
        -----
        This method appends the provided force components (x, y, z) to the list 'force' representing the forces acting on the atom.

        Examples
        --------
        >>> atom_instance = atom("example", 12.0, True, 1)
        >>> atom_instance.add_force(1.0, 2.0, 3.0)
        >>> print(atom_instance.force)
        array([[1.0, 2.0, 3.0]])
        """
        self.force = np.vstack([self.force, np.array([x, y, z], dtype=float)])


    def generate_velocity(self, dt):
        """
        Generate the bidimensional array representing the velocity of each atom.

        Parameters
        ----------
        dt : float
            Time interval.

        Notes
        -----
        Velocities are calculated as the difference between the current position and the position from the previous time step, divided by the time interval.

        Examples
        --------
        >>> atom_instance = atom("example", 12.0, True, 1)
        >>> atom_instance.add_position_past(1.0, 2.0, 3.0)
        >>> atom_instance.add_position(2.0, 4.0, 6.0)
        >>> atom_instance.generate_velocity(0.5)
        >>> print(atom_instance.velocity)
        array([[2.0, 4.0, 6.0]])
        """
        self.velocity = (np.array(self.position, dtype=float) - np.array(self.position_past, dtype=float)) / dt

