import numpy as np

class Constants(object):
    """
    Constants()

    Constructor for constants used in the script


    Attributes
    ----------
    bond_length_limits: dict
        Dictionary of loose and tight distance limits on bond lengths
    number_bonds: dict
        Number of bonds the supported atom types can form
    monovalent: list
        Monovalent atom types

    """
    # TODO replace this with a hidden markov model or similar
    def __init__(self):

        # Found from analyzing the CCDC 2016 database.
        # loose_lower, lower, loose_upper, upper
        self.bond_length_limits = {#("As","As"): (2.20, 2.30, 2.65, 2.80),
                                   #("As","Br"): (2.20, 2.30, 3.30, 3.40),
                                   #("As","Cl"): (2.10, 2.15, 2.40, 3.30),
                                   #("As","C" ): (1.75, 1.80, 2.10, 2.20),
                                   #("As","F" ): (1.50, 1.60, 1.80, 2.10),
                                   #("As","I" ): (2.30, 2.40, 3.70, 3.80),
                                   #("As","N" ): (1.60, 1.70, 2.05, 2.70),
                                   #("As","O" ): (1.40, 1.60, 2.00, 2.30),
                                   #("As","P" ): (2.10, 2.20, 2.40, 2.60),
                                   #("As","S" ): (1.90, 2.00, 2.40, 3.10),
                                   ("Br","Br"): (2.20, 2.20, 2.80, 4.00),
                                   ("Br","C" ): (1.75, 1.80, 2.00, 2.10),
                                   ("Br","I" ): (2.50, 2.55, 3.00, 3.50),
                                   ("Br","P" ): (2.00, 2.10, 2.70, 3.00),
                                   ("C" ,"Cl"): (1.60, 1.65, 1.85, 2.00),
                                   ("C" ,"C" ): (1.10, 1.25, 1.70, 1.80),
                                   ("C" ,"F" ): (1.20, 1.25, 1.45, 1.50),
                                   ("C" ,"H" ): (0.85, 0.95, 1.15, 1.25),
                                   ("C" ,"I" ): (1.90, 2.00, 2.20, 2.25),
                                   ("C" ,"N" ): (1.00, 1.10, 1.60, 1.70),
                                   ("C" ,"O" ): (1.00, 1.15, 1.50, 1.60),
                                   ("C" ,"P" ): (1.45, 1.65, 1.95, 2.00),
                                   ("C" ,"S" ): (1.50, 1.60, 1.90, 2.00),
                                   ("Cl","I" ): (2.30, 2.35, 2.75, 3.10),
                                   ("Cl","N" ): (1.60, 1.65, 1.80, 1.90),
                                   ("Cl","O" ): (1.20, 1.30, 1.50, 1.60),
                                   ("Cl","P" ): (1.90, 1.95, 2.20, 2.40),
                                   ("Cl","S" ): (1.90, 1.95, 2.45, 3.10),
                                   ("F" ,"I" ): (1.80, 1.90, 2.15, 3.00),
                                   ("F" ,"P" ): (1.40, 1.50, 1.65, 1.70),
                                   ("F" ,"S" ): (1.40, 1.45, 1.75, 1.80),
                                   ("H" ,"N" ): (0.80, 0.95, 1.10, 1.25),
                                   ("H" ,"O" ): (0.70, 0.85, 1.10, 1.30),
                                   ("I" ,"I" ): (2.60, 2.70, 3.20, 3.60),
                                   ("I" ,"N" ): (1.90, 1.95, 2.55, 2.65),
                                   ("I" ,"O" ): (1.55, 1.60, 2.60, 3.00),
                                   ("I" ,"P" ): (2.30, 2.35, 2.60, 3.00),
                                   ("I" ,"S" ): (2.30, 2.40, 2.95, 3.25),
                                   ("N" ,"N" ): (1.00, 1.10, 1.50, 1.60),
                                   ("N" ,"O" ): (1.10, 1.15, 1.50, 1.60),
                                   ("N" ,"P" ): (1.40, 1.50, 1.80, 2.10),
                                   ("N" ,"S" ): (1.40, 1.50, 1.75, 1.85),
                                   ("O" ,"O" ): (1.20, 1.25, 1.55, 1.60),
                                   ("O" ,"P" ): (1.35, 1.40, 1.80, 1.90),
                                   ("O" ,"S" ): (1.35, 1.40, 1.65, 1.80),
                                   ("P" ,"P" ): (1.95, 2.00, 2.35, 2.60),
                                   ("P" ,"S" ): (1.75, 1.85, 2.20, 2.30),
                                   ("S" ,"S" ): (1.90, 2.00, 2.40, 2.60)
                             }
        # hydrogen bond lengths were taken from neutron diffraction data that
        # didn't have many S-H hits, so guesstimate them
        self.bond_length_limits[("H","S")] = (1.2,1.3,1.4,1.5)

        ## make inverse atom order
        #for key, value in self.bond_length_limits.items():
        #    self.bond_length_limits[key[::-1]] = value

        # Number of bonds each atom type commonly form
        # Doubles as list of atom types implemented
        self.number_bonds = {#"As": np.asarray( [3,4]    , dtype=int) ,
                             "Br": np.asarray( [1]      , dtype=int) ,
                             "C" : np.asarray( [2,3,4]  , dtype=int) ,
                             "Cl": np.asarray( [1]      , dtype=int) ,
                             "F" : np.asarray( [1]      , dtype=int) ,
                             "H" : np.asarray( [1]      , dtype=int) ,
                             "I" : np.asarray( [1]      , dtype=int) ,
                             "N" : np.asarray( [1,2,3,4], dtype=int) ,
                             "O" : np.asarray( [1,2,3]  , dtype=int) ,
                             "P" : np.asarray( [3]      , dtype=int) ,
                             "S" : np.asarray( [1,2,3,4], dtype=int)
                             }
        # monovalent atoms
        self.monovalent = ["Br","Cl","F","H","D","I"]

        # Properties of different sybyl atom types for bonding
        # This is for very internal use. Basically the properties are
        # {num_bonds: (no. pi electrons, participating lone pairs, base charge)}
        # (num_bonds + no_pi - charge) should equal the valency
        # monovalent atoms is {1: (0,0,0)}
        # TODO find a molecule with S.2 so I can understand how it works
        self.sybyl_bonds = {
                            "C.3"  : {4: (0,0,0)},
                            "C.2"  : {3: (1,0,0)},
                            "C.1"  : {2: (2,0,0)},
                            "O.co2": {1: (1,1,0)},
                            "O.3"  : {2: (0,0,0), 3: (0,0,+1)},
                            "O.2"  : {1: (1,1,0)},
                            "N.4"  : {4: (0,0,+1)},
                            "N.3"  : {3: (0,0,0)},
                            "N.2"  : {2: (1,0,0)},
                            "N.am" : {3: (1,1,+1)},
                            "N.pl3": {3: (1,1,+1)},
                            "N.1"  : {1: (2,0,0), 2: (2,0,+1)},
                            "S.o"  : {3: (0,0,+1)},
                            "S.o2" : {4: (2,0,0)},
                            "S.3"  : {2: (0,0,0)},
                            "Br"   : {1: (0,0,0)},
                            "Cl"   : {1: (0,0,0)},
                            "F"    : {1: (0,0,0)},
                            "I"    : {1: (0,0,0)},
                            "H"    : {1: (0,0,0)},
                            }

    def get_bond_length_limits(self, element1, element2):
        if element1 > element2:
            element1, element2 = element2, element1

        if (element1, element2) in self.bond_length_limits:
            return self.bond_length_limits[(element1, element2)]

        # Return limits such that no bonds are formed
        # if the elements are not in self.bond_length_limits
        return (0,0,0,0)

class Settings(object):
    """
    Settings()

    General settings


    Attributes
    ----------
    reactant_filenames: list
        Coordinate files for reactants
    product_filenames: list
        Coordinate files for products
    #TODO

    """
    def __init__(self):
        self.reactant_filenames = [None]
        self.product_filenames = [None]
        self.print_level = 1
        self.file_format = None
        self.method = None
        self.create_atoms = False
        self.rotation_objective = False
        self.bond_objective = False

    def update(self, args):
        """
        Update the object from parsed command line arguments

        """

        self.reactant_filenames = args.reactants
        self.product_filenames = args.products
        self.print_level = args.print_level
        self.file_format = args.format
        self.method = args.method

        # Sets flags needed to define the pipeline of the method
        self.construct_pipeline()

        self.atomic_sybyl_weight = args.weight
        self.annealing_method = "multiplication" # multiplication/addition
        self.initial_inverse_temperature = 1e-3
        self.final_inverse_temperature = 1e3
        self.max_annealing_iterations = 10**4
        self.max_relaxation_iterations = 10**4
        self.max_softassign_iterations = 10**4
        self.annealing_convergence_threshold = 1e-3
        self.relaxation_convergence_threshold = 1e-4
        self.softassign_convergence_threshold = 1e-4

    def construct_pipeline(self):
        """
        Sets flags needed to define the pipeline of the method

        """

        if self.method == "rotate":
            self.rotation_objective = True
            self.atomic_objective = False
            self.bond_objective = False

        if self.method == "no-bond":
            self.rotation_objective = True
            self.atomic_objective = True
            self.bond_objective = False

        elif self.method == "full":
            self.rotation_objective = True
            self.atomic_objective = True
            self.bond_objective = True

        elif self.method == "info":
            self.create_atoms = True
            self.rotation_objective = False
            self.atomic_objective = False
            self.bond_objective = False

        # sanity check override
        if self.bond_objective == True or self.atomic_objective == True:
            self.create_atoms = True

        # TODO more bond


    # TODO add parameters in ordering algorithm
