

# Found from analyzing the CCDC 2016 database.
# loose_lower, lower, loose_upper, upper
bond_length_limits = {("As","As"): (2.20, 2.30, 2.65, 2.80),
                      ("As","Br"): (2.20, 2.30, 3.30, 3.40),
                      ("As","Cl"): (2.10, 2.15, 2.40, 3.30),
                      ("As","C" ): (1.75, 1.80, 2.10, 2.20),
                      ("As","F" ): (1.50, 1.60, 1.80, 2.10),
                      ("As","I" ): (2.30, 2.40, 3.70, 3.80),
                      ("As","N" ): (1.60, 1.70, 2.05, 2.70),
                      ("As","O" ): (1.40, 1.60, 2.00, 2.30),
                      ("As","P" ): (2.10, 2.20, 2.40, 2.60),
                      ("As","S" ): (1.90, 2.00, 2.40, 3.10),
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
bond_length_limits[("H","S")] = (1.2,1.3,1.4,1.5)

# make inverse atom order
for key, value in bond_length_limits.items():
    bond_length_limits[key[::-1]] = value

# number of bonds each atom type commonly form
number_bonds = {"As": [3,4],
                "Br": [1],
                "C" : [2,3,4],
                "Cl": [1],
                "F" : [1],
                "H" : [1],
                "I" : [1],
                "N" : [1,2,3,4],
                "O" : [1,2],
                "P" : [3],
                "S" : [1,2,3,4]
                }

# monovalent atoms
monovalent = ["Br","Cl","F","H","I"]

def get_csd_atom_type(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
    Comments refer to the definitions given from the above site

    Parameters
    ----------
    atom: object
        atom object

    Returns
    -------
    sybyl: string
        sybyl atom type

    """

    # TODO name->element_symbol

    # 2.2 If element_symbol is D then atom_type is H
    if atom.element_symbol == "D":
        sybyl = "H"
    # 2.3 If element_symbol is P then atom_type is P.3
    elif atom.element_symbol == "P":
        sybyl = "P.3"
    # 2.5 If element_symbol is C then
    elif atom.element_symbol == "C":
        sybyl = get_csd_atom_type_C(atom)
    # 2.6 If element_symbol is O then
    elif atom.element_symbol == "O":
        sybyl = get_csd_atom_type_O(atom)
    # 2.7 If element_symbol is N then
    elif atom.element_symbol == "N":
        sybyl = get_csd_atom_type_N(atom)
    # 2.8 If element_symbol is S then
    elif atom.element_symbol == "S":
        sybyl = get_csd_atom_type_S(atom)
    # 2.10 If element_symbol is none of the above then atom_type is element_symbol
    else:
        sybyl = atom.element_symbol

    return sybyl

def get_csd_atom_type_C(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching for carbon
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object
    """

    # 2.5.1 If num_bond .ge. 4 then atom_type is C.3
    if atom.num_bond >= 4: return "C.3"
    # 2.5.2 If num_bond .eq. 1 then calculate bond_distance
    if atom.num_bond == 1:
        # TODO calc bond distance
        bond_distance = None

        # 2.5.2.1 If bond_distance .gt. 1.41A then atom_type is C.3
        if bond_distance > 1.41: return "C.3"
        # 2.5.2.2 If bond_distance .le. 1.22A then atom_type is C.1
        if bond_distance <= 1.22: return "C.1"
        # 2.5.2.3 If bond_distance is none of the above then atom_type is C.2
        return "C.2"
    # 2.5.3 If element_symbol is C and none of the above then calculate average_angle about C
    # TODO calc average angle
    average_angle = None
    # 2.5.3.1 If average_angle .le. 115 deg then atom_type is C.3
    if average_angle <= 115: return "C.3"
    # 2.5.3.2 If average_angle .gt. 160 deg then atom_type is C.1
    if average_angle > 160: return "C.1"
    # 2.5.3.3 If average_angle is none of the above then atom_type is C.2
    return "C.2"

def get_csd_atom_type_O(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object

    """

    # 2.6.1 If num_nonmet .eq. 1 then
    if atom.num_bond == 1:
        # TODO get bonding partner this

        # 2.6.1.1 If bond is to carbon .AND. carbon forms a total of 3
        #         bonds, 2 of which are to an oxygen forming only 1
        #         non-metal bond then atom_type is O.co2
        if this.element_symbol == "C" and this.num_bonds == 3:
            for this_neighbor in this.bonded: # TODO ?
                # TODO ?
                if this_neighbor.element_symbol == "O" and this_neighbor.id != this.id \
                        and this_neighbor.num_bond == 1:
                            return "O.co2"
        # 2.6.1.2 If bond is to phosphorus .AND. phosphorous forms at 
        #         least 2 bonds to an oxygen forming only 1 non-metal 
        #         bond then atom_type is O.co2
        if this.element_symbol == "P":
            for this_neighbor in this.bonded:
                if this_neighbor.element_symbol == "O" and this_neighbor.id != this.id \
                        and this_neighbor.num_bond == 1:
                            return "O.co2"
    # 2.6.3 If num_bond .ge. 2 then atom_type is O.3
    if atom.num_bond >= 2: return "O.3"

    # 2.6.4 If element_symbol is O and none of the above then atom_type is O.2
    return "O.2"

def get_csd_atom_type_N(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object

    """

    # 2.7.1 If num_nonmet .eq. 4 then atom_type is N.4
    if atom.num_bond == 4: return "N.4"

    # 2.7.2 If num_nonmet .eq. 1 then calculate bond_distance
    if atom.num_bond == 1:
        # TODO calc bond distance
        bond_distance = None

        # 2.7.2.1 If bond_distance .le. 1.2A then atom_type is N.1
        if bond_distance >= 1.2: return "N.1"

        # 2.7.2.2 If bond_distance .gt. 1.2A then atom_type is N.3
        return "N.3"

    # 2.7.3 If num_nonmet .eq. 3 .AND. one bond is to C--O or C--S then atom_type is N.am
    if atom.num_bond == 3:
        for this in atom.bonded:
            if this.element_symbol == "C":
                for this_neighbor in this.bonded:
                    # At this stage it's not known if a bond is double or single
                    # so just check that O and S only bonds to the parent carbon
                    if this_neighbor.element_symbol in ("O","S") and \
                            this_neighbor.num_bond == 1:
                                return "N.am"
            # 2.7.4 If num_nonmet .eq. 3 otherwise then calculate sum_of_angles around N
            average_angle = None # TODO

            # 2.7.4.1 If sum_of_angles .ge. 350 deg then atom_type is N.pl3
            if average_angle >= 350: return "N.pl3"
            # 2.7.4.2 If sum_of_angles .lt. 350 deg then atom_type is N.3
            return "N.3"

    # 2.7.5 If element_symbol is N and none of the above then calculate average_angle about N
    # check if average_angles is already calculated
    if average_angles not in locals():
        average_angle = None # TODO

    # 2.7.5.1 If average_angle .gt. 160 deg then atom_type is N.1
    if average_angle > 160: return "N.1"
    # 2.7.5.2 If average_angle .le. 160 deg then atom_type is N.2
    return "N.2"



def get_csd_atom_type_S(atom):
    """
    CSD non-matched (3d) deterministic sybyl atom type matching
    from http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html

    Parameters
    ----------
    atom: object
        atom object

    """

    # 2.8.1 If num_nonmet .eq. 3 .AND. 1 bond is to an oxygen with only one
    #       non-metal bond then atom_type is S.o
    if atom.num_bond == 3:
        for this in atom.bonds:
            if this.element_symbol == "O" and this.num_bond == 1: return "S.o"
    # 2.8.2  If num_nonmet .eq. 4 .AND. 2 bonds are to an oxygen with only 
    #        one non-metal bond then atom_type is S.o2
    if atom.num_bond == 4:
        count = 0
        for this in atom.bonds:
            if this.element_symbol == "O" and this.num_bond == 1: count += 1
            if count == 2: return "S.o2"
    # 2.8.3 If num_bond .ge. 2 then atom_type is S.3
    if atom.num_bond >= 2: return "S.3"
    # 2.8.4 If element_symbol is S and none of the above then atom_type is S.2
    return "S.2"

if __name__ == "__main__":
    pass
