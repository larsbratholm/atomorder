"""
atomorder/parsers.py

File parsers for coordinate files

"""

import numpy as np

def read_coordinates(filename, format_ = "guess"):
    """
    Parses coordinates file

    Parameters
    ---------
    filename: string
        file name
    format_: string
        file format

    """
    # Guess from extension if file format is not given.
    if format_ == "guess":
        format_ = filename.split(".")[-1]

    if format_ == "xyz":
        return read_coordinates_xyz(filename)
    elif self._format == "pdb":
        return read_coordinates_pdb(filename)
    else:
        quit("Error: Unknown file format %s" % format_)

def read_coordinates_xyz(filename):
    """
    Read and parse xyz coordinate file

    Parameters
    ----------
    filename: string
        Coordinates filename

    Returns
    -------
    atoms: array-like
        N-sized array of element symbols, e.g. ['H','C']
    coordinates: array-like
        Nx3-sized array of euclidian coordinates

    """
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        try:
            n_atoms = int(lines[0])
        except ValueError:
            quit("Error reading XYZ-file. Expected no. of atoms in first line \
                  of filefile. Got %s" % lines[0])
        coordinates = np.zeros((n_atoms, 3), dtype = float)
        for l, line in enumerate(lines[2:]):
            try:
                tokens = line.split()
                atoms.append(tokens[0])
                coordinates[l] = np.asarray(tokens[1:4], dtype=float)
            except IndexError, ValueError:
                quit("Error reading line %d in inputfile %s: \n line %s" % (l+3, filename, line))
    atoms = np.asarray(atoms)
    #coordinates -= coordinates.mean(0)
    return atoms, coordinates

def get_coordinates_pdb(filename):
    """
    get_coordinates_pdb(filename)

    Get coordinates from the first chain in a pdb file

    Parameters
    ----------
    filename: string
        PDB3 filename

    Returns
    -------
    atoms: array-like
        N-sized array of element_symbols, e.g. ['H','C']
    coordinates: array-like
        Nx3-sized array of euclidian coordinates

    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    coordinates = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the element symbol can be inferred from the element symbol given in column 3.
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the element symbol
                try:
                    atom = tokens[2][0]
                    if atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                    else:
                        atom = tokens[2][1]
                        if atom in ["H", "C", "N", "O", "S", "P"]:
                            atoms.append(atom)
                        else:
                            exit("Error parsing element symbol for the following line: \n%s" % line)
                except IndexError:
                        exit("Error parsing element symbol for the following line: \n%s" % line)

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n%s" % line)
                # Try to read the coordinates
                try:
                    coordinates.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        coordinates.append(np.asarray([x,y,z],dtype=float))
                    except IndexError, TypeError:
                        exit("Error parsing input for the following line: \n%s" % line)

    coordinates = np.asarray(V)
    atoms = np.asarray(atoms)
    if coordinates.shape[0] != atoms.size:
        error("Mismatch in number of parsed element symbols (%d) and number of parsed coordinates (%d) from PDB file: %s" \
                % (coordinates.shape[0], atoms.size, filename))
    return atoms, coordinates

def write_xyz(coordinates, elements, filename, title = ""):
    """
    Writes the given elements and coordinates to XYZ formated file

    """
    N = elements.size

    with open(filename, "w") as f:
        f.write(str(N)+"\n"+title+"\n")
        for i in xrange(N):
            line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}\n".format(elements[i], *coordinates[i])
            f.write(line)

def write_mol2(coordinates, atoms, filename, split):
    N = coordinates.shape[0]
    origin = []
    for i, n in enumerate(split):
        origin += [i+1]*n

    lines = []
    for i,atom in enumerate(atoms):
        line = "\t{:d} {:s}\t{:.4f}\t{:.4f}\t{:4f} {:s}\t {:d} LIG{:d}\t0.0\n".format(i+1, atom.element_symbol, coordinates[i,0], coordinates[i,1], coordinates[i,2], atom.sybyl, origin[i], origin[i])
        lines.append(line)
    lines.append("@<TRIPOS>BOND\n")
    count = 0
    for atom1 in atoms:
        for atom2 in atom1.bonded_atoms:
            if atom1.mixture_index >= atom2.mixture_index:
                continue
            count += 1
            # TODO print sybyl bond 
            line = "\t{:d}\t{:d}\t{:d}\t1\n".format(count, atom1.mixture_index+1, atom2.mixture_index+1)
            lines.append(line)
    with open(filename, "w") as f:
        f.write("@<TRIPOS>MOLECULE\n\n%d %d 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n" % (N,count))
        for line in lines:
            f.write(line)
# 23 23 0 0 0
#SMALL
#GASTEIGER
#
#@<TRIPOS>ATOM
#      1 C          -1.7175    1.4267   -0.0723 C.ar    1  LIG1        0.0594
#      2 C          -2.9787    0.8162   -0.1439 C.ar    1  LIG1       -0.0477
#      3 C          -3.1199   -0.5684   -0.1047 C.ar    1  LIG1       -0.0609
#      4 C          -1.9971   -1.3773    0.0127 C.ar    1  LIG1       -0.0582
#      5 C          -0.7362   -0.7967    0.0861 C.ar    1  LIG1       -0.0185
#      6 C          -0.5838    0.5925    0.0336 C.ar    1  LIG1        0.1420
#      7 C          -1.6765    2.9418   -0.1207 C.2     1  LIG1        0.1636
#      8 C          -0.3674    3.7031   -0.0636 C.3     1  LIG1       -0.0013
#      9 O          -2.7031    3.5888   -0.2123 O.2     1  LIG1       -0.2922
#     10 H          -3.8773    1.4220   -0.2316 H       1  LIG1        0.0626
#     11 H          -4.1094   -1.0157   -0.1620 H       1  LIG1        0.0618
#     12 H          -2.1033   -2.4587    0.0498 H       1  LIG1        0.0619
#     13 H           0.1370   -1.4366    0.1875 H       1  LIG1        0.0654
#     14 O           0.6987    1.0989    0.1355 O.3     1  LIG1       -0.4253
#     15 H           0.2699    3.4283   -0.9065 H       1  LIG1        0.0309
#     16 H           0.1580    3.4909    0.8698 H       1  LIG1        0.0309
#     17 H          -0.5477    4.7796   -0.1129 H       1  LIG1        0.0309
#     18 C           1.7296    0.7539   -0.6667 C.2     1  LIG1        0.3092
#     19 O           1.6365    0.0338   -1.6412 O.2     1  LIG1       -0.2507
#     20 C           3.0060    1.4639   -0.2836 C.3     1  LIG1        0.0336
#     21 H           2.8360    2.5412   -0.2266 H       1  LIG1        0.0342
#     22 H           3.7897    1.2702   -1.0194 H       1  LIG1        0.0342
#     23 H           3.3467    1.1170    0.6943 H       1  LIG1        0.0342
#@<TRIPOS>BOND
#     1    19    18    2
#     2    22    20    1
#     3    15     8    1
#     4    18    20    1
#     5    18    14    1
#     6    20    21    1
#     7    20    23    1
#     8    10     2    1
#     9     9     7    2
#    10    11     3    1
#    11     2     3   ar
#    12     2     1   ar
#    13     7     1    1
#    14     7     8    1
#    15    17     8    1
#    16     3     4   ar
#    17     1     6   ar
#    18     8    16    1
#    19     4    12    1
#    20     4     5   ar
#    21     6     5   ar
#    22     6    14    1
#    23     5    13    1

