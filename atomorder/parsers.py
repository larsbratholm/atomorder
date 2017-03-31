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

def write_mol2(coordinates, atoms, filename, split, match = None):
    
    N = coordinates.shape[0]
    order = np.arange(N, dtype=int)
    rev_order = np.arange(N, dtype=int)
    if match != None:
        largest_row_indices = match.argmax(1)
        unique_indices, counts = np.unique(largest_row_indices, return_counts=True)
        if (counts == 1).all() == False:
            pass
        else:
            order = largest_row_indices
            rev_order = np.asarray([np.where(order == i)[0][0] for i in range(N)])
            #assert(np.allclose(rev_order, match.argmax(0)))

    origin = []
    for i, n in enumerate(split):
        origin += [i+1]*n

    lines = []
    for i in range(N):
        j = order[i]
        line = "\t{:d} {:s}\t{:.4f}\t{:.4f}\t{:4f} {:s}\t {:d} LIG{:d}\t0.0\n".format(i+1, atoms[j].element_symbol, coordinates[j,0], coordinates[j,1], coordinates[j,2], atoms[j].sybyl, origin[j], origin[j])
        lines.append(line)
    lines.append("@<TRIPOS>BOND\n")
    count = 0
    for i in range(N):
        j = order[i]
        atom1 = atoms[j]
        for atom2 in atom1.bonded_atoms:
            if i >= rev_order[atom2.mixture_index]:
                continue
            count += 1
            # TODO print sybyl bond 
            line = "\t{:d}\t{:d}\t{:d}\t1\n".format(count, i+1, rev_order[atom2.mixture_index]+1)
            lines.append(line)
    with open(filename, "w") as f:
        f.write("@<TRIPOS>MOLECULE\n\n%d %d 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n" % (N,count))
        for line in lines:
            f.write(line)
