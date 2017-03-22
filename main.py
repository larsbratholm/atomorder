#!/usr/bin/env python2

from atomorder import parse_args, settings, Reaction, Ordering


if __name__ == "__main__":
    # Create the structure objects of the
    # reactants and products
    reaction = Reaction()
    # Do the atom assignment / ordering
    Ordering(reaction)
