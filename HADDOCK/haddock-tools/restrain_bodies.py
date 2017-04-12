#!/usr/bin/env python

"""
Creates distance restraints to lock several chains together. Useful to avoid unnatural
flexibility or movement due to sequence/numbering gaps.
"""

from __future__ import print_function

import logging
import sys
import random
import itertools

# Set random seed to have reproducibility
random.seed(917)

# Functions/Methods

def calc_euclidean(i, j):
    return ( (j[0]-i[0])**2 + (j[1]-i[1])**2 + (j[2]-i[2])**2 )**0.5

def read_structure(pdbf, exclude=None):
    """
    Reads a PDB file and returns a list of parsed atoms
    """
    _atoms = set(['CA', 'P']) # Alpha-Carbon (Prot), Backbone Phosphorous (DNA)
    _altloc = set([' ', 'A'])

    if not exclude:
        exclude = set()
    else:
        exclude = set(exclude)

    res_list = []
    with open(pdbf, 'r') as pdb_handle:
        for line in pdb_handle:
            field = line[0:4]
            if field != 'ATOM':
                continue

            aname = line[12:16].strip()
            chain = line[21] if line[21].strip() else line[72:76].strip() # chain ID or segID
            if chain not in exclude and aname in _atoms and line[16] in _altloc:
                resi = int(line[22:26])
                coords = map(float, (line[30:38], line[38:46], line[46:54]))
                res_list.append( (chain, resi, aname, coords) )

    if not res_list:
        logging.critical('[!] PDB File seems empty or no CA/P atoms found: {0}'.format(pdbf))
        sys.exit(1)

    return res_list

def get_bodies(atom_lst, prot_threshold=4.0, dna_threshold=7.5):
    """
    Determines gaps in an atom list following simple distance based criteria.
    Returns continuous fragments.
    """

    bodies = []

    threshold = {'CA': prot_threshold, 'P': dna_threshold}

    body_start = 0
    for i, atom in enumerate(atom_lst[1:], start=1):
        p_atom = atom_lst[i-1]

        chain, resi, aname, xyz = atom
        p_chain, p_resi, p_aname, p_xyz = p_atom

        if (chain == p_chain) and (aname == p_aname): # Internal Gap
            d_xyz = calc_euclidean(xyz, p_xyz)
            if d_xyz >= threshold[aname]:
                logging.debug('[+++] (Internal) Body: {0}:{1}'.format(body_start, i-1))
                bodies.append( (body_start, i-1) )
                body_start = i # Set new beginning

        elif (chain != p_chain) or (aname != p_name): # Different molecules/types
            logging.debug('[+++] Body: {0}:{1}'.format(body_start, i-1))
            bodies.append( (body_start, i-1) )
            body_start = i # Set new beginning

    if not bodies: # Single continuous molecule
        bodies.append( (0, len(atom_lst)) )
    else:
        logging.debug('[+++] Body: {0}:{1}'.format(body_start, i))
        bodies.append( (body_start, i) ) # Last body

    logging.info('[++] Found {0} bodies'.format(len(bodies)))

    return bodies

def build_restraints(bodies):
    """
    Generates distance restraints to maintain the relative
    orientation of the different bodies during the simulations.

    Generates two unique restraints per pair of bodies.

    Each restraint is created using two random atoms on each body
    and using their exact euclidean distance as target distance.
    """

    def pick_residues(body, max_trials=10):
        # Pick two random residues in each body
        # Make sure they are far apart from each other
        n_trials = 0
        while 1:
            try:
                res_i, res_ii = random.sample(body, 2)
            except ValueError:
                # Likely, sample size is 1
                logging.warning('[!] One-sized body found. This may lead to problems..')
                return (body[0], body[0])

            logging.debug('[+++] Trial {0}: {1} & {2}'.format(n_trials, res_i, res_ii))
            if abs(res_i - res_ii) > 3:
                logging.info('[++] Picked residues {0} & {1}'.format(res_i, res_ii))
                return (res_i, res_ii)
            n_trials += 1
            if n_trials == max_trials:
                msg = '[!] Could not pick two unique distant residues in body after {0} tries'
                logging.info(msg.format(max_trials))
                return (res_i, res_ii)

    restraints = []

    n_bodies = range(len(bodies))
    combinations = itertools.combinations(n_bodies, 2)

    for pair_bodies in combinations:
        body_i, body_j = pair_bodies
        logging.debug('[+++] Restraining body {0} to body {1}'.format(body_i, body_j))

        st_body_i, en_body_i = bodies[body_i]
        st_body_j, en_body_j = bodies[body_j]
        res_i, res_ii = pick_residues(range(st_body_i, en_body_i+1))
        res_j, res_jj = pick_residues(range(st_body_j, en_body_j+1))

        logging.info('[++] Created restraint: {0}:{1} <--> {2}:{3}'.format(body_i, res_i, body_j, res_j))
        restraints.append( (res_i, res_j) )
        logging.info('[++] Created restraint: {0}:{1} <--> {2}:{3}'.format(body_i, res_ii, body_j, res_jj))
        restraints.append( (res_ii, res_jj) )

    return restraints

def generate_tbl(atom_lst, restraints):
    """
    Makes a list of TBL-formatted restraints.
    """

    for r in restraints:
        i, j = r
        atom_i, atom_j = atom_lst[i], atom_lst[j]
        dist_ij = calc_euclidean(atom_i[3], atom_j[3])

        tbl = "assign (segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_i)
        tbl+= "(segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_j)
        tbl+= "{0:3.3f} 0.0 0.0".format(dist_ij)
        print(tbl)

def generate_pml(atom_lst, restraints):
    """
    Makes a list of restraints in Pymol format
    """

    for n, r in enumerate(restraints, start=1):
        i, j = r
        atom_i, atom_j = atom_lst[i], atom_lst[j]
        dist_ij = calc_euclidean(atom_i[3], atom_j[3])

        pml = "distance restraint_{0}, ".format(n)
        pml+= "(chain {0[0]} and resi {0[1]} and name {0[2]}), ".format(atom_i)
        pml+= "(chain {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_j)
        print(pml)


if __name__ == '__main__':
    from argparse import ArgumentParser
    
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('structures', nargs='+', 
                   help='PDB structures to restraint')
    ap.add_argument('--exclude', '-e', nargs='+',
                   help='Chains to exclude from the calculation')
    ap.add_argument('--verbose', '-v', action='count')
    args = ap.parse_args()
    
    # Set Logger
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose > 1:
        level = logging.DEBUG
    else:
        level = logging.WARNING
    logging.basicConfig(level=level, format='%(message)s')
    logger = logging.getLogger(__name__)

    # Main logic
    atom_lst = []
    for s in args.structures:
        atom_lst += read_structure(s, exclude=args.exclude)

    bodies = get_bodies(atom_lst)
    restraints = build_restraints(bodies)
    generate_tbl(atom_lst, restraints)
