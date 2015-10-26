#!/usr/bin/env python

import methods
import sqlite3

filesdir="/home/vesso/Code/MC_dih/current/files"

ff_nonbonded="/home/vesso/appstack/gromacs-5.1-icc/share/gromacs/top/amber99sb.ff/ffnonbonded.itp"
ff_bonded="/home/vesso/appstack/gromacs-5.1-icc/share/gromacs/top/amber99sb.ff/ffbonded.itp"
topol_file=filesdir+"/"+"topol.top"
pdb_file=filesdir+"/"+"opt.pdb"

output_atoms=filesdir+"/"+"atoms.txt"
output_bonds=filesdir+"/"+"bonds.txt"
output_pdhdr=filesdir+"/"+"pdhdr.txt"
output_pdhdr_params=filesdir+"/"+"pdhdr-params.txt"

connection=sqlite3.connect(':memory:')
cursor=connection.cursor()

methods.get_ffnonbonded(cursor,ff_nonbonded)
methods.get_ffbonded(cursor,ff_bonded)
methods.get_atoms_from_top(cursor,topol_file)
methods.get_bonds_from_top(cursor,topol_file)
methods.get_coords_from_pdb(cursor,pdb_file)
methods.get_prop_dihedrals_from_top(cursor,topol_file)

methods.dump_bonds(cursor,output_bonds)
methods.dump_atoms(cursor,output_atoms)

methods.dump_proper_dihedrals(cursor,output_pdhdr,output_pdhdr_params)

