# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
__doc__="""
This module contains the Molecule class. It allows modifications
of structure files on the residue level:
Basic Usage:
     >>> model = Model().read(args['-f']) # read structure file
     >>> res = model.residues[0]   # pick first residue

     The Molecule instance contains:
     - res.atoms   -> list of atoms
     - res.chain   -> the chain the residue belongs to
     .
     .
     Some methods:
     >>> res.set_resname('CALA')  # change residue name
     >>> a = Atom(name = 'DUM', x=[1,2,3])  # create dummy atom
     >>> res.insert_atom(0,a)    # insert atom in residue at pos 0
     >>> res.append(a)           # append atom at the end
     >>> phi = res.get_phi()     # calculate phi angle
     >>> al = res.fetchm(['CA','C','N']) select multiple atoms (in that order)
     >>> res.remove_atom(al[0]) delete atom
     >>> del res['CA']   delete CA atom
     >>> del res[-1]     delete last atom
     .
     .
"""
import sys
from atomselection import *
import library, copy
from rotamer import _aa_chi

class Molecule(Atomselection):
    """ Storage class for a Molecule/residue"""
    
    def __init__(self, **kwargs):
        Atomselection.__init__(self)
        self.natoms = 0
        self.previous = None
        self.next = None
        self.resname = ''
        self.orig_id = 0
        self.real_resname = ''
        self.chain = None
        self.model = None
        self.chain_id = ''
        self.id = 0
        for key, val in kwargs.items():
            setattr(self,key,val)

    def __str__(self):
        s = '<Molecule: id = %d name = %s chain_id = %s '\
            % (self.id, self.resname, self.chain_id)
        s+=' natoms = %d>' % len(self.atoms)
        return s

    def __getitem__(self, item):
        return self.fetch(item)[0]

    def __delitem__(self,item):
        """ delete atom """
        if item.isdigit():
            atom = self.atoms[item]
            self.remove_atom(atom)
        else:
            # assume its the atom name
            atom = self.fetch(item)[0]
            if atom:
                self.remove_atom(atom)

    def has_atom( self, atom_name ):
        if self.fetch(atom_name):
            return True
        else:
            return False
                
    def new_aa(self, aa, hydrogens = True):
        aa = aa.upper()
        if len(aa) == 1:
            resname = library._aacids_dic[aa]
        else:
            resname = aa
        self.resname = resname
        self.unity = 'A'
        self.atoms = []
        for i, entry in enumerate(library._aacids[resname]):
            if hydrogens == False and entry[0][0] == 'H':
                continue
            else:
                a = Atom()
                a.id = i+1
                a.name = entry[0]
                a.symbol = a.name[0]
                a.x = entry[1]
                a.occ = 1.
                a.resname = resname
                a.m = library._atommass[a.symbol]
                a.unity = 'A'
                self.atoms.append(a)
        self.set_resid(1)
        return self

    def get_real_resname(self):
        dic = {'LYP':'LYS','LYSH':'LYS','LYN':'LYS','LSN':'LYS','CYM':'CYS',
               'CYS2':'CYS','CYN':'CYS','HIE':'HIS','HIP':'HIS','CYSH':'CYS',
               'HID':'HIS','HISA':'HIS','HISB':'HIS','HSE':'HIS','HSP':'HIS',
               'HSD':'HIS','HISH':'HIS','HISD':'HIS','ASH':'ASP','ASPP':'ASP','ASPH':'ASP',
               'GLH':'GLU','GLUH':'GLU','GLUP':'GLU','SEP':'SEP','SEQ':'SEQ',
               }
        if dic.has_key(self.resname): self.real_resname =  dic[self.resname]
        else: self.real_resname = self.resname

    def get_psi(self, degree = False):
        chidx = self.chain.residues.index(self)
        if chidx == len(self.chain.residues)-1:
            return -999.
        else:
            N, CA, C = self.fetchm(['N','CA','C'])
        next_mol = self.chain.residues[chidx+1]
        N2 = next_mol['N']
        dih = N.dihedral(CA,C,N2)
        if not degree:
            return dih
        else:
            return dih*180/pi

    def set_psi(self, degree, propagate = True ):
        psi = self.get_psi()
        diff = degree*pi/180. - psi
        CA, C = self.fetchm(['CA','C'])
        R = Rotation( CA.x, C.x )
        rot_atoms = self.fetch_atoms(['O','OXT','O1','O2','OC1','OC2'])
        if propagate:
            chidx = self.chain.residues.index(self)
            if chidx != len(self.chain.residues)-1:
                for r in self.chain.residues[chidx+1:]:
                    for atom in r.atoms:
                        rot_atoms.append( atom )
        for atom in rot_atoms:
            atom.x = R.apply(atom.x, diff)

    def set_psi_down( self, degree, propagate = True ):
        psi = self.get_psi()
        diff = degree*pi/180. - psi
        CA, C = self.fetchm(['CA','C'])
        R = Rotation( C.x, CA.x )
        rot_atoms = self.fetch_atoms(['O','OXT','O1','O2','OC1','OC2','CA','C'], inv = True)
        if propagate:
            chidx = self.chain.residues.index(self)
            if chidx != 0:
                for r in self.chain.residues[:chidx]:
                    for atom in r.atoms:
                        rot_atoms.append( atom )
        for atom in rot_atoms:
            atom.x = R.apply(atom.x, diff)
        


    def get_phi(self,degree=False):
        """calculate phi angle (C-1)-N-CA-C"""
        chidx = self.chain.residues.index(self)
        if chidx == 0:
            return -999
        else:
            previous = self.chain.residues[chidx-1]
        C = previous.fetchm(['C'])[0]
        N,CA,C2 = self.fetchm(['N','CA','C'])
        dih = C.dihedral(N,CA,C2) 
        if not degree:
            return dih
        else:
            return dih*180/pi

    def set_phi(self, degree, propagate = True ):
        if self.resname == 'PRO': return # does not work
        phi = self.get_phi()
        diff = degree*pi/180. - phi
        N, CA = self.fetchm(['N','CA'])
        R = Rotation( N.x, CA.x )
        rot_atoms = self.fetch_atoms(['N','H','CA'],inv = True)
        if propagate:
            chidx = self.chain.residues.index(self)
            if chidx != len(self.chain.residues)-1:
                for r in self.chain.residues[chidx+1:]:
                    for atom in r.atoms:
                        rot_atoms.append( atom )
        for atom in rot_atoms:
            atom.x = R.apply(atom.x, diff)

    def set_phi_down(self, degree, propagate = True ):
        if self.resname == 'PRO': return # does not work
        phi = self.get_phi()
        diff = degree*pi/180. - phi
        N, CA = self.fetchm(['N','CA'])
        R = Rotation( CA.x, N.x )
        rot_atoms = self.fetch_atoms('H')
        if propagate:
            chidx = self.chain.residues.index(self)
            if chidx != 0:
                for r in self.chain.residues[:chidx]:
                    for atom in r.atoms:
                        rot_atoms.append( atom )
        for atom in rot_atoms:
            atom.x = R.apply(atom.x, diff)
        


    def get_omega(self,degree=False):
        chidx = self.chain.residues.index(self)
        if chidx == len(self.chain.residues)-1:
            return -999.
        else:
            next_mol = self.chain.residues[chidx+1]
        CA,C = self.fetchm(['CA','C'])
        N,CA2 = next_mol.fetchm(['N','CA'])
        dih = CA.dihedral(C,N,CA2) 
        if not degree:
            return dih
        else:
            return dih*180/pi

    def set_omega(self, degree ):
        phi = self.get_omega()
        diff = degree*pi/180. - phi
        chidx = self.chain.residues.index(self)
        if chidx == len(self.chain.residues)-1:
            return
        next_mol = self.chain.residues[chidx+1]
        C = self['C']
        N = next_mol['N']
        R = Rotation( C.x, N.x )
        rot_atoms = next_mol.fetch_atoms('N', inv = True)
        chidx = self.chain.residues.index(next_mol)
        if chidx != len(self.chain.residues)-1:
            for r in self.chain.residues[chidx+1:]:
                for atom in r.atoms:
                    rot_atoms.append( atom )
        for atom in rot_atoms:
            atom.x = R.apply(atom.x, diff)

    def set_omega_down(self, degree):
        phi = self.get_omega()
        diff = degree*pi/180. - phi
        chidx = self.chain.residues.index(self)
        if chidx == 0:
            return
        next_mol = self.chain.residues[chidx+1]
        C = self['C']
        N = next_mol['N']
        R = Rotation( N.x, C.x )
        rot_atoms = self.fetch_atoms('C', inv = True)
        chidx = self.chain.residues.index(self)
        if chidx != 0:
            for r in self.chain.residues[:chidx]:
                for atom in r.atoms:
                    rot_atoms.append( atom )
        for atom in rot_atoms:
            atom.x = R.apply(atom.x, diff)
        

    def nchi(self):
        self.get_real_resname()
        try:
            return len(_aa_chi[self.real_resname])
        except:
            return 0

    def get_chi(self, chi, degree = False ):
        if chi > self.nchi(): return -999
        dih_atoms = self.fetchm( _aa_chi[self.real_resname][chi][0] )
        dih = dih_atoms[0].dihedral( dih_atoms[1], dih_atoms[2], dih_atoms[3])
        if not degree:
            return dih
        else:
            return dih*180./pi

    def set_chi(self, chi, phi):
        if chi > self.nchi(): return 
        ang = self.get_chi(chi)
        dih_atoms = self.fetchm( _aa_chi[self.real_resname][chi][0] )
        rot_atoms = self.fetch_atoms( _aa_chi[self.real_resname][chi][1] )
        delta = phi/180*pi - ang
        r = Rotation( dih_atoms[1].x, dih_atoms[2].x )
        for atom in rot_atoms:
            atom.x = r.apply( atom.x, delta )

    def set_conformation(self, rotamer):
        self.get_real_resname()
        nchi = len(_aa_chi[self.real_resname] )
        for chi in range(nchi):
            self.set_chi(chi+1, rotamer[chi+1])

        
    def set_resname(self,resname):
        self.resname = resname
        for atom in self.atoms:
            atom.resname = resname
            
    def set_resid(self,resid):
        self.id = resid
        for atom in self.atoms:
            atom.resnr = resid

    def set_orig_resid(self,resid):
        if self.orig_id == 0:
            self.orig_id = resid

    def set_chain(self,chain):
        self.chain = chain
        for atom in self.atoms:
            atom.chain = chain
            
    def set_molecule(self):
        for atom in self.atoms:
            atom.molecule = self
            
    def set_chain_id(self,chain_id):
        self.chain_id = chain_id
        for atom in self.atoms:
            atom.chain_id = chain_id
#        self.chain.id = chain_id

    def insert_atom(self,pos,atom,id=True):
        """ insert atom at a certain position"""
        if pos not in range(len(self.atoms)+1):
            print 'Molecule has only %d atoms' % len(self.atoms)
            return
        else:
            if id:
                atom.resnr = self.id
                atom.resname = self.resname
                atom.chain_id = self.chain_id
            if pos == len(self.atoms):
                at = self.atoms[-1]
                if self.model is not None:
                    idx_model = self.model.atoms.index(at)+1
                if self.chain is not None:
                    idx_chain = self.chain.atoms.index(at)+1
                
            else:
                at = self.atoms[pos]
                if self.model is not None:
                    idx_model = self.model.atoms.index(at)
                if self.chain is not None:
                    idx_chain = self.chain.atoms.index(at)
            atom.molecule = self
            atom.chain = self.chain
            atom.model = self.model
            if pos == len(self.atoms):
                self.atoms.append(atom)
            else:
                self.atoms.insert(pos,atom)
            if self.model is not None:
                self.model.atoms.insert(idx_model,atom)
            if self.chain is not None:
                self.chain.atoms.insert(idx_chain,atom)
            if self.model is not None:
                self.model.renumber_atoms()

    def fetch(self,key,how='byname',wildcard=False):
        """ select atoms by name or by element"""
        result = []
        if how == 'byname':
            if not wildcard:
                for atom in self.atoms:
                    if atom.name==key:
                        result.append(atom)
            else:
                for atom in self.atoms:
                    if key in atom.name:
                        result.append(atom)
        elif how == 'byelem':
            for atom in self.atoms:
                if atom.symbol == key:
                    result.append(atom)
        return result
    
    def fetchm(self,keys,how='byname'):
        """select list of atom by name or element"""
        result = []
        if how=='byname':
            for key in keys:
                for atom in self.atoms:
                    if atom.name == key:
                        result.append(atom)
        elif how=='byelem':
            for atom in self.atoms:
                if atom.symbol in keys:
                    result.append(atom)
        return result

    def remove_atom(self,atom):
        """ delete atom from molecule"""
        if self.chain is not None:
            have_chain = True
        else: have_chain = False
        if self.model is not None:
            have_model = True
        else: have_model = False
        
        aidx = self.atoms.index(atom)
        if have_chain:
            chidx = self.chain.atoms.index(atom)
        if have_model:
            modidx = self.model.atoms.index(atom)
        del self.atoms[aidx]
        if have_chain:
            del self.chain.atoms[chidx]
        if have_model:
            del self.model.atoms[modidx]
            self.model.renumber_atoms()

    def append(self,atom):
        """ attach atom at the end"""
        if not isinstance(atom,Atom):
            raise TypeError, "%s is not an Atom instance" % str(atom) 
        else:
            n = len(self.atoms)
            if n == 0:
                self.atoms.append(atom)
            else:
                self.insert_atom(n,atom)

    def copy(self):
        """ copy molecule"""
        return copy.deepcopy(self)


    def get_bonded(self):
        bl = library._bonds[self.resname]
        for i, atom in enumerate(self.atoms[:-1]):
            for j, at in enumerate(self.atoms[i:]):
                n1 = atom.long_name.strip()
                n2 = at.long_name.strip()
                if (n1,n2) in bl or (n2,n1) in bl:
                    atom.bonds.append(at)
                    at.bonds.append(atom)
                    

    def get_mol2_types(self, nterminus = False):
        if not library._mol2_types.has_key(self.resname):
            print 'No mol2 lib entry for residue %s' % self.resname
            sys.exit(1)
        dic = library._mol2_types[self.resname]
        for atom in self.atoms:
            if atom.symbol == 'H':
                atom.atype = 'H'
                atom.q = 0
            else:
                if atom.name in ['OC1','OC2','O1','O2','OXT']:
                    atom.atype = 'O.co2'
                    atom.q = -.5
                elif atom.name == 'N' and nterminus:
                    atom.atype = 'N.4'
                    atom.q = 1
                else:
                    atom.atype = dic[atom.name][0]
                    atom.q = dic[atom.name][1]
                
    def is_protein_residue(self):
        if self.resname in library._protein_residues:
            return True
        else:
            return False


###################################################################################
# SMALL MOLECULE HANDLING
###################################################################################

# SD file

class SDMolecule:

    def __init__(self, lst):

        self.molfile = ''
        self.properties = {}
        self.name = ''
        self.name2 = ''
        self.read( lst )
        
    def read( self, lst ):
        self.name = lst[0].strip()
        self.name2 = lst[1].rstrip()
        s = ''
        for line in lst[2:]:
            if line.startswith('M  END'):
                s+=line.rstrip()
                break
            else:
                s+=line
        self.molfile = s
        idx = lst.index('M  END\n')
        for i, line in enumerate(lst[idx:]):
            if line.startswith('>'):
                prop = line.split()[1]
                self.properties[prop] = ''
                for line in lst[idx+i+1:]:
                    if line.startswith('>'): break
                    else:
                        self.properties[prop]+=line

    def write( self, fp ):
        print >>fp, self.name
        print >>fp, self.name2,
        print >>fp
        print >>fp, self.molfile
        for k, v in self.properties.items():
            print >>fp, '>', k
            print >>fp, v,
        print >>fp, '$$$$'
        

class SDFile:

    def __init__(self, filename=None):

        self.molecules = []
        if filename:
            self.read( filename )

    def read(self, fn ):
        #print >>sys.stderr, 'Reading SDFile:', fn
        l = open(fn).readlines()
        lst = []
        for line in l:
            if line.startswith('$$$$'):
                mol = SDMolecule( lst )
                self.molecules.append( mol )
                lst = []
            else:
                lst.append( line )

    def write(self, fp):
        for m in self.molecules:
            m.write( fp )


# mol2 file

class Mol2Molecule:

    def __init__(self, lines):
        self.atom_lines = []
        self.atoms = []
        self.bonds = []
        self.bond_lines = []
        self.header = []
        self.name = ''
        self.counts = []
        self.mol_type = ''
        self.charge_type = ''
        self.footer = []
        self.keys = []
        self.__get_keys(lines)
        self.read(lines)

        
    def __get_keys(self,lines):
        for line in lines:
            if line.startswith('@'):
                entr = line.strip()
                if entr not in ['@<TRIPOS>MOLECULE', '@<TRIPOS>ATOM', '@<TRIPOS>BOND']:
                    self.keys.append( entr )
    def read(self, lines):
        self.header = readUntil( lines, '@<TRIPOS>ATOM')
        self.atom_lines = readSection(lines,'@<TRIPOS>ATOM','@')
        self.bond_lines = readSection(lines,'@<TRIPOS>BOND','@')
        self.__parse_molecule()
        self.__parse_atoms()
        self.__parse_bonds()
        if self.keys:
            self.footer = [self.keys[0]]+readSection(lines, self.keys[0], 'XX')

    def __parse_molecule(self):
        lines = readSection(self.header, '@<TRIPOS>MOLECULE', '@')
        self.name = lines[0].strip()
        self.counts = [int(x) for x in lines[1].strip().split()]
        self.num_atoms = self.counts[0]
        self.num_bonds = self.counts[1]
        self.num_substr = self.counts[2]
        self.num_feat = self.counts[3]
        self.num_sets = self.counts[4]
        
        self.mol_type = lines[2].strip()
        self.charge_type = lines[3].strip()
        

    def __parse_atoms(self):
        for line in self.atom_lines:
            a = Atom( mol2line = line )
            self.atoms.append( a )

    def atom_by_id( self, idx ):
        for atom in self.atoms:
            if atom.id == idx: return atom
        return None

    def __parse_bonds(self):
        entries = parseList('iiis', self.bond_lines )
        for e in entries:
            atom1 = self.atom_by_id( e[1] )
            atom2 = self.atom_by_id( e[2] )
            bond_type = e[3]
            if atom1 is None or atom2 is None:
                print >>sys.stderr,'Mol2Molecule: Error in bond parsing'
                sys.exit(1)
            self.bonds.append( [atom1, atom2, bond_type] )

    def write(self, fp = sys.stdout):
        print >>fp, '@<TRIPOS>MOLECULE'
        print >>fp, self.name
        print >>fp, self.num_atoms, self.num_bonds, self.num_substr, self.num_feat, self.num_sets
        print >>fp, self.mol_type
        print >>fp, self.charge_type
        print >>fp
        print >>fp, '@<TRIPOS>ATOM'
        for atom in self.atoms:
            print >>fp, "%7d %-8s %9.4f %9.4f %9.4f %-6s %3d %-8s %9.4f" %\
                  (atom.id, atom.symbol, atom.x[0], atom.x[1], atom.x[2], atom.atype,\
                   atom.resnr, atom.resname, atom.q)
        print >>fp, '@<TRIPOS>BOND'
        for i, b in enumerate(self.bonds):
            print >>fp, "%6d %6d %6d %6s" % ((i+1), b[0].id, b[1].id, b[2])
        for line in self.footer:
            print >>fp, line
        
    def add_atom( self, atom ):
        n = len(self.atoms)+1
        atom.id = n
        atom.resnr = self.atoms[0].resnr
        atom.rename = self.atoms[0].resname
        self.num_atoms+=1
        self.atoms.append( atom )
        
    def add_bond(self, bond ):
        self.num_bonds+=1
        self.bonds.append( bond )
        

class Mol2File:

    def __init__(self, filename = None):

        self.molecules = []
        if filename:
            self.read( filename )


    def read(self,fn):
        l = open(fn).readlines()
        lst = []
        for i, line in enumerate(l):
            if line.startswith('@<TRIPOS>MOLECULE'):
                lines = readUntil(l[i+1:],'@<TRIPOS>MOLECULE')
                mol = Mol2Molecule( [line]+lines )
                self.molecules.append( mol )


    def write(self,fn=None):
        file_opened = False
        if fn:
            if hasattr(fn,"write"):
                fp = fn
            else:
                try:
                    fp = open(fn,'w')
                    file_opened = True
                except:
                    fp = sys.stdout
        else:
            fp = sys.stdout
        for m in self.molecules:
            m.write( fp )
        if file_opened:
            fp.close()
            
            
