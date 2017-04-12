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
This module contains the Chain class. It allows modifications
of structure files on the chain level:
Basic Usage:
     >>> model = Model().read(args['-f']) # read model
     >>> ch = model.chdic['A'] # select chain A

     The Chain instance contains:
     - ch.residues  -> list of residues
     - ch.atoms     -> list of atoms
     .
     .
     
     Some methods:
     >>> ch.get_sequence()    # return sequence in one-letter code
     >>> first_res = ch.residues[0].copy() # copy first residue
     >>> ch.append(first_res)    # attach residue at the end
     >>> ch.insert(5,first_res)  # insert residue at position 5
     .
     .
     
"""
from atomselection import *
from molecule import *
import copy, library
#import builder

class Chain(Atomselection):

    def __init__(self, seq = None, **kwargs):
        Atomselection.__init__(self)
        self.residues = []
        self.residue_tree_ok = False
        self.model = None
        self.nres = 0
        self.id = ''
        for key, val in kwargs.items():
            setattr(self,key,val)
        if self.residues:
            self.al_from_resl()
        if seq:
            self.create( seq )
            
    def __str__(self):
        s = '< Chain: id = %s nres = %d natoms = %d >' % \
            (self.id, len(self.residues), len(self.atoms))
        return s

    def __delitem__(self,item):
        res = self.residues[item]
        self.remove_residue(res)

        
    def get_sequence(self):
        """ returns the sequence as string (for fasta format)"""
        seq = ''
        for i, res in enumerate(self.residues):
            try:
                seq+=library._one_letter[res.resname]
                if (i+1) % 70 == 0: seq+='\n'
            except:
                seq+='-'+res.resname+'-'
        return seq

    def sequence(self):
        seq = ''
        for i, res in enumerate(self.residues):
            try:
                seq+=library._one_letter[res.resname]
            except:
                seq+='-'+res.resname+'-'
        return seq

    def al_from_resl(self):
        self.atoms = []
        for r in self.residues:
            r.chain = self
            for atom in r.atoms:
                self.atoms.append(atom)
                atom.chain = self
                
    def insert_residue(self,pos,mol,newResNum=False):
        if self.model is not None:
            have_model = True
        else:
            have_model = False
        if pos not in range(len(self.residues)+1):
            raise ValueError, 'Chain has only %d residues' % \
                  len(self.residues)
        else:
            mol.set_resid(-999)
	    if newResNum != False:
		mol.set_resid(newResNum)
            mol.set_chain_id(self.id)
            mol.chain = self
            if pos == len(self.residues):
                m = self.residues[-1]
                if have_model:
                    idx_model = self.model.residues.index(m)+1
            else:
                m = self.residues[pos]
                if have_model:
                    idx_model = self.model.residues.index(m)
                    mol.model = self.model
            if have_model:
                self.model.residues.insert(idx_model,mol)
                self.residues.insert(pos,mol)
		if newResNum==False:
                    self.model.renumber_residues()
                self.model.al_from_resl()
                self.model.renumber_atoms()
                self.al_from_resl()
            else:

                if pos == len(self.residues):
                    self.residues.append(mol)
                else:
                    self.residues.insert(pos,mol)
		if newResNum==False:
                    self.renumber_residues()
                self.al_from_resl()
                self.renumber_atoms()
        self.make_residue_tree()
        
    def renumber_residues(self):
        for i, res in enumerate(self.residues):
            res.set_resid(i+1)


    def insert_chain(self,pos,chain):
        idx_model = -1
        if self.model is not None:
            have_model = True
        else:
            have_model = False
        if hasattr(chain,"residues"):
            resl = chain.residues
        else:
            resl = chain
        if pos not in range(len(self.residues)+1):
            raise ValueError, 'Chain has only %d residues' % \
                  len(self.residues)
        if pos == len(self.residues):
            m = self.residues[-1]
            if have_model:
                idx_model = self.model.residues.index(m)+1
        else:
            m = self.residues[pos]
            if have_model:
                idx_model = self.model.residues.index(m)
        
        first = self.residues[:pos]
        last = self.residues[pos:]
        for res in resl:
            res.set_chain(self)
            res.set_chain_id(self.id)
        self.residues = first+resl+last
        self.atoms = []
        for r in self.residues:
            for atom in r.atoms:
                self.atoms.append(atom)
        if idx_model != -1:
            self.model.atoms = []
            self.model.residues = []
            for chain in self.model.chains:
                for atom in chain.atoms:
                    self.model.atoms.append(atom)
                for r in chain.residues:
                    self.model.residues.append(r)
            self.model.renumber_atoms()
            self.model.renumber_residues()
        else:
            self.renumber_atoms()
            self.renumber_residues()
        self.make_residue_tree()

    def remove_residue(self,residue,bKeepResNum=False):
        idx = self.residues.index(residue)
        try:
            midx = self.model.residues.index(residue)
        except:
            midx = -1
        del self.residues[idx]
        if midx != -1:
            del self.model.residues[midx]
            self.model.atoms = []
            for r in self.model.residues:
                for atom in r.atoms:
                    self.model.atoms.append(atom)
        self.atoms = []
        for r in self.residues:
            for atom in r.atoms:
                self.atoms.append(atom)
        if midx!=-1:
            self.model.renumber_atoms()
	    if bKeepResNum==False:
                self.model.renumber_residues()
        self.make_residue_tree()
        
    def replace_residue(self,residue, new, bKeepResNum=False):
        idx = self.residues.index(residue)
	if bKeepResNum==True:
            self.insert_residue(idx,new,residue.id)
	else:
	    self.insert_residue(idx,new)
        self.remove_residue(residue,bKeepResNum)

    def remove_atom(self,atom):
        m = atom.molecule
        m.remove_atom(atom)
        
    def fetch_residues(self,key, inv=False):
        if not hasattr(key,"append"):
            key = [key]
        result = []
        if not inv:
            for r in self.residues:
                if r.resname in key:
                    result.append(r)
        else:
            for r in self.residues:
                if r.resname not in key:
                    result.append(r)
        return result
    
    def set_chain_id(self,chain_id):
        old_id = self.id
        self.id = chain_id
        for r in self.residues:
            r.set_chain_id(chain_id)
        if self.model:
            self.model.chdic[chain_id] = self
            del self.model.chdic[old_id]
            

    def append(self,mol):
        if not isinstance(mol,Molecule):
            raise TypeError, "%s is not a Molecule instance" % str(mol)
        else:
            n = len(self.residues)
            self.insert_residue(n,mol)

    def copy(self):
        return copy.deepcopy(self)


    def get_bonded(self):
        for r in self.residues:
            r.get_bonded()
        for i, r in enumerate(self.residues[:-1]):
            r2 = self.residues[i+1]
            c = r.fetch('C')[0]
            n = r2.fetch('N')[0]
            d = n-c
            if self.unity=='A':
                d*=.1
            if d > 0.45:
                print 'Warning: Long Bond %d-%s <-> %d-%s' %\
                      (r.id,r.resname,r2.id,r2.resname)
            else:
                c.bonds.append(n)
                n.bonds.append(c)

    def get_mol2_types(self):
        r0 = self.residues[0]
        r0.get_mol2_types(nterminus = True)
        for r in self.residues[1:]:
            r.get_mol2_types()

    def get_mol2_resname(self):
        r0 = self.residues[0]
        r0.mol2_resname = r0.resname.lower()+'n+'
        r1 = self.residues[-1]
        r1.mol2_resname = r1.resname.lower()+'c-'
        for r in self.residues[1:-1]:
            if r.resname in ['LYS','LYP','LYSH','LYN']:
                r.mol2_resname = 'lys+'
            elif r.resname == 'ARG': r.mol2_resname = 'arg+'
            elif r.resname == 'GLU': r.mol2_resname = 'glu-'
            elif r.resname == 'ASP': r.mol2_resname = 'asp-'
            else: r.mol2_resname = r.resname.lower()
            

    def add_nterm_cap(self):
        if self.unity =='nm':
            self.nm2a()
            changed = True
        else:
            changed = False
        self.nbuild('GLY')
        m = self.nterminus()
        del m['H']
        n,ca,ha1,ha2  = m.fetchm(['N','CA','HA1','HA2'])
        ca.name = 'CH3'
        ha1.name = 'HH31'
        ha2.name = 'HH32'
        n.name = 'HH33'
        vec = array(ca.x)-array(n.x)
        x = 1./linalg.norm(vec)
        n.x = array(ca.x)-vec*x
        m.set_resname('ACE')        
        if changed:
            self.a2nm()

    def add_cterm_cap(self):
        if self.unity =='nm':
            self.nm2a()
            changed = True
        else:
            changed = False
        self.cbuild('GLY')
        m = self.cterminus() # new terminus
        del m['O']
        c,ca,ha1,ha2  = m.fetchm(['C','CA','HA1','HA2'])
        ca.name = 'CH3'
        ha1.name = 'HH31'
        ha2.name = 'HH32'
        c.name = 'HH33'
        vec = array(ca.x)-array(c.x)
        x = 1./linalg.norm(vec)
        c.x = array(ca.x)-vec*x
        m.set_resname('NME')
        if changed:
            self.a2nm()

        
        
    def attach_chain(self, newchain, phi = -139., psi = 135.):
        if self.unity =='nm':
            self.nm2a()
            changed = True
        else:
            changed = False

        cterm = self.residues[-1]
        ch = builder.attach_chain(cterm, newchain, phi = phi, psi = psi)
        self.insert_chain(len(self.residues), ch)
        if changed:
            self.a2nm()
        

    def __prepare_nterm_for_extension(self):
        nterm = self.nterminus()
        h = nterm.fetch('H')
        if h: return # ready
        else:
            if nterm.resname != 'PRO':
                try:
                    h1 = nterm.fetch('1H')[0] # will become H
                    h1.name = 'H'
                except:
                    h1 = nterm.fetch('H1')[0] # will become H
                    h1.name = 'H'
            else:
                try:
                    del nterm['H1']
                except:
                    try:
                        del nterm['1H']
                    except:
                        pass
            try:
                del nterm['2H']
            except:
                try:
                    del nterm['H2']
                except:
                    pass
            try:
                del nterm['3H']
            except:
                pass
            try:
                del nterm['H3']
            except:
                pass
            
    def __prepare_cterm_for_extension(self):
        cterm = self.cterminus()
        if  cterm.resname not in library._protein_residues:
            print " Cannot attach to this residue! ", cterm.resname
            sys.exit(1)
        a = cterm.fetch_atoms('OXT')
        if a:
            del cterm['OXT']
        a = cterm.fetch_atoms('O2')
        if a:
            del cterm['O2']
        a = cterm.fetch_atoms('O\'\'')
        if a:
            del cterm['O\'\'']
        a = cterm.fetch_atoms('O1')
        if a:
            a[0].name = 'O'
        a = cterm.fetch_atoms('O\'')
        if a:
            a[0].name = 'O'
        

    def cbuild(self, resn, next_phi = -139., psi = 135.):
        if len(resn) == 1:
            resn = library._aacids_dic[resn]
#        new = builder.make_residue( resn )
        new = Molecule().new_aa(resn)
        self.__prepare_cterm_for_extension()
        cterm = self.cterminus()
        Ca, C, O = cterm.fetchm(['CA','C','O'])
        N, CA2 = new.fetchm(['N','CA'])
        if self.unity != 'A':
            self.nm2a()
            changed = True
        else:
            changed = False
        l = 1.33 # length of peptide bond
        v = array(C.x)-array(Ca.x)
        x = 1./linalg.norm(v)
        vec = v*x*l
        newpos = array(C.x)+vec
        t = array(N.x)-newpos
        # shift to new position
        for a in new.atoms:
            a.x = array(a.x)-t
        # gen rotation vector
        # this rotation yields the correct
        # O-C-N angle
        v1 = array(C.x)-array(Ca.x)
        v2 = array(O.x)-array(C.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        cr = C.x + cr
        # do first rotation
        r = Rotation(cr,C.x)

        for atom in new.atoms:
            atom.x = r.apply(atom.x,pi/3.)

            
        # do next rotation
        # here we correct the C-N-Ca angle


        # special correction for proline
        if new.resname == 'PRO':
             N, CA2,H = new.fetchm(['N','CA','CD'])
             Ca, C, O, Nn = cterm.fetchm(['CA','C','O','N'])
             v1 = array(N.x)-array(CA2.x)
             v2 = array(H.x)-array(N.x)
             cr = cross(v1,v2)
             x = 1./linalg.norm(cr)
             cr = cr*x
             v3 = array(Ca.x)-array(C.x)
             v4 = array(C.x)-array(O.x)
             cr2 = cross(v3,v4)
             x = 1./linalg.norm(cr2)
             cr2 = cr2*x
             d = dot(cr,cr2)
             dd = arccos(dot(cr,cr2))
             v5 = cross(cr,cr2)
             x = 1./linalg.norm(v5)
             v5 = v5*x
             rv = N.x + v5
             r = Rotation(rv,N.x)
             for atom in new.atoms:
                 if atom.name !='N':
                     atom.x = r.apply(atom.x,-dd)


             
        Ca, C, O = cterm.fetchm(['CA','C','O'])
        N, CA2 = new.fetchm(['N','CA'])
        v1 = array(N.x)-array(C.x)
        v2 = array(CA2.x)-array(N.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        cr = N.x + cr
        an = N.angle(C,CA2)

        delta = an - 120*pi/180

        # do second rotation
        r = Rotation(cr,N.x)
        for atom in new.atoms:
            if atom.name !='N':
                atom.x = r.apply(atom.x,-delta)

        an = N.angle(C,CA2)
         
        if new.resname == 'PRO':
            CD = new.fetch('CD')[0]
            aa = N.angle(CD,C, degree=True)
            if aa < 100:
                r = Rotation( N.x, CA2.x )
                for atom in new.atoms:
                    if atom.name not in ['N','CA']:
                        atom.x = r.apply( atom.x , pi )

        # correct H if we have one
        H = new.fetch('H')
        if H:
            H=H[0]
            # we place the H-atom along the
            # N-Ca vector.
            # then we move it
            v1 = array(CA2.x)-array(N.x)
            x = 1./linalg.norm(v1)
            H.x = N.x + v1*x
            H.x = r.apply(H.x,2*pi/3.)
        self.append( new )
        cterm.set_omega( 180 )
        cterm.set_psi( psi, propagate = True )
        new.set_phi( next_phi )
        if changed:
            self.a2nm()


        
    def nbuild(self, resn, phi = -139., psi = 135.): 

        if len(resn) == 1:
            resn = library._aacids_dic[resn]
#        new = builder.make_residue(resn,hydrogens = True)
        new = Molecule().new_aa(resn)
        
        # we need Ca, C and O from nterm
        Ca, C, O, Nn = new.fetchm(['CA','C','O','N'])
        # and N from new
        self.__prepare_nterm_for_extension()
        nterm = self.nterminus()

        if nterm.resname=='PRO':
            N, CA2,H = nterm.fetchm(['N','CA','CD'])
        else:
            try:
                N, CA2,H = nterm.fetchm(['N','CA','H'])
            except:
                N, CA2,H = nterm.fetchm(['N','CA','CN'])

        l = 1.33 # length of peptide bond

        v = array(N.x)-array(CA2.x)
        x = 1./linalg.norm(v)
        vec = v*x*l
        newpos = array(N.x)+vec
        t = array(C.x)-newpos
        
        # shift to new position
        for a in new.atoms:
            a.x = array(a.x)-t
        # get Ca,C,O and N,H,CA2 in plane
        
        v1 = array(N.x)-array(CA2.x)
        v2 = array(H.x)-array(N.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        
        v3 = array(Ca.x)-array(C.x)
        v4 = array(C.x)-array(O.x)
        cr2 = cross(v3,v4)
        x = 1./linalg.norm(cr2)
        cr2 = cr2*x
        dd = arccos(dot(cr,cr2))
        v5 = cross(cr,cr2)
        x = 1./linalg.norm(v5)
        v5 = v5*x
        rv = C.x + v5
        r = Rotation(rv,C.x)
        for atom in new.atoms:
            if atom.name !='C':
                atom.x = r.apply(atom.x,dd)
        # correct angle fo N,H,C
        v1 = array(N.x)-array(CA2.x)
        v2 = array(H.x)-array(N.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        cr = N.x + cr
        ang = N.angle(H,C)
        dd = 2*pi/3. - ang
        # do first rotation
        r = Rotation(cr,N.x)
        for atom in new.atoms:
            atom.x = r.apply(atom.x,dd)
        # correct angle fo C,N,Ca
        v1 = array(C.x)-array(N.x)
        v2 = array(Ca.x)-array(C.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        cr = C.x + cr
        ang = C.angle(Ca,N)
        dd = 2*pi/3. - ang

        r = Rotation(cr,C.x)
        for atom in new.atoms:
            if atom.name !='C':
                atom.x = r.apply(atom.x,dd)
        dih = CA2.dihedral(N,C,Ca) 
        delta = pi-dih
        r = Rotation(N.x,C.x)
        for atom in new.atoms:
            if atom.name !='C':
                atom.x = r.apply(atom.x,delta)
        ang = Nn.dihedral(Ca,C,N) 
        dd = pi - ang
        r = Rotation(C.x,Ca.x)
        for atom in new.atoms:
            if atom.name !='C':
                atom.x = r.apply(atom.x,dd)


        dih = H.dihedral(N,C,O) 
        delta = pi-dih
        r = Rotation(N.x,C.x)
        O.x = r.apply(O.x,delta)

        self.insert_residue(0, new )
        new.set_omega_down( 180 )
        nterm.set_phi_down( phi, True)
        new.set_psi_down( psi, True)


    def create(self, seq, phi_psi = [], ss = []):
##         print 'seq = ', seq
##         print 'phi_psi', phi_psi
##         if ss:
##             assert len(seq) == len(ss)
##             phi_psi = []
##             for i in range(len(seq)):
##                 if ss[i] == 'H':
##                     phi_psi.append( [-57,-47] )
##                 elif ss[i] == 'E':
##                     phi_psi.append( [-139, 135] )
                    
##         if not phi_psi:
##             print seq, len(seq)
##             for i in range(len(seq)):
##                 print i
##                 phi_psi.append( [-139, 135] )
##         else:
##             print phi_psi
##             print seq, phi_psi, len(seq), len(phi_psi)
##             assert len(phi_psi) == len(seq)
        first_res = seq[0]
        resn = library._aacids_dic[first_res]
#        new = builder.make_residue(resn,hydrogens = True)
        new = Molecule().new_aa(resn)
        
        self.residues.append( new )
        new.chain = self
##         new.set_phi( phi_psi[0][0] )
##         new.set_psi( phi_psi[0][1] )
        new.unity = 'A'
        self.al_from_resl()
        for i, aa in enumerate(seq[1:]):
            self.cbuild(aa) #, phi_psi[i+1][0], phi_psi[i+1][1])
        return self
        

    def fuse(self, new, phi=-139, psi=135 ):
        if new.unity != 'A': newchain.nm2a()
        nterm = new.nterminus()
        self.__prepare_cterm_for_extension()
        cterm = self.cterminus()
        Ca, C, O = cterm.fetchm(['CA','C','O'])
        N, CA2 = nterm.fetchm(['N','CA'])
        if self.unity != 'A':
            self.nm2a()
            changed = True
        else:
            changed = False
        l = 1.33 # length of peptide bond
        v = array(C.x)-array(Ca.x)
        x = 1./linalg.norm(v)
        vec = v*x*l
        newpos = array(C.x)+vec
        t = array(N.x)-newpos
        # shift to new position
        for a in new.atoms:
            a.x = array(a.x)-t
        # gen rotation vector
        # this rotation yields the correct
        # O-C-N angle
        v1 = array(C.x)-array(Ca.x)
        v2 = array(O.x)-array(C.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        cr = C.x + cr
        # do first rotation
        r = Rotation(cr,C.x)

        for atom in new.atoms:
            atom.x = r.apply(atom.x,pi/3.)
        # do next rotation
        # here we correct the C-N-Ca angle


        # special correction for proline
        if nterm.resname == 'PRO':
             N, CA2,H = nterm.fetchm(['N','CA','CD'])
             Ca, C, O, Nn = cterm.fetchm(['CA','C','O','N'])
             v1 = array(N.x)-array(CA2.x)
             v2 = array(H.x)-array(N.x)
             cr = cross(v1,v2)
             x = 1./linalg.norm(cr)
             cr = cr*x
             v3 = array(Ca.x)-array(C.x)
             v4 = array(C.x)-array(O.x)
             cr2 = cross(v3,v4)
             x = 1./linalg.norm(cr2)
             cr2 = cr2*x
             d = dot(cr,cr2)
             dd = arccos(dot(cr,cr2))
             v5 = cross(cr,cr2)
             x = 1./linalg.norm(v5)
             v5 = v5*x
             rv = N.x + v5
             r = Rotation(rv,N.x)
             for atom in new.atoms:
                 if atom.name !='N' or atom.molecule != nterm:
                     atom.x = r.apply(atom.x,-dd)

        Ca, C, O = cterm.fetchm(['CA','C','O'])
        N, CA2 = nterm.fetchm(['N','CA'])
        v1 = array(N.x)-array(C.x)
        v2 = array(CA2.x)-array(N.x)
        cr = cross(v1,v2)
        x = 1./linalg.norm(cr)
        cr = cr*x
        cr = N.x + cr
        an = N.angle(C,CA2)

        delta = an - 120*pi/180
        # do second rotation
        r = Rotation(cr,N.x)
        for atom in new.atoms:
            if atom.name !='N' or atom.molecule != nterm:
                atom.x = r.apply(atom.x,-delta)

        if nterm.resname == 'PRO':
            CD = nterm.fetch('CD')[0]
            aa = N.angle(CD,C, degree=True)
            if aa < 100:
                if delta < 0:
                    aa_new = -2*pi/3.
                else:
                    aa_new = 2*pi/3.
                for atom in new.atoms:
                    if atom.name !='N' or atom.molecule != nterm:
                        atom.x = r.apply(atom.x,-aa_new)
        # correct H if we have one
        H = nterm.fetch('H')
        if H:
            H=H[0]
            # we place the H-atom along the
            # N-Ca vector.
            # then we move it
            v1 = array(CA2.x)-array(N.x)
            x = 1./linalg.norm(v1)
            H.x = N.x + v1*x
            H.x = r.apply(H.x,2*pi/3.)
        self.insert_chain( len(self.residues), new )
        cterm.set_omega( 180 )
        cterm.set_psi( psi, propagate = True )
        nterm.set_phi( phi, True )
        if changed:
            self.a2nm()

    def make_residue_tree(self):
        self.residue_tree_ok = True
        for i, r in enumerate(self.residues[:-1]):
            if r.is_protein_residue():
                next_res = self.residues[i+1]
                if next_res.is_protein_residue():
                    C = r['C']
                    N = next_res['N']
                    d = C - N
                    if N.unity == 'nm':
                        d*=10
                    if d < 1.7: # bond
                        r.next = next_res
                        next_res.previous = r
                    else:
                        print >>sys.stderr, 'Gap between residues ', r, '< - >', next_res, 'dist = ', d
                        self.residue_tree_ok = False
                            
    def cterminus(self):
        i = len(self.residues) - 1
        r = self.residues[i]
        if r.is_protein_residue(): return r
        while i != 0:
            i-=1
            r = self.residues[i]
            if r.is_protein_residue(): return r
        return None

    def nterminus(self):
        i = 0
        r = self.residues[0]
        if r.is_protein_residue(): return r
        while i<len(self.residues):
            i+=1
            r = self.residues[i]
            if r.is_protein_residue(): return r
        return None
    
    def rename_atoms(self):
        for atom in self.atoms:
            atom.make_long_name()
            atom.name = atom.long_name.strip()
        cterm = self.cterminus()
        if cterm is not None:
            for atom in cterm.atoms:
                if atom.name == 'O1' or atom.name == 'OC1':
                    atom.name = 'O'
                    atom.make_long_name()
                elif atom.name == 'O2' or atom.name == 'OC2':
                    atom.name = 'OXT'
                    atom.make_long_name()

    def residue(self, idx):
        return self.residues[idx-1]



