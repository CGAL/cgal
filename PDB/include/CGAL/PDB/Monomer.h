/* Copyright 2004
   Stanford University

   This file is part of the DSR PDB Library.

   The DSR PDB Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DSR PDB Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA. */

#ifndef CGAL_DSR_PDB_RESIDUE_H
#define CGAL_DSR_PDB_RESIDUE_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Atom.h>
#include <CGAL/PDB/small_map.h>
#include <CGAL/Tools/Label.h>
#include <vector>
#include <map>

namespace CGAL { namespace PDB {


//! The class representing a residue or nucleotide.
/*!  All the information concerning atoms and bonds for each monomer
  is stored here. To add atoms to monomers, new monomers, or bonds
  to monomers, look in Monomer_data.cpp. There is documentation there
  of what you need to do.
*/
class Monomer {
  struct Monomer_type_tag{};
  struct Atom_type_tag{};
public:
  
  //! The labels for the types of residues.
  /*!
    \note These values must be packed, start from 0 and end with INV.
  */
  enum Type { GLY=0, ALA, VAL, LEU, ILE,
	      SER, THR, CYS, MET, PRO,
	      ASP, ASN, GLU, GLN, LYS,
	      ARG, HIS, PHE, TYR, TRP,
	      ACE, NH2, ADE, URA, CYT, GUA, THY, INV };


  //! The labels of atoms within residues
  /*!  These are the labels for each atom in each residue. The
    identifiers are attempting to following the PDB specs. Feel free to add more if needed.

    AL_N must be before AL_CA which must be before AL_C to get the backbone order correct.
  */
  enum Atom_key {AL_OTHER, AL_INVALID,
		 AL_N, AL_CA, AL_C, AL_O,

		 AL_H, AL_1H, AL_2H, AL_3H,
		      
		 AL_HA, AL_1HA, AL_2HA,
		     
		 AL_CB, AL_HB, AL_1HB, AL_2HB, AL_3HB,
		     
		 AL_OXT, AL_CH3,

		 AL_CG, AL_CG1, AL_CG2, AL_HG, AL_1HG, AL_2HG, //AL_HG1, 
		 AL_1HG1, AL_2HG1, AL_3HG1, AL_1HG2,
		     
		 AL_2HG2, AL_3HG2, AL_OG, AL_OG1, AL_SG,
		     
		 // AL_HD1, AL_HD2,
		 AL_CD, AL_CD1, AL_CD2,  AL_HD, AL_1HD, AL_2HD, AL_3HD, AL_1HD1,
		 AL_2HD1, AL_3HD1, AL_1HD2, AL_2HD2, AL_3HD2, AL_SD,
		 AL_OD1, AL_OD2, AL_ND1, AL_ND2,
		     
		 AL_CE, AL_CE1, AL_CE2, AL_CE3, AL_HE, AL_1HE, AL_2HE,
		 AL_3HE, //AL_HE1, AL_HE2, AL_HE3, 
		 AL_1HE2, AL_2HE2,
		 AL_OE1, AL_OE2, AL_NE, AL_NE1, AL_NE2,
		     
		 AL_CZ, AL_CZ2, AL_CZ3, AL_NZ, AL_HZ, AL_1HZ, AL_2HZ,
		 AL_3HZ, // AL_HZ2, AL_HZ3,
		     
		 AL_CH2, AL_NH1, AL_NH2, AL_OH, AL_HH, AL_1HH1,
		 AL_2HH1, AL_HH2, AL_1HH2, AL_2HH2, 
		 AL_1HH3, AL_2HH3, AL_3HH3,
    
		 AL_P, AL_OP1, AL_OP2,
		 AL_O5p, AL_C5p, AL_H5p, AL_H5pp,
		 AL_C4p, AL_H4p, AL_O4p, AL_C1p, AL_H1p,
		 AL_C3p, AL_H3p, AL_O3p, AL_C2p, AL_H2p,
		 AL_H2pp, AL_O2p, AL_HO2p,
		 
		 AL_N9, AL_C8, AL_H8,
		 AL_N7, AL_C5, AL_C4, AL_N3,
		 AL_C2, AL_H2, AL_N1, AL_C6, AL_N6,
		 AL_H61, AL_H62,

		 AL_O6, AL_H1, AL_N2, AL_H21, AL_H22,

		 AL_H6, AL_H5, AL_O2, AL_N4, AL_H41, AL_H42,
		 AL_H3, AL_O4, AL_C7, AL_H71, AL_H72, AL_H73,
		 AL_LAST_LABEL
  };



  //! This class defines the value_type of the Atom_iterator
  class Atom_pair: public small_map_value_type<Atom_key, Atom> {
    typedef small_map_value_type<Atom_key, Atom> P;
  public:
    Atom_pair(Atom_key k, const Atom &a=Atom()): P(k,a){}
    Atom_pair(){}
    const Atom &atom() const {return P::data();}
    Atom &atom() {return P::data();}
  };
  typedef small_map<Atom_pair> AtomsMap;


  //! Default constructor. Makes and invalid monomer.
  Monomer(){}

  //! Make a monomer of a given type
  Monomer(Type al);

  CGAL_COPY_CONSTRUCTOR(Monomer);

  void copy_from(const Monomer &o);

  void swap_with(Monomer &o);

  //! The label for the monomer
  Type type() const{
    return label_;
  }

  CGAL_ITERATOR(Atom, atom,  AtomsMap::const_iterator,
                AtomsMap::iterator, 
                atoms_.begin(),
		atoms_.end());

  CGAL_FIND(Atom,
            atoms_.find(fix_atom_key(k)),
            atoms_.end());

  CGAL_INSERT(Atom, insert_internal(k,m));

  Monomer::Atoms::iterator insert_internal(Atom_key k, const Atom &a);

  //! Return true if monomer of this type can have atoms of that type
  bool can_have_atom(Atom_key al) const;


  //! Set an atom using a string as a key
  void set_atom_from_string(const char *str, const Atom &a) {
    insert(atom_key(str), a);
  }

  //! Remove an atom from the monomer
  void erase_atom(Atom_key al);
  

  //! Write it for debugging
  void dump(std::ostream &out) const;

  //! Write it for debugging
  std::ostream& write(std::ostream &out) const {dump(out);return out;}


  //! Write the lines for a pdb file
  /*!
    Indices start at the start_index and the new start_index is returned.
  */
  int write(char chain, int monomer_index, char insert_code, int start_index, std::ostream &out) const;
   
    


 //! The Bond_iterator value_type is a pair of these
  class Bond_endpoint {
    Atoms::const_iterator aci_;
    Bond_endpoint(Atoms::const_iterator aci): aci_(aci){}
    friend class Monomer;
  public:
    Bond_endpoint(){}
    Atom_key key() const {return aci_->key();}
    const Atom &atom() const {return aci_->atom();}
  };



  //! A bond between two atoms in a monomer.
  /*!
    The ints refer the the atom index.
  */
  typedef std::pair<Bond_endpoint, Bond_endpoint> Bond;

  //! Return a list of all the bonds in the monomer
  CGAL_CONST_ITERATOR(Bond, bond, 
                      std::vector<Bond>::const_iterator,
                      bonds_.begin(),
                      bonds_.end());


  //! Return a point representing the sidechain
  /*!  If the sidechain is empty this returns the CA
    location. Otherwise it returns the location of some atom or the
    average of some atom locations.
  */
  Point sidechain_point() const;

    
  //! Set whether all the inter-atom bonds are present or not.
  /*!
    This must be true before Monomer::bonds_begin() is called.
  */
  void set_has_bonds(bool tf);
    
  //! Return whether this monomer has the inter-atom bonds computed.
  bool has_bonds() const {
    return !bonds_.empty();
  }


  //! Return the label of the first atom along the backbone in this monomer
  Atom_consts::const_iterator front_atom() const {
    if (is_amino_acid()) return find(AL_N);
    else return atoms().end();
  }
  //! Return the label of the last atom along the backbone in this monomer
  Atom_consts::const_iterator back_atom() const {
    if (is_amino_acid()) return find(AL_C);
    else return atoms().end();
  }

  //! Return the label of the first atom along the backbone in this monomer
  Atoms::iterator front_atom() {
    if (is_amino_acid()) return find(AL_N);
    else return atoms().end();
  }
  //! Return the label of the last atom along the backbone in this monomer
  Atoms::iterator back_atom() {
    if (is_amino_acid()) return find(AL_C);
    else return atoms().end();
  }

  //----- Static functions

  //! Convert a string for an amino acid type into a tag
  static Type type(const std::string &st);
  //! A string so you can write an amino acid type
  static std::string type_string(Type rl);

  //! Return the element corresponding to an atom label
  static Atom::Type element(Atom_key al);

  //! return the string corresponding to an atom key
  static std::string atom_key_string(Atom_key al);
  //! Return an atom type from a string
  /*!  Note, this type may be adjusted when the atoms is added to a
    monomer to address naming inconsistencies.
  */
  static Atom_key atom_key(const char *c);

  bool is_amino_acid() const {
    return label_ < ADE;
  }

private:

  Atom_key fix_atom_key(Atom_key atom_label) const;

  /*--------------- Data members-----------------------------------*/
  AtomsMap atoms_;
  std::vector<Bond> bonds_;
  Type label_;
};


//! Assign unique indices to all atoms in the monomer, starting at optional start value
/*!
  This returns the next unused index. 
*/
inline int index_atoms(const Monomer &m, int start=0) {
  CGAL_PDB_FOREACH(const Monomer::Atom_pair& a, m.atoms()) {
    a.atom().set_index(Atom::Index(start++));
  }
  return start;
}


CGAL_SWAP(Monomer)
CGAL_OUTPUT(Monomer)

}}
#endif
