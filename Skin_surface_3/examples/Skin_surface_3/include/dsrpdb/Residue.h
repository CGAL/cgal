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
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef DSR_PDB_RESIDUE_H
#define DSR_PDB_RESIDUE_H

#include <dsrpdb/Point.h>
//#include <dsr/pdb/Protein.h>
#include <dsrpdb/Atom.h>
#include <vector>
#include <map>
#include <dsrpdb/small_map.h>
//#include <dsrpdb/label.h>

namespace dsrpdb {

  class Protein;

  //! The class representing a residue.
  /*!  All the information concerning atoms and bonds for each residue
    is stored here. To add atoms to residues, new residues, or bonds
    to residues, look in Residue_data.cc. There is documentation there
    of what you need to do.
  */
  class Residue {
    friend class Protein;
    struct Residue_type_tag{};
    struct Atom_type_tag{};
  public:
    //! The type for storing residue indices in the PDB
    typedef PDB_index<Residue> Index;

    //! The labels for the types of residues.
    enum Type { GLY=0, ALA, VAL, LEU, ILE,
		SER, THR, CYS, MET, PRO,
		ASP, ASN, GLU, GLN, LYS,
		ARG, HIS, PHE, TYR, TRP,
		ACE, NH2, INV };


    //! The labels of atoms within residues
    /*!  These are the labels for each atom in each residue. The
      identifiers are attempting to following the PDB specs. Feel free to add more if needed.

      AL_N must be before AL_CA which must be before AL_C to get the backbone order correct.
    */
    enum Atom_label {AL_OTHER, AL_INVALID,
		     AL_N, AL_H, AL_1H, AL_2H, AL_3H,
		     
		     
		     AL_CA, AL_HA, AL_1HA, AL_2HA,
		     
		     AL_CB, AL_HB, AL_1HB, AL_2HB, AL_3HB,
		     
		     AL_C, AL_O, AL_OXT, AL_CH3,

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
		     AL_1HH3, AL_2HH3, AL_3HH3};



    //! A bond between two atoms in a residue.
    /*!
      The ints refer the the atom index.
    */
    typedef std::pair<Atom::Index,Atom::Index> Bond;

    typedef small_map<Atom_label, Atom> Atoms;

    //! Default constructor. Makes and invalid residue.
    Residue(){}

    //! Make a residue of a given type
    Residue(Type al);

    //! The label for the residue
    Type type() const;


    //! Return a list of all the bonds in the residue
    typedef std::vector<Bond>::const_iterator Bonds_iterator;
 
    //! Begin iterating through the list of all the bonds
    /*!
      Note that the iterator will be invalidated if the residue is changed.

      Note that if Residue::has_bonds() is false, this returns an empty sequence.
     */
    Bonds_iterator bonds_begin() const {
      return bonds_.begin();
    }
    //! End bond iteration
    Bonds_iterator bonds_end() const {
      return bonds_.end();
    }
    //! The number of atoms present in the residue
    unsigned int number_of_bonds() const;


    //! An iterator to list all the atoms
    typedef Atoms::iterator Atoms_iterator;

    //! Return a list of the labels of all the atoms which are present
    Atoms_iterator atoms_begin() {
      return atoms_.begin();
    }

    //! End iterating through the atoms
    Atoms_iterator atoms_end() {
      return atoms_.end();
    }

    //! An iterator to list all the atoms
    typedef Atoms::const_iterator Const_atoms_iterator;

    //! Return a list of the labels of all the atoms which are present
    Const_atoms_iterator atoms_begin() const {
      return atoms_.begin();
    }

    //! End iterating through the atoms
    Const_atoms_iterator atoms_end() const {
      return atoms_.end();
    }



    //! The number of atoms present in the residue
    unsigned int number_of_atoms() const;


    //! Return true if the atom is in the atoms() list
    bool has_atom(Atom_label al) const;

    //! Return true if residues of this type can have atoms of that type
    bool can_have_atom(Atom_label al) const;

    
    //! Return the data for an atom
    const Atom& atom(Atom_label al) const;

    //! Return the label of the atom with this index.
    Atom_label atom_label(Atom::Index model_index) const;

    //! Set an atom
    /*!
      If the atom is not already there, then the bonds iterators are invalidated. 
    */
    void set_atom(Atom_label al, const Atom &a);

    //! Set an atom using a string as a label
    void set_atom_from_string(const char *str, const Atom &a) {
      set_atom(atom_label(str), a);
    }
    //! The index of the last atom in the residue
    Atom::Index last_atom_index() const;
   
    //! The index for the residue
    /*!
      This is  0 based index so it is the PDB index -1.
    */
    Index index() const {
      return index_;
    }

    //! Set the index for the residue
    void set_index(Index i);
    

    //! Write it for debugging
    void dump(std::ostream &out) const;
    //! Write the lines for a pdb file
    void write(char chain, std::ostream &out) const;
   
    

    //! Return a point representing the sidechain
    /*!  If the sidechain is empty this returns the CA
      location. Otherwise it returns the location of some atom or the
      average of some atom locations.
    */
    Point sidechain_point() const;

    
    //! Set whether all the inter-atom bonds are present or not.
    /*!
      This must be true before Residue::bonds_begin() is called.
    */
    void set_has_bonds(bool tf);
    
    //! Return whether this residue has the inter-atom bonds computed.
    bool has_bonds() const {
      return !bonds_.empty();
    }


    //----- Static functions

    //! Convert a string for an amino acid type into a tag
    static Type type(const std::string &st);
    //! A string so you can write an amino acid type
    static std::string type_string(Type rl);

    //! Return the element corresponding to an atom label
    static Atom::Type element(Atom_label al);

    //! return the string corresponding to an atom label
    static std::string atom_label_string(Atom_label al);
    //! Return an atom label from a string
    /*!  Note, this label may be adjusted when the atoms is added to a
      residue to address naming inconsistencies.
    */
    static Atom_label atom_label(const char *c);

   


  protected:
    /*Residue( Label al, 
	     std::vector<Residue::Atom_label> * atoms,
	     std::vector<std::pair<Residue::Atom_label,
	     Residue::Atom_label> > *bonds);*/
    Atom::Index index(Residue::Atom_label al) const;

    Atom::Index min_atom_index() const {
      return min_atom_index_;
    }

    Atoms_iterator atoms_iterator_from_index(Atom::Index ind);
    Const_atoms_iterator atoms_iterator_from_index(Atom::Index ind) const;
  private:

    /*--------------- Data members-----------------------------------*/
    Atoms atoms_;
    std::vector<Bond> bonds_;
    Type label_;
    Index index_;
    Atom::Index min_atom_index_;
  };
}
#endif
