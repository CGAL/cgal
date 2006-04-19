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

#ifndef DSR_PDB_PROTEIN_H
#define DSR_PDB_PROTEIN_H

#include <dsrpdb/Point.h>
#include <dsrpdb/Atom.h>
#include <dsrpdb/Residue.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>


namespace dsrpdb {

  class Model;
  
  //! A class representing a single chain of a protein.  
  /*!
     It probably should be called chain, but it is not :-)
  */
  class Protein {
    friend class Model;
  public:

    /*!
       \example check_protein.cc
    */

    //! A chemical bond within the protein
    typedef std::pair<Atom::Index,Atom::Index> Bond;

    /*! Build a protein from a pdb file.
      See check_protein.cc for an example of using this to read a
      pdb file.

      The optional bool controls whether errors (such as unparsable
      PDB lines). Set it to false to disable printing errors.
    */
    Protein(std::istream &in, bool print_errors=false);

    //! Default
    Protein();

    //! Return the chain id
    char chain() const;

    //! Set the chain id. 
    void set_chain(char c);
    
    //! An iterator to go through all the residues
    typedef std::vector<Residue>::iterator Residues_iterator;
    //! Iterate through all the residues
    Residues_iterator residues_begin() {
      return residues_.begin();
    }
    //! End iterator through the residues.
    Residues_iterator residues_end() {
      return residues_.end();
    }
    //! Iterator for residues in order of addition.
    typedef std::vector<Residue>::const_iterator Const_residues_iterator;
    //! Begin iterating residues.
    Const_residues_iterator residues_begin() const {
      return residues_.begin();
    }
    //! End iterating residues.
    Const_residues_iterator residues_end() const {
      return residues_.end();
    }
    unsigned int number_of_residues() const {
      return residues_.size();
    }

    void new_residue(const Residue &res);
    
    //! An iterator to iterate through all the atoms of the protein
    class Atoms_iterator;
    //! Begin iterating through the atoms
    Atoms_iterator atoms_begin();
    //! End iterating through the atoms.
    Atoms_iterator atoms_end();

    //! An iterator to iterate through all the atoms of the protein
    class Const_atoms_iterator;
    //! Begin iterating through the atoms
    Const_atoms_iterator atoms_begin() const;
    //! End iterating through the atoms.
    Const_atoms_iterator atoms_end() const;
    //! This is non-const time.
    unsigned int number_of_atoms() const;

    //! An iterator to iterate through all the bonds of the protein
    class Bonds_iterator;
    /*! Begin iterating through all the bonds
      \note These bonds are indexed by the atom numbers in the pdb
      which will not correspond to the atom sequence numbers in the
      atoms_begin() sequence when the PDB atoms do not start from 1 or
      have missing residues.  This is a bug (I think) and will be
      fixed. However, doing so is slightly complicated as I need to
      handle the case of missing residues.
    */
    Bonds_iterator bonds_begin() const;
    //! End iterating through all the bonds.
    Bonds_iterator bonds_end() const;
    //! This is non-const time.
    unsigned int number_of_bonds() const;

    //! Return the dsrpdb::Residue which contains the atoms with a particular index.
    /*!  This operation is currently linear in the number of
      residues. Making it logrighmic should be easy, but it isn't
      done.
    */
    const Residue& residue_containing_atom(Atom::Index atom_index) const;
    //! Return the dsrpdb::Residue which contains the atom witha particular index.
    Residue& residue_containing_atom(Atom::Index atom_index);

    //! Return the atom which has the pdb index passed.
    const Atom& atom(Atom::Index atom_index) const;

    //! Return the atom which has the pdb index passed.
    void set_atom(Atom::Index atom_index, const Atom &a);


    //! The sequence of residue types.
    std::vector<Residue::Type> sequence() const;

    //! Write as part of pdb file.
    void write(std::ostream &out) const;
    //! Write a pdb file. 
    /*!
      See check_protein.cc for an example of using this to write a
     pdb file.
    */
    void write_pdb(std::ostream &out) const;
    //! Dump as human readable.
    void dump(std::ostream &out) const;


    //! Return the residue with a particular index.
    /*!  \note This function has changed. Use the iterators if you
      just want to go through the list.
    */
    const Residue& residue(Residue::Index i) const;

    //! Return true if there is a residue with that index
    bool has_residue(Residue::Index i) const;

    /*const Atom &atom(unsigned int i) const {
      for (Const_residues_iterator it = residues_begin(); it != residues_.end(); ++it){
	
      }
      }*/


    //! Return whether bonds have been computed for this protein. 
    bool has_bonds() const {
      return residues_[0].has_bonds();
    }

    //! Set whether the protein has bonds or not.
    void set_has_bonds(bool tf) {
      for (unsigned int i=0; i< residues_.size(); ++i){
	residues_[i].set_has_bonds(tf);
      }
    }


#if 0
    //! An interator through the backbone coordinates.
    /*!
      The value_type is a dsr::Point. 
    */
    class Backbone_coordinates_iterator;
    //! Begin iterating through the backbone coordinates.
    Backbone_coordinates_iterator backbone_coordinates_begin() const;
    //! End iterating through the backbone coordinates.
    Backbone_coordinates_iterator backbone_coordinates_end() const;
#endif
    
  protected:
    void process_line(const char *line);
    
    unsigned int residue_offset_of_atom_index(Atom::Index i) const;

    unsigned int residue_offset(Residue::Index i) const;
    
    std::vector<Residue> residues_;
    std::vector<std::string> header_;
    //static Residue dummy_residue_;
    char chain_;
  };







  //! An iterator through the atoms of a dsrpdb::Protein.
  class Protein::Atoms_iterator {
    friend class Protein;
    friend class Protein::Const_atoms_iterator;
  public:
    //! The value_type is a dsrpdb::Atom
    typedef Atom value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t difference_type;
    typedef Residue::Atoms_iterator::reference reference;
    typedef Residue::Atoms_iterator::pointer pointer;

    reference operator*() {
      return *ait_;
    }
    pointer operator->() {
      return ait_.operator->();
    }
    Atoms_iterator operator++() {
      ++ait_;
      if (ait_== aend_) {
	++rit_;
	if (rit_!= rend_){
	  ait_=rit_->atoms_begin();
	  aend_= rit_->atoms_end();
	}
      }
      return *this;
    }
    
    bool operator==(const Atoms_iterator& o) const {
      if (rit_ == rend_) return rit_==o.rit_;
      else return rit_== o.rit_ && ait_ == o.ait_;
    }
    bool operator!=(const Atoms_iterator& o) const {
      if (rit_== rend_) return rit_!= o.rit_;
      else return rit_!= o.rit_ || ait_ != o.ait_;
    }

  protected:
    Atoms_iterator(std::vector<Residue>::iterator b, 
		   std::vector<Residue>::iterator e): rit_(b), rend_(e){
      if (b != e) {
	ait_= rit_->atoms_begin();
	aend_= rit_->atoms_end();
      }
    }
    std::vector<Residue>::iterator rit_, rend_;
    Residue::Atoms_iterator ait_, aend_;
  };


  //! An iterator through the atoms of a dsrpdb::Protein.
  class Protein::Const_atoms_iterator {
    friend class Protein;
  public:
    //! The value_type is a dsrpdb::Atom
    typedef Atom value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t difference_type;
    typedef Residue::Const_atoms_iterator::reference reference;
    typedef Residue::Const_atoms_iterator::pointer pointer;

    reference operator*() const {
      return *ait_;
    }
    pointer operator->() const {
      return ait_.operator->();
    }
    Const_atoms_iterator operator++() {
      ++ait_;
      if (ait_== aend_) {
	++rit_;
	if (rit_!= rend_){
	  ait_=rit_->atoms_begin();
	  aend_= rit_->atoms_end();
	}
      }
      return *this;
    }
    
    bool operator==(const Const_atoms_iterator& o) const {
      if (rit_ == rend_) return rit_==o.rit_;
      else return rit_== o.rit_ && ait_ == o.ait_;
    }
    bool operator!=(const Const_atoms_iterator& o) const {
      if (rit_== rend_) return rit_!= o.rit_;
      else return rit_!= o.rit_ || ait_ != o.ait_;
    }

    Const_atoms_iterator(Protein::Atoms_iterator it){
      rit_= it.rit_;
      rend_= it.rend_;
      ait_= it.ait_;
      aend_= it.aend_;
    }

  protected:
    Const_atoms_iterator(std::vector<Residue>::const_iterator b, 
			 std::vector<Residue>::const_iterator e): rit_(b), rend_(e){
      if (b != e) {
	ait_= rit_->atoms_begin();
	aend_= rit_->atoms_end();
      }
    }
    std::vector<Residue>::const_iterator rit_, rend_;
    Residue::Const_atoms_iterator ait_, aend_;
  };


  //! An iterator through the bonds of a dsrpdb::Protein
  class Protein::Bonds_iterator {
    friend class Protein;
  public:
    //! The value_type is a dsrpdb::Protein::Bond
    typedef Bond value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t difference_type;
    typedef const Bond& reference;
    typedef const Bond* pointer;
    
    reference operator*() const {
      return cur_bond();
    }
    pointer operator->() const {
      return &cur_bond();
    }
    Bonds_iterator operator++() {
      if (nc_) {
	nc_=false;
      } else {
	++ait_;
	if (ait_== aend_) {
	  std::vector<Residue>::const_iterator orit= rit_;
	  ++rit_;
	  if (rit_!= rend_){
	    if (static_cast<unsigned int>(orit->index()) +1 == static_cast<unsigned int>(rit_->index())
		&& orit->has_atom(Residue::AL_C) && rit_->has_atom(Residue::AL_N)){
	      nc_=true;
	    } /*else {
	      std::cout << "Skipping bond between residues " << orit->index() << " and " << rit_->index() << "\n";
	      std::cout << static_cast<unsigned int>(orit->index()) << " " << static_cast<unsigned int>(rit_->index()) << std::endl;
	      rit_->dump(std::cout);
	      orit->dump(std::cout);
	      }*/
	    ait_=rit_->bonds_begin();
	    aend_= rit_->bonds_end();
	  }
	}
      }
      return *this;
    }
     bool operator==(const Bonds_iterator& o) const {
       if (rit_ == rend_) return rit_==o.rit_;
       else return rit_== o.rit_ && ait_ == o.ait_ && nc_== o.nc_;
    }
    bool operator!=(const Bonds_iterator& o) const {
      if (rit_== rend_) return rit_!= o.rit_;
      else return rit_!= o.rit_ || ait_ != o.ait_ || nc_ != o.nc_;
    }
    
    Bonds_iterator(){}
  protected:
    Bonds_iterator(std::vector<Residue>::const_iterator b, 
		   std::vector<Residue>::const_iterator e): rit_(b), rend_(e), nc_(false){
      if (b != e) {
	ait_= rit_->bonds_begin();
	aend_= rit_->bonds_end();
      }
    }

    const Bond &cur_bond() const {
      if (nc_) {
	assert(rit_ != rend_);
	static Bond b;
	b= Bond((rit_-1)->atom(Residue::AL_C).index(),
		rit_->atom(Residue::AL_N).index());
	return b;
      } else return *ait_;
    }

    std::vector<Residue>::const_iterator rit_, rend_;
    Residue::Bonds_iterator ait_, aend_;
    bool nc_;
  };
};


#endif
