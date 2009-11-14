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

#ifndef CGAL_DSR_PDB_CHAIN_H
#define CGAL_DSR_PDB_CHAIN_H

#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Atom.h>
#include <CGAL/PDB/Monomer.h>
#include <CGAL/PDB/small_map.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>
#include <CGAL/PDB/internal/dummies.h>
#include <CGAL/PDB/internal/Nested_iterator.h>
#include <boost/tuple/tuple.hpp>

namespace CGAL { namespace PDB {

/*!
  \file Chain.h 

  This file contains the class Chain.
*/

class Model;
  
//! A class representing a single chain of a protein.  
/*!
    
 */
class Chain {
  friend class Model;

public:

  //! The type for storing residue indices in the PDB
  typedef CGAL::Label<Chain> Monomer_key;

  struct IR_tag{};
  typedef CGAL::Label<IR_tag> IR_key;

  //! The value_type returned by the Monomer_iterators
  class Monomer_pair: public small_map_value_type<Monomer_key, Monomer> {
    typedef small_map_value_type<Monomer_key, Monomer> P;
  public:
    Monomer_pair(Monomer_key k, Monomer a): P(k,a){}
    Monomer_pair(Monomer_key k): P(k){}
    Monomer_pair(){}
    const Monomer &monomer() const {return P::data();}
    Monomer &monomer() {return P::data();}
  };

  typedef small_map<Monomer_pair> Container;
  
  //! Default
  Chain();

  CGAL_ITERATOR(Monomer, monomer, Container::const_iterator,
                Container::iterator,
                residues_.begin(),
		residues_.end());


  //! A unique identified of an atom in the chain
  class Atom_key {
    Monomer_key mk_;
    Monomer::Atom_key ak_;
   public:
    Atom_key(Monomer_key mk, Monomer::Atom_key ak): mk_(mk), ak_(ak){}
    Atom_key(){}
    operator Monomer::Atom_key() const {
      return ak_;
    }
    operator Monomer_key() const {
      return mk_;
    }
    Monomer::Atom_key atom_key() const {
      return ak_;
    }
    Monomer_key monomer_key() const {
      return mk_;
    }
  };


  //! The endpoint of a bond
  class Bond_endpoint {
    const Atom *atom_;
    Atom_key ak_;
    friend class Bond_it;
  public:
    Bond_endpoint(Atom_key ak, const Atom *atom): atom_(atom), ak_(ak){}    
    Bond_endpoint():atom_(NULL){}
    const Atom &atom() const {return *atom_;}
    Atom_key key() const {return ak_;}
  };

  //! A chemical bond within the protein
  typedef std::pair<Bond_endpoint, Bond_endpoint> Bond; 

  CGAL_INSERT(Monomer, insert_internal(k,m));

  class Atom_pair: public small_map_value_type<Atom_key, Atom> {
      Atom_key index_;
      Atom *atom_;
    public:
      Atom_pair(Atom_key f, Atom* s): index_(f), atom_(s){}
      Atom_key key() const {return index_;}
      Atom &atom() const {return *atom_;}
      Atom_pair():atom_(NULL){}
    };

  class Atom_const_pair: public small_map_value_type<Atom_key, Atom> {
    Atom_key index_;
    const Atom *atom_;
  public:
    Atom_const_pair(Atom_key f, const Atom* s): index_(f), atom_(s){}
    Atom_key key() const {return index_;}
    const Atom &atom() const {return *atom_;}
    Atom_const_pair():atom_(NULL){}
  };
private:

  //! \cond
  Monomers::iterator insert_internal(Monomer_key k, const Monomer &m);

  struct Iterator_traits {
    typedef Monomers Outer;
    typedef Monomer::Atoms Inner;
    typedef Atom_pair value_type;
    struct Get_inner{
      Inner operator()(Outer::iterator it) const {
	return it->monomer().atoms();
      }
    };
    struct Make_value{
      value_type operator()(Outer::iterator oit,
                            Inner::iterator iit) const {
	return value_type(Atom_key(oit->key(), iit->key()), &iit->atom());
      }
    };
  };


  struct Iterator_const_traits {
    typedef Monomer_consts Outer;
    typedef Monomer::Atom_consts Inner;
    typedef Atom_const_pair value_type;
    struct Get_inner{
      Inner operator()(Outer::iterator it) const {
	return it->monomer().atoms();
      }
    };
    struct Make_value{
      value_type operator()(Outer::iterator oit,
                            Inner::iterator iit) const {
	return value_type(Atom_key(oit->key(), iit->key()), &iit->atom());
      }
    };
  };
  //! \endcond
public:
  
  
  //! An iterator to iterate through all the atoms of the protein  
  CGAL_ITERATOR(Atom, atom,
                internal::Nested_iterator<Iterator_const_traits >,
                internal::Nested_iterator<Iterator_traits >,
                boost::make_iterator_range(residues_.begin(), residues_.end()),
                boost::make_iterator_range(residues_.end(), residues_.end()));

  
  //! This is non-const time.
  unsigned int number_of_atoms() const;

  //! \cond
  void swap_with(Chain &o);

  class Bond_it {
    friend class Chain;
  public:
    typedef Bond value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t difference_type;
    typedef const Bond& reference;
    typedef const Bond* pointer;
    
    reference operator*() const {
      return ret_;
    }
    pointer operator->() const {
      return &ret_;
    }
    Bond_it operator++() {
      if (between_) {
	between_=false;
	++rit_;
	if (rit_ != rend_) {
	  ait_= rit_->monomer().bonds().begin();
	  make_bond();
	}
      } else {
	++ait_;
	if (ait_== rit_->monomer().bonds().end()) {
          Monomer_consts::iterator nrit= rit_;
	  ++nrit;
	  between_=true;
	  if (nrit == rend_ 
	      || nrit->key().index() != rit_->key().index()+1){
	    operator++();
	  } else {
	    make_bond();
	  }
	} else {
	  make_bond();
	}
      }
      return *this;
    }
    bool operator==(const Bond_it& o) const {
      if (between_ != o.between_) return false;
      if (rit_ == rend_) return rit_==o.rit_;
      else return rit_== o.rit_ && ait_ == o.ait_;
    }
    bool operator!=(const Bond_it& o) const {
      return !operator==(o);
    }
    
    Bond_it(){}

    CGAL_COPY_CONSTRUCTOR(Bond_it);

    Bond_it(Monomer_consts r): rit_(r.begin()), rend_(r.end()), 
				       between_(false){
      if (!r.empty()) {
	ait_= rit_->monomer().bonds().begin();
	make_bond();
      }
    }
  private:
    void copy_from(const Bond_it &o) {
      rit_= o.rit_;
      rend_= o.rend_;
      if (rit_!= rend_) {
	ait_= o.ait_;
	ret_= o.ret_;
      }
      between_= o.between_;
    }

    void make_bond() {
      if (between_) {
	CGAL_assertion(rit_+1 != rend_);
	ret_= Bond(Bond_endpoint(Atom_key(rit_->key(), rit_->monomer().back_atom()->key()),
				 &rit_->monomer().back_atom()->atom()),
		   Bond_endpoint(Atom_key((rit_+1)->key(), (rit_+1)->monomer().front_atom()->key()),
				 &(rit_+1)->monomer().front_atom()->atom()));
      } else {
	ret_= Bond(Bond_endpoint(Atom_key(rit_->key(), ait_->first.key()),
				 &ait_->first.atom()),
		   Bond_endpoint(Atom_key(rit_->key(), ait_->second.key()),
				 &ait_->second.atom()));
      } 
    }

    Monomer_consts::iterator rit_, rend_;
    Monomer::Bonds::const_iterator ait_;
    bool between_;
    Bond ret_;
  };
  //! \endcond

  CGAL_CONST_ITERATOR(Bond, bond, Bond_it,
                      boost::make_iterator_range(residues_.begin(), residues_.end()),
                      boost::make_iterator_range(residues_.end(), residues_.end()));
 

  //! This is non-const time.
  unsigned int number_of_bonds() const;
  
  /* //! Return the spherical_coordinates of the atom relative to its parent.

     Spherical_point spherical_coordinates(Atom::Index atom_index) const;

     //! Return the parent of an atom with regards to computing spherical coordinates.
     Atom::Index parent_atom(Atom::Index atom_index) const;*/


  //! The sequence of residue types.
  std::vector<Monomer::Type> sequence() const;

  //! Write as part of pdb file.
  int  write(char chain, int start_index, std::ostream &out) const;

  //! Write a pdb file. 
  /*!
    See check_protein.cpp for an example of using this to write a
    pdb file.
  */
  void write_pdb(std::ostream &out) const;

  //! Dump as human readable.
  void dump(std::ostream &out) const;

  //! Dump as human readable.
  std::ostream& write(std::ostream &out) const{dump(out); return out;}

  CGAL_FIND(Monomer, residues_.find(k), residues_.end());

  /*const Atom &atom(unsigned int i) const {
    for (Const_residues_iterator it = residues_begin(); it != residues_.end(); ++it){
	
    }
    }*/


  //! Return whether bonds have been computed for this protein. 
  bool has_bonds() const;

  //! Set whether the protein has bonds or not.
  void set_has_bonds(bool tf);

  CGAL_GETSET(std::string, name, name_);
private:
    
  //unsigned int residue_offset_of_atom_key(Atom::Index i) const;

  //unsigned int residue_offset(Residue::Index i) const;
    
  Container residues_;
  std::vector<std::string> header_;
  //static Residue dummy_residue_;
  typedef small_map<small_map_value_type<Monomer_key,
					 small_map<small_map_value_type<IR_key, Monomer> > > > IR_Map;
  IR_Map insert_residues_;
  std::string name_;
};


CGAL_SWAP(Chain)
CGAL_OUTPUT(Chain)




//! Assign unique indices to all atoms in the Chain, starting at optional start value
/*!
  This returns the next unused index. 
*/
inline int index_atoms(const Chain &c, int start=0) {
  CGAL_PDB_FOREACH(std::iterator_traits<Chain::Monomer_consts::iterator>::reference m, c.monomers()) {
    start= index_atoms(m.monomer(), start);
  }
  return start;
}
}}

#endif
