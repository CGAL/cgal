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

#ifndef CGAL_DSR_PDB_ITERATOR_H
#define CGAL_DSR_PDB_ITERATOR_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Chain.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>

CGAL_PDB_BEGIN_NAMESPACE

struct Get_atom {
  typedef Atom result_type;
  template <class Vt>
  const Atom &operator()(const Vt &v) const {
    return v.atom();
  }
  template <class Vt>
  Atom &operator()(Vt &v) const {
    return v.atom();
  }
};


//! Return an interator which returns an Atom object
template <class It>
boost::transform_iterator<Get_atom, It> make_atom_iterator(It it){
  return boost::transform_iterator<Get_atom, It>(it, Get_atom());
}


struct Get_bond_indices {
  typedef std::pair<unsigned int,unsigned int> result_type;
  template <class B>
  const result_type& operator()(const B & b) const {
    static result_type ret;
    //Chain::Monomer_iterator ma= b.first.first;
    //Monomer::Atom_label la= b.first.second;

    ret= std::make_pair(b.first.atom().index().index(),
			b.second.atom().index().index());
    return ret;
  }
    
};


//! Return an iterator which returns a pair of indices for a bond
template <class It>
boost::transform_iterator<Get_bond_indices, It>
make_bond_indices_iterator(It it){
  return boost::transform_iterator<Get_bond_indices, It>(it,
							 Get_bond_indices());
}


//! Return true if an atom is a backbone atoms.
struct Is_backbone {
  typedef bool result_type;
  template <class A>
  bool operator()(const A &a) {
    return static_cast<Monomer::Atom_key>(a.key()) == Monomer::AL_C 
      || static_cast<Monomer::Atom_key>(a.key()) == Monomer::AL_CA 
      || static_cast<Monomer::Atom_key>(a.key()) == Monomer::AL_N;
  }
};


//! Return an iterator which returns skips non-backbone atoms
template <class It>
boost::filter_iterator<Is_backbone, It>
make_backbone_iterator(It itb, It ite){
  return boost::make_filter_iterator(Is_backbone(), 
				     itb, ite);
}

//! Return true if an atom is a AC.
struct Is_CA {
  typedef bool result_type;
  template <class A>
  bool operator()(const A &a) {
    return static_cast<Monomer::Atom_key>(a.key()) == Monomer::AL_CA;
  }
};

//! Return an iterator which returns skips non-backbone atoms
template <class It>
boost::filter_iterator<Is_CA, It>
make_ca_iterator(It itb, It ite){
  return boost::make_filter_iterator(Is_CA(), 
				     itb, ite);
}

template <class OK_atom>
struct Is_ok_bond {
  typedef bool result_type;
  Is_ok_bond(OK_atom ao): ao_(ao){}

  struct Atom_pair{
    Monomer::Atom_key k_;
    const Atom *a_;
    Atom_pair(Monomer::Atom_key k, const Atom *a): k_(k), a_(a){}
    const Atom& atom() const {return a_;}
    Monomer::Atom_key key() const {return k_;}
  };

  template <class B>
  bool operator()(const B &b) {
    Monomer::Atom_key ka= static_cast<Monomer::Atom_key>(b.first.key());
    Monomer::Atom_key kb= static_cast<Monomer::Atom_key>(b.second.key());
    const Atom *pa =&b.first.atom();
    const Atom *pb =&b.second.atom();
    if ( ao_(Atom_pair(ka, pa)) && ao_(Atom_pair(kb, pb))) {
      return true;
    } else {
      /*std::cout << "Rejected " << ka << "(" << pa << ") to " << kb 
	<< "(" << pb << ")" << std::endl;*/
      return false;
    }
  }
  OK_atom ao_;
};

//! Return an iterator which returns skips non-backbone atoms
template <class It, class OKA>
boost::filter_iterator<Is_ok_bond<OKA> , It>
make_ok_bond_iterator( OKA oka, It itb, It ite){
  return boost::make_filter_iterator(Is_ok_bond<OKA>(oka), 
				     itb, ite);
}



/*! \example extracting_geometry.cpp
  This example shows how to use the various iterator adaptors for extracting geometry and connectivity.
 */

CGAL_PDB_END_NAMESPACE
#endif
