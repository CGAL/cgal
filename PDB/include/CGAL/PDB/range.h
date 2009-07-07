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

#ifndef CGAL_DSR_PDB_RANGE_H
#define CGAL_DSR_PDB_RANGE_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/Monomer.h>
#include <CGAL/PDB/Heterogen.h>
#include <CGAL/PDB/PDB.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range.hpp>
#include <CGAL/PDB/internal/rangelib/filter.hpp>
#include <CGAL/PDB/internal/rangelib/transform.hpp>
#include <CGAL/PDB/internal/rangelib/algo.hpp>
#include <algorithm>


namespace CGAL { namespace PDB {

    using internal::rangelib::rng::for_each;
    using internal::rangelib::rng::copy;

//! boost changes from size to distance a 1.34, so use this
template <class Range>
unsigned int distance(const Range &r) {
  return std::distance(r.begin(), r.end());
}
    
#define CGAL_PDB_GET_MAP(ucname, lcname)                \
    template <class From>                               \
    struct Get_##lcname {                               \
      typedef ucname& result_type;                      \
      result_type operator()(From v) const {            \
        return v.lcname();                              \
      }                                                 \
    };                                                  \
    template <class From2>                               \
    struct Get_##lcname<const From2> {                   \
      typedef const ucname& result_type;                      \
      result_type operator()(From2 v) const {                 \
        return v.lcname();                              \
      }                                                 \
    };                                                  \

CGAL_PDB_GET_MAP(Atom, atom);
CGAL_PDB_GET_MAP(Monomer, monomer);
CGAL_PDB_GET_MAP(Chain, chain);
CGAL_PDB_GET_MAP(Model, model);
CGAL_PDB_GET_MAP(Heterogen, heterogen);

template <class A>
struct Get_point {
  typedef const Point & result_type;
  result_type operator()(A a) const {
    return a.point();
  }
};

template <class A>
struct Get_index {
  typedef Atom::Index result_type;
  result_type operator()(A a) const {
    return a.index();
  }
};

template <class A>
struct Get_key {
  typedef typename A::Key result_type;
  result_type operator()(const A& a) const {
    return a.key();
  }
};

#define CGAL_PDB_MAKE_RANGE(ucname, lcname)             \
    template <class Range>                                              \
    internal::rangelib::transformed_range<Range, Get_##lcname<typename std::iterator_traits<typename Range::iterator>::reference> > \
    make_##lcname##_range( Range r){                                   \
      return internal::rangelib::transformed(r, Get_##lcname<typename std::iterator_traits<typename Range::iterator>::reference>()); \
    }                                                                   \

    
    /**
    template <class Range>                                              \
    internal::rangelib::transformed_range<Range, Get_##lcname<const typename Range::iterator::reference> > \
    make_##lcname##_range(const Range& r){                                    \
      return internal::rangelib::transformed(r, Get_##lcname<const typename Range::iterator::reference>()); \
    }                                                                   \

     */
CGAL_PDB_MAKE_RANGE(Atom, atom);
CGAL_PDB_MAKE_RANGE(Monomer, monomer);
CGAL_PDB_MAKE_RANGE(Chain, chain);
CGAL_PDB_MAKE_RANGE(Model, model);
CGAL_PDB_MAKE_RANGE(Heterogen, heterogen);
CGAL_PDB_MAKE_RANGE(Point, point);
CGAL_PDB_MAKE_RANGE(Index, index);
CGAL_PDB_MAKE_RANGE(Key, key);


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


//! Return an iterator range which returns a pair of indices for a bond
template <class Range>
internal::rangelib::transformed_range<Range, Get_bond_indices>
make_bond_indices_range(Range r){
  return internal::rangelib::transformed(r,Get_bond_indices());
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


//! Return an iterator range which returns skips non-backbone atoms
template <class Range>
internal::rangelib::filtered_range<Range, Is_backbone>
make_backbone_range(Range r){
  return internal::rangelib::filtered(r, Is_backbone());
}

//! Return true if an atom is a AC.
struct Is_CA {
  typedef bool result_type;
  template <class A>
  bool operator()(const A &a) {
    return static_cast<Monomer::Atom_key>(a.key()) == Monomer::AL_CA;
  }
};

//! Return an iterator range which returns skips non-backbone atoms
template <class Range>
internal::rangelib::filtered_range<Range, Is_CA>
make_ca_range(Range r){
  return internal::rangelib::filtered(r, Is_CA());
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

//! Return an iterator range which returns skips non-backbone atoms
    template <class Range, class OKA>
internal::rangelib::filtered_range<Range, Is_ok_bond<OKA> >
make_filtered_bond_range(OKA oka, Range r){
  return internal::rangelib::filtered(r, Is_ok_bond<OKA>(oka));
}



/*! \example extracting_geometry.cpp
  This example shows how to use the various iterator adaptors for extracting geometry and connectivity.
 */

}}
#endif
