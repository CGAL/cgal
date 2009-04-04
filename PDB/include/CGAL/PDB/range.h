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
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range.hpp>
#include <algorithm>


namespace CGAL { namespace PDB {

template <class Range>
unsigned int size(Range r) {
  return std::distance(r.begin(), r.end());
}

template <class Range, class F>
const F& for_each(const Range& r, const F &f) {
  CGAL_PDB_FOREACH(typename Range::iterator::reference v, r) {
    f(v);
  }
  return f;
}
template <class Range, class F>
const F& for_each(Range& r, const F &f) {
  CGAL_PDB_FOREACH(typename Range::iterator::reference v, r) {
    f(v);
  }
  return f;
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
      typedef const ucname& result_type;                       \
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
  typedef unsigned int result_type;
  result_type operator()(A a) const {
    return a.index();
  }
};

template <class A>
struct Get_key {
  typedef typename A::Key result_type;
  result_type operator()(A a) const {
    return a.key();
  }
};

#define CGAL_PDB_MAKE_RANGE(ucname, lcname)             \
    template <class Range>                                              \
    boost::iterator_range<boost::transform_iterator< Get_##lcname<typename Range::iterator::reference> , typename Range::iterator> > \
    make_##lcname##_range(const Range& r){                              \
      typedef boost::transform_iterator<Get_##lcname<typename Range::iterator::reference>, typename Range::iterator> Tr; \
      return boost::make_iterator_range(Tr(r.begin(), Get_##lcname<typename Range::iterator::reference>()), \
                                        Tr(r.end(), Get_##lcname<typename Range::iterator::reference>())); \
    }                                                                   \
    template <class Range>                                              \
    boost::iterator_range<boost::transform_iterator< Get_##lcname<typename Range::iterator::reference> , typename Range::iterator> > \
    make_##lcname##_range(Range& r){                                    \
      typedef boost::transform_iterator<Get_##lcname<typename Range::iterator::reference>, typename Range::iterator> Tr; \
      return boost::make_iterator_range(Tr(r.begin(), Get_##lcname<typename Range::iterator::reference>()), \
                                        Tr(r.end(), Get_##lcname<typename Range::iterator::reference>())); \
    }                                                                   \
    
    
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
boost::iterator_range<boost::transform_iterator<Get_bond_indices, typename Range::iterator> >
make_bond_indices_range(const Range &r){
  typedef boost::transform_iterator<Get_bond_indices, typename Range::iterator> Tr;
  return boost::make_iterator_range(Tr(r.begin(),Get_bond_indices()),
                                    Tr(r.end(), Get_bond_indices()));
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
boost::iterator_range<boost::filter_iterator<Is_backbone, typename Range::iterator> >
make_backbone_range(const Range& r){
  return boost::make_iterator_range(boost::make_filter_iterator(Is_backbone(), 
                                                                r.begin(), r.end()),
                                    boost::make_filter_iterator(Is_backbone(), 
                                                                r.end(), r.end()));                                 
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
boost::iterator_range<boost::filter_iterator<Is_CA, typename Range::iterator> >
make_ca_range(const Range& r){
  return boost::make_iterator_range(boost::make_filter_iterator(Is_CA(), 
                                                                r.begin(),
                                                                r.end()),
                                    boost::make_filter_iterator(Is_CA(), 
                                                                r.end(),
                                                                r.end()));
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
boost::iterator_range<boost::filter_iterator<Is_ok_bond<OKA> , typename Range::iterator> >
make_ok_bond_range( OKA oka, const Range& r){
  return boost::make_iterator_range(boost::make_filter_iterator(Is_ok_bond<OKA>(oka), 
                                                                r.begin(), r.end()),
                                    boost::make_filter_iterator(Is_ok_bond<OKA>(oka), 
                                                                r.end(), r.end()));
}



/*! \example extracting_geometry.cpp
  This example shows how to use the various iterator adaptors for extracting geometry and connectivity.
 */

}}
#endif
