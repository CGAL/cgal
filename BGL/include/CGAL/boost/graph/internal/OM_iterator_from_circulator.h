// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Philipp Moeller

#ifndef CGAL_OM_ITERATOR_FROM_CICULATOR_H
#define CGAL_OM_ITERATOR_FROM_CICULATOR_H

#include <iostream>
#include <cstddef>

#include <CGAL/circulator.h>

#include <boost/operators.hpp>
#include <boost/concept_check.hpp>

#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

namespace CGAL {

// adapted from circulator.h, does not support
// random_access_circulators and returns the underlying circulator
// instead of dereferencing it
template <class C, bool Prevent_deref = true>
class OM_iterator_from_circulator {
private:
  // The m_anchor is normalized to be a minimal circulator.
  C         m_anchor;
  C         current;
  bool past_the_end;

  typedef  std::iterator_traits<C>                       I_traits;
  typedef  typename  I_traits::iterator_category         I_Iter_cat;
  typedef  I_Iterator_from_circulator_traits<I_Iter_cat> I__traits;

public:
  typedef C                                               Circulator;
  typedef OM_iterator_from_circulator<C, Prevent_deref> Self;

  typedef typename I__traits::iterator_category iterator_category;

  typedef typename 
  boost::mpl::if_c<  Prevent_deref
                   , C
                   , typename C::value_type
                  >::type             value_type;

  typedef typename C::difference_type difference_type;
  typedef typename 
  boost::mpl::if_c<  Prevent_deref
                   , C&
                   , typename C::reference
                  >::type             reference;
  typedef typename 
  boost::mpl::if_c<  Prevent_deref
                   , C*
                   , typename C::reference
                  >::type             pointer;

  OM_iterator_from_circulator(){}

  OM_iterator_from_circulator(const C circ, int n)
    : m_anchor(circ), current(circ), past_the_end(n == 1) {}


  bool done() const
  {
    return past_the_end || (! current);
  }


  bool operator==( const Self& i) const {
    CGAL_assertion( m_anchor == i.m_anchor);  // same anchor?
    return (done() && i.done()) || (((!done()) && (!i.done())) && ( current == i.current));
  }
  
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }


// we cannot enable_if on operator* and operator-> because they do not
// depend on template parameters directly and default template
// arguments for functions are C++11. we redirect on helper member
// templates as a work-around.
private:
  template <bool Prevent_deref_>
  typename boost::enable_if_c<Prevent_deref_, reference>::type
  indirection() const {
    return const_cast<Self*>(this)->current;
  }
  template <bool Prevent_deref_>
  typename boost::disable_if_c<Prevent_deref_, reference>::type
  indirection() const {
    return *current;
  }
public:
  reference operator*() const {
    return indirection<Prevent_deref>();
  }
  
private:
  template <bool Prevent_deref_>
  typename boost::disable_if_c<Prevent_deref_, pointer>::type
  structure_dereference() {
    return &(*current);
  }
  template <bool Prevent_deref_>
  typename boost::enable_if_c<Prevent_deref_, pointer>::type
  structure_dereference() {
    return &current;
  }
public:
  pointer operator->() const {
    return structure_dereference<Prevent_deref>();
  }

  Self& operator++() {
    ++current;
    return *this;
  }
  Self  operator++(int) {
    std::cerr << "operator++(int)" << std::endl;
    Self tmp = *this;
    ++*this;
    return tmp;
  }
  Self& operator--() {
    --current;
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }
  /*
  Self& operator+=( difference_type n) {
    if ( n < 0 && current == m_anchor)  // We are leaving the anchor.
      --m_winding;
    current += n;
    if ( n > 0 && current == m_anchor)  // Back again at the anchor.
      ++m_winding;
    return *this;
  }

  bool operator<( const Self& i) const {
    CGAL_assertion( m_anchor  != NULL);
    CGAL_assertion( m_anchor  == i.m_anchor);
    return (     (m_winding < i.m_winding)
                 || (    (m_winding == i.m_winding)
                         && (current - m_anchor) < (i.current - m_anchor)
                   )
      );
  }
  bool operator> ( const Self& i) const { return i < *this; }
  bool operator<=( const Self& i) const { return !(i < *this); }
  bool operator>=( const Self& i) const { return !(*this < i); }
  */
  const C*    anchor()             const { return m_anchor;}

  Circulator  current_circulator() const { return current;}
};

} // CGAL

#endif /* CGAL_OM_ITERATOR_FROM_CICULATOR_H */
