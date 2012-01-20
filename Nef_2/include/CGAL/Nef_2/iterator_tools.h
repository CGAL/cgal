// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_ITERATORTOOLS_H
#define CGAL_ITERATORTOOLS_H

#include <CGAL/basic.h>
#include <CGAL/circulator.h>

namespace CGAL {

template <typename Iter, typename Move> 
class CircFromIt : public Iter {
    // Ptr  node;    // The internal node ptr inherited from It.
    typedef CircFromIt<Iter,Move> Self;
public:
    typedef typename Iter::iterator_category Icategory;
    typedef I_Circulator_from_iterator_traits<Icategory> CTraits;
    typedef typename CTraits::iterator_category iterator_category;

    CircFromIt() : Iter(0) {}
    CircFromIt(Iter i) : Iter(i) {}

// OPERATIONS Forward Category
// ---------------------------

    bool operator==( Nullptr_t CGAL_assertion_code(p) ) const {
      CGAL_assertion( p == NULL );
      return Iter::operator==( Iter(NULL) );
    }
    bool operator!=( Nullptr_t p ) const {
      return !(*this == p);
    }
    bool operator==( const Self& i ) const {
      return Iter::operator==(i);
    }
    bool operator!=( const Self& i) const {
        return !(*this == i);
    }

    Self& operator++() {
      Move move;
      move.forward(*this);
      return *this;
    }
    Self  operator++(int) {
      CircFromIt tmp = *this;
      ++*this;
      return tmp;
    }

// OPERATIONS Bidirectional Category
// ---------------------------------

    Self& operator--() {
      Move move;
      move.backward(*this);
      return *this;
    }
    Self  operator--(int) {
      CircFromIt tmp = *this;
      --*this;
      return tmp;
    }

};

template <typename Iter, typename Pnt> 
class PntItFromVertIt : public Iter {
public:
  typedef PntItFromVertIt<Iter,Pnt> Self;
  typedef Iter Base;
  typedef Pnt  value_type;
  typedef const Pnt* pointer;
  typedef const Pnt& reference;

  PntItFromVertIt() : Base() {}
  PntItFromVertIt(Iter it) : Base(it) {}
  PntItFromVertIt(const Self& it) : Base(it) {}

  reference operator*() const 
  { return Base::operator*().point(); }
  pointer operator->() const 
  { return &(operator*()); }
  Self& operator++() { return (Self&)Base::operator++(); }
  Self operator++(int) { Self tmp=*this; ++*this; return tmp; }

};

template <class H>
std::string PH(H h)
{ if (h == H()) return "nil"; return h->debug(); }

} //namespace CGAL
#endif // CGAL_ITERATORTOOLS_H
