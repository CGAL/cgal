// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/iterator_tools.h
// package       : Nef_2
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// ============================================================================
#ifndef CGAL_ITERATORTOOLS_H
#define CGAL_ITERATORTOOLS_H

#include <CGAL/basic.h>
#include <CGAL/circulator.h>

CGAL_BEGIN_NAMESPACE

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

    bool operator==( CGAL_NULL_TYPE p ) const {
      CGAL_assertion( p == NULL );
      return Iter::operator==( Iter(NULL) );
    }
    bool operator!=( CGAL_NULL_TYPE p ) const {
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

CGAL_END_NAMESPACE
#endif // CGAL_ITERATORTOOLS_H
