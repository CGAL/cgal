// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_CIRCULATOR_PROJECT_H
#define CGAL_CIRCULATOR_PROJECT_H 1

#include <CGAL/circulator.h>

namespace CGAL {

template < class C,
           class Fct,
           class Ref = typename C::reference,
           class Ptr = typename C::pointer>
class Circulator_project {
protected:
  C        nt;    // The internal circulator.
public:
  typedef C  Circulator;
  typedef Circulator_project<C,Fct,Ref,Ptr> Self;

  typedef typename  C::iterator_category   iterator_category;
  typedef typename  Fct::result_type       value_type;
  typedef typename  C::difference_type     difference_type;
  typedef typename  C::size_type           size_type;
  typedef           Ref                    reference;
  typedef           Ptr                    pointer;

  // CREATION
  // --------

  Circulator_project() {}
  Circulator_project( Circulator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Circulator  current_circulator() const { return nt;}
  Ptr ptr() const {
    Fct fct;
    return &(fct(*nt));
  }

  bool operator==( std::nullptr_t CGAL_assertion_code(p) ) const {
    CGAL_assertion( p == 0);
    return ( nt == 0);
  }
  bool  operator!=( std::nullptr_t p) const { return !(*this == p); }
  bool  operator==( const Self& i) const { return ( nt == i.nt); }
  bool  operator!=( const Self& i) const { return !(*this == i); }
  Ref   operator*()  const { return *ptr(); }
  Ptr   operator->() const { return ptr(); }
  Self& operator++() {
    ++nt;
    return *this;
  }
  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  // OPERATIONS Bidirectional Category
  // ---------------------------------

  Self& operator--() {
    --nt;
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

  // OPERATIONS Random Access Category
  // ---------------------------------

  Self  min_circulator() const {
    return Self( nt.min_circulator());
  }
  Self& operator+=( difference_type n) {
    nt += n;
    return *this;
  }
  Self  operator+( difference_type n) const {
    Self tmp = *this;
    return tmp += n;
  }
  Self& operator-=( difference_type n) {
    return operator+=( -n);
  }
  Self  operator-( difference_type n) const {
    Self tmp = *this;
    return tmp += -n;
  }
  difference_type  operator-( const Self& i) const {
    return nt - i.nt;
  }
  Ref  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
  }
};

template < class Dist, class Fct, class C, class Ref, class Ptr>
inline
Circulator_project<C,Fct,Ref,Ptr>
operator+( Dist n, Circulator_project<C,Fct,Ref,Ptr> i) { return i += n; }

} //namespace CGAL
#endif // CGAL_CIRCULATOR_PROJECT_H //
// EOF //
