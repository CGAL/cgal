// Copyright (c) 2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_ITERATOR_IDENTITY_H
#define CGAL_ITERATOR_IDENTITY_H 1

#include <CGAL/circulator.h>

namespace CGAL {

template < class I,
           class Ref  = typename std::iterator_traits<I>::reference,
           class Ptr  = typename std::iterator_traits<I>::pointer,
           class Val  = typename std::iterator_traits<I>::value_type,
           class Dist = typename std::iterator_traits<I>::difference_type,
           class Ctg  = typename std::iterator_traits<I>::iterator_category>
class Iterator_identity {
protected:
  I        nt;    // The internal iterator.
public:
  typedef I      Iterator;
  typedef Iterator_identity<I,Ref,Ptr,Val,Dist,Ctg>   Self;
  typedef Ctg    iterator_category;
  typedef Val    value_type;
  typedef Dist   difference_type;
  typedef Ref    reference;
  typedef Ptr    pointer;

  // CREATION
  // --------

  Iterator_identity() {}
  Iterator_identity( Iterator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------

  Iterator  current_iterator() const { return nt;}

  bool operator==( const Self& i) const {
    return ( nt == i.nt);                                    //###//
  }
  bool operator!=( const Self& i) const {
    return !(*this == i);
  }
  Ref  operator*() const {
    return *nt;                                              //###//
  }
  Ptr  operator->() const {
    return nt.operator->();                                  //###//
  }
  Self& operator++() {
    ++nt;                                                    //###//
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
    --nt;                                                    //###//
    return *this;
  }
  Self  operator--(int) {
    Self tmp = *this;
    --*this;
    return tmp;
  }

  // OPERATIONS Random Access Category
  // ---------------------------------

  Self& operator+=( difference_type n) {
    nt += n;                                                 //###//
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
    return nt - i.nt;                                        //###//
  }
  Ref  operator[]( difference_type n) const {
    Self tmp = *this;
    tmp += n;
    return tmp.operator*();
  }
  bool operator<( const Self& i) const {
    return ( nt < i.nt);                                     //###//
  }
  bool operator>( const Self& i) const {
    return i < *this;
  }
  bool operator<=( const Self& i) const {
    return !(i < *this);
  }
  bool operator>=( const Self& i) const {
    return !(*this < i);
  }
};

template < class I, class Ref, class Ptr, class Val,
           class Dist, class Ctg>
inline
Iterator_identity<I,Ref,Ptr,Val,Dist,Ctg>
operator+( Dist n, Iterator_identity<I,Ref,Ptr,Val,Dist,Ctg> i)
{ return i += n; }

} //namespace CGAL
#endif // CGAL_ITERATOR_IDENTITY_H //
// EOF //
