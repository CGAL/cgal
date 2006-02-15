// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp <fransw@inf.ethz.ch>
//                 Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_QP_SOLVER_ITERATOR_H
#define CGAL_QP_SOLVER_ITERATOR_H

#ifndef CGAL_ITERATOR
#  include <CGAL/iterator.h>
#endif

CGAL_BEGIN_NAMESPACE

template < typename T >
class Const_oneset_iterator {
public:
  
  // types
  typedef  std::random_access_iterator_tag    iterator_category;
  typedef  ptrdiff_t                          difference_type;
  typedef  T                                  value_type;
  typedef  value_type*                        pointer;
  typedef  value_type&                        reference;
  
  typedef  Const_oneset_iterator<T>           Self;
  typedef  difference_type                    Diff;
  typedef  value_type                         Val;
  typedef  pointer                            Ptr;
  typedef  reference                          Ref;
  
  // construction
  Const_oneset_iterator( const T& t = T(), Diff n = 0)
    : value( t), index( n)
  { }
  
  // access
  Ref        operator *  ( )       { return  value; }
  const Ref  operator *  ( ) const { return  value; }
  Ptr        operator -> ( )       { return &value; }
  const Ptr  operator -> ( ) const { return &value; }
  
  // equality operator
  bool       operator == ( const Self& x) const { return ( index==x.index); }
  bool       operator != ( const Self& x) const { return ( index!=x.index); }
  
  // forward operations
  // ------------------
  Self&      operator ++ (    ) {                   ++index; return *this; }
  Self       operator ++ ( int) { Self tmp = *this; ++index; return tmp;   }
  
  // bidirectional operations
  // ------------------------
  Self&      operator -- (    ) {                   --index; return *this; }
  Self       operator -- ( int) { Self tmp = *this; --index; return tmp;   }
  
  // random access operations
  // ------------------------
  // access
  Ref        operator [] ( Diff i)       { return value;}
  const Ref  operator [] ( Diff i) const { return value;}
  
  // less operator
  bool       operator <  ( const Self& x) const { return ( index < x.index);}
  
  // arithmetic operations
  Self&      operator += ( Diff n) { index += n; return *this; }
  Self&      operator -= ( Diff n) { index -= n; return *this; }
  
  Self       operator +  ( Diff n) const { Self tmp = *this; return tmp+=n; }
  Self       operator -  ( Diff n) const { Self tmp = *this; return tmp-=n; }
  
  Diff       operator -  ( const Self& x) const { return index - x.index; }
  
private:
  
  // data members
  Val   value;
  Diff  index;
};

CGAL_END_NAMESPACE

#endif // CGAL_QP_SOLVER_ITERATOR_H

// ===== EOF ==================================================================
