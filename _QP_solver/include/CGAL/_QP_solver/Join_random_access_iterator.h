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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
#define CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H

#include <CGAL/basic.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class RndAccIt, class Operation >
class Join_random_access_iterator_1 {
  public:
    typedef  typename Operation::result_type  value_type;
    typedef  ptrdiff_t                        difference_type;
    typedef  value_type*                      pointer;
    typedef  value_type&                      reference;
    typedef  std::random_access_iterator_tag  iterator_category;

    typedef  Join_random_access_iterator_1<RndAccIt,Operation>  Self;
    typedef  value_type                                         Val;
    typedef  difference_type                                    Dist;
    typedef  reference                                          Ref;

    // forward operations
    Join_random_access_iterator_1( ) { }
    Join_random_access_iterator_1( const RndAccIt& it_1) : it1( it_1) { }
    Join_random_access_iterator_1( const RndAccIt&  it_1,
				   const Operation& operation)
	: it1( it_1), op( operation) { }

    bool   operator == ( const Self& it) const { return ( it1 == it.it1); }
    bool   operator != ( const Self& it) const { return ( it1 != it.it1); }

    Val    operator *  ( ) const { return op( *it1); }

    Self&  operator ++ (    ) {                   ++it1; return *this; }
    Self   operator ++ ( int) { Self tmp = *this; ++it1; return tmp;   }

    // bidirectional operations
    Self&  operator -- (    ) {                   --it1; return *this; }
    Self   operator -- ( int) { Self tmp = *this; --it1; return tmp;   }

    // random access operations
    Self&  operator += ( Dist n) { it1 += n; return *this; }
    Self&  operator -= ( Dist n) { it1 -= n; return *this; }

    Self   operator +  ( Dist n) const { Self tmp = *this; return tmp += n; }
    Self   operator -  ( Dist n) const { Self tmp = *this; return tmp -= n; }

    Dist   operator -  ( const Self& it) const { return it1 - it.it1; }

    Val    operator [] ( int i) const { return op( it1[ i]); }

    bool   operator <  ( const Self& it) const { return ( it1 <  it.it1); }
    bool   operator >  ( const Self& it) const { return ( it1 >  it.it1); }
    bool   operator <= ( const Self& it) const { return ( it1 <= it.it1); }
    bool   operator >= ( const Self& it) const { return ( it1 >= it.it1); }

  private:
    RndAccIt   it1;
    Operation  op;
};

CGAL_END_NAMESPACE
  
#endif

// ===== EOF ==================================================================
