// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/_QP_solver/Join_random_access_iterator.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.3
// revision_date : 2000/09/11
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: join random access iterator
// ============================================================================

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

    bool   operator <  ( const Self&) const { return ( it1 <  it.it1); }
    bool   operator >  ( const Self&) const { return ( it1 >  it.it1); }
    bool   operator <= ( const Self&) const { return ( it1 <= it.it1); }
    bool   operator >= ( const Self&) const { return ( it1 >= it.it1); }

  private:
    RndAccIt   it1;
    Operation  op;
};

CGAL_END_NAMESPACE
  
#endif

// ===== EOF ==================================================================
