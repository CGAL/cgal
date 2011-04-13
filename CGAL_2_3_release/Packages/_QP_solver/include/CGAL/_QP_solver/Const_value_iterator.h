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
// file          : include/CGAL/_QP_solver/Const_value_iterator.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.2
// revision_date : 2000/08/29
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: sparse vector iterator
// ============================================================================

#ifndef CGAL_CONST_VALUE_ITERATOR_H
#define CGAL_CONST_VALUE_ITERATOR_H

#include <CGAL/basic.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class T >
class Const_value_iterator
#ifdef __BORLANDC__
    : public CGAL_STD::
        iterator<std::random_access_iterator_tag,T,ptrdiff_t,T*,T&>
#endif
                           {
  public:
    typedef  T                                value_type;
    typedef  ptrdiff_t                        difference_type;
    typedef  value_type*                      pointer;
    typedef  value_type&                      reference;
    typedef  std::random_access_iterator_tag  iterator_category;

    typedef  Const_value_iterator<T>  Self;
    typedef  value_type               Val;
    typedef  difference_type          Dist;
    typedef  reference                Ref;
    typedef  pointer                  Ptr;

    // forward operations
    Const_value_iterator( const T& t = T(), Dist i = 0)
	: index( i), value( t) { }

    bool       operator == ( const Self& it) const { return index == it.index;}
    bool       operator != ( const Self& it) const { return index != it.index;}

    Val        operator *  ( ) const { return  value; }
    Ptr        operator -> ( ) const { return &value; }

    Self&      operator ++ (    ) { ++index; return *this; }
    Self       operator ++ ( int) { Self tmp = *this; ++index; return tmp; }

    // bidirectional operations
    Self&      operator -- (    ) { --index; return *this; }
    Self       operator -- ( int) { Self tmp = *this; --index; return tmp; }

    // random access operations
    Self&  operator += ( Dist i) { index += i; return *this; }
    Self&  operator -= ( Dist i) { index -= i; return *this; }

    Self   operator +  ( Dist i) const
	{ Self tmp = *this; tmp += i; return tmp; }
    Self   operator -  ( Dist i) const
	{ Self tmp = *this; tmp -= i; return tmp; }

    Dist   operator -  ( const Self& it) const { return index - it.index; }

    Val    operator [] ( int) const { return value; }

    bool   operator <  ( const Self& it) const { return index <  it.index; }
    bool   operator >  ( const Self& it) const { return index >  it.index; }
    bool   operator <= ( const Self& it) const { return index <= it.index; }
    bool   operator >= ( const Self& it) const { return index >= it.index; }

  private:
    Dist  index;
    T     value;
};

CGAL_END_NAMESPACE
  
#endif

// ===== EOF ==================================================================
