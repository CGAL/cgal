// ============================================================================
//
// Copyright (c) 1997-2004 The CGAL Consortium
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
// file          : include/CGAL/QP_solver/iterator.h
// package       : $CGAL_Package: QP_solver $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Quadratic Programming Engine - Iterators
// ============================================================================

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
