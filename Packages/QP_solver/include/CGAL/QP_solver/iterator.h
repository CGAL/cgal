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
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/QP_engine/iterator.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2002/01/27
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Quadratic Programming Engine - Iterators
// ============================================================================

#ifndef CGAL_QP_ENGINE_ITERATOR_H
#define CGAL_QP_ENGINE_ITERATOR_H

// includes
#ifndef CGAL_QP_ENGINE_BASIC_H
#  include <CGAL/QP_engine/basic.h>
#endif
#ifndef CGAL_ITERATOR
#  include <CGAL/iterator.h>
#endif

CGAL_BEGIN_NAMESPACE

// ==================
// class declarations
// ==================
template < class T >
class QPE_const_value_iterator;

template < class Iterator, class Operation >
class QPE_transform_iterator_1;

template < class Iterator1, class Iterator2, class Operation >
class QPE_transform_iterator_2;

// =====================
// class implementations
// =====================

// ------------------------
// QPE_const_value_iterator
// ------------------------
template < class T >
class QPE_const_value_iterator {

  public:

    // types
    typedef  std::random_access_iterator_tag    iterator_category;
    typedef  ptrdiff_t                          difference_type;
    typedef  T                                  value_type;
    typedef  value_type*                        pointer;
    typedef  value_type&                        reference;

    typedef  QPE_const_value_iterator<T>        Self;
    typedef  difference_type                    Diff;
    typedef  value_type                         Val;
    typedef  pointer                            Ptr;
    typedef  reference                          Ref;

    // construction
    QPE_const_value_iterator( const T& t = T(), Diff n = 0)
	: value( t), index( n)
    { }

    // access
    Ref        operator *  ( )       { return  value; }
    const Ref  operator *  ( ) const { return  value; }
    Ptr        operator -> ( )       { return &value; }
    const Ptr  operator -> ( ) const { return &value; }

    // equality operator
    bool       operator == ( const Self& x) const { return ( index==x.index); }

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
    Ref        operator [] ( Diff i)       { return value; }
    const Ref  operator [] ( Diff i) const { return value; }

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

// ------------------------
// QPE_transform_iterator_1
// ------------------------
template < class Iterator, class Operation >
class QPE_transform_iterator_1 {

  public:

    // types
    typedef  typename std::iterator_traits<Iterator>::iterator_category
                                                iterator_category;
    typedef  typename std::iterator_traits<Iterator>::difference_type
                                                difference_type;
    typedef  typename Operation::result_type    value_type;
    typedef  value_type*                        pointer;
    typedef  value_type&                        reference;

    typedef  QPE_transform_iterator_1<Iterator,Operation>
                                                Self;
    typedef  difference_type                    Diff;
    typedef  value_type                         Val;

    // construction
    QPE_transform_iterator_1( const Iterator&   it_ = Iterator(),
			      const Operation&  op_ = Operation())
	: it( it_), op( op_)
    { }

    // access
    Val    operator *  ( ) const { return op( *it); }

    // equality operators
    bool   operator == ( const Self& x) const { return ( it == x.it); }
    bool   operator != ( const Self& x) const { return ( it != x.it); }

    // forward operations
    // ------------------
    Self&  operator ++ (    ) {                   ++it; return *this; }
    Self   operator ++ ( int) { Self tmp = *this; ++it; return tmp;   }

    // bidirectional operations
    // ------------------------
    Self&  operator -- (    ) {                   --it; return *this; }
    Self   operator -- ( int) { Self tmp = *this; --it; return tmp;   }

    // random access operations
    // ------------------------
    // access
    Val    operator [] ( Diff i) const { return op( it[ i]); }

    // less operator
    bool   operator <  ( const Self& x) const { return ( it < x.it); }

    // arithmetic operations
    Self&  operator += ( Diff n) { it += n; return *this; }
    Self&  operator -= ( Diff n) { it -= n; return *this; }

    Self   operator +  ( Diff n) const { Self tmp = *this; return tmp += n; }
    Self   operator -  ( Diff n) const { Self tmp = *this; return tmp -= n; }

    Diff   operator -  ( const Self& x) const { return it - x.it; }

  private:

    // data members
    Iterator   it;
    Operation  op;
};

template < class Iterator, class Operation > inline
QPE_transform_iterator_1<Iterator,Operation>
transformer( const Iterator& it, const Operation& op)
{
    return QPE_transform_iterator_1<Iterator,Operation>( it, op);
}

// ------------------------
// QPE_transform_iterator_2
// ------------------------
template < class Iterator1, class Iterator2, class Operation >
class QPE_transform_iterator_2 {

  public:

    // types
    typedef  typename std::iterator_traits<Iterator1>::iterator_category
                                                iterator_category;
    typedef  typename std::iterator_traits<Iterator1>::difference_type
                                                difference_type;
    typedef  typename Operation::result_type    value_type;
    typedef  value_type*                        pointer;
    typedef  value_type&                        reference;

    typedef  QPE_transform_iterator_2<Iterator1,Iterator2,Operation>
                                                Self;
    typedef  difference_type                    Diff;
    typedef  value_type                         Val;

    // construction
    QPE_transform_iterator_2( const Iterator1&  it1_ = Iterator1(),
			      const Iterator2&  it2_ = Iterator2(),
			      const Operation&  op_  = Operation())
	: it1( it1_), it2( it2_), op( op_)
    { }

    // access
    Val    operator *  ( ) const { return op( *it1, *it2); }

    // equality operator
    bool   operator == ( const Self& x) const
                              { return ( ( it1 == x.it1) && ( it2 == x.it2)); }

    // forward operations
    // ------------------
    Self&  operator ++ (    ) {                   ++it1; ++it2; return *this; }
    Self   operator ++ ( int) { Self tmp = *this; ++it1; ++it2; return tmp;   }

    // bidirectional operations
    // ------------------------
    Self&  operator -- (    ) {                   --it; --it2; return *this; }
    Self   operator -- ( int) { Self tmp = *this; --it; --it2; return tmp;   }

    // random access operations
    // ------------------------
    // access
    Val    operator [] ( Diff i) const { return op( it1[ i], it2[ i]); }

    // less operator
    bool   operator <  ( const Self& x) const { return ( it < x.it); }

    // arithmetic operations
    Self&  operator += ( Diff n) { it1 += n; it2 += n; return *this; }
    Self&  operator -= ( Diff n) { it1 -= n; it2 -= n; return *this; }

    Self   operator +  ( Diff n) const { Self tmp = *this; return tmp += n; }
    Self   operator -  ( Diff n) const { Self tmp = *this; return tmp -= n; }

    Diff   operator -  ( const Self& x) const { return it1 - x.it1; }

  private:

    // data members
    Iterator1  it1;
    Iterator2  it2;
    Operation  op;
};

template < class Iterator1, class Iterator2, class Operation > inline
QPE_transform_iterator_2<Iterator1,Iterator2,Operation>
transformer( const Iterator1& it1, const Iterator2& it2, const Operation& op)
{
    return QPE_transform_iterator_2<Iterator1,Iterator2,Operation>(it1,it2,op);
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_ENGINE_ITERATOR_H

// ===== EOF ==================================================================
