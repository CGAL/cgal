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
// file          : include/CGAL/QP_engine/functors.h
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Quadratic Programming Engine - Functors
// ============================================================================

#ifndef CGAL_QP_ENGINE_FUNCTORS_H
#define CGAL_QP_ENGINE_FUNCTORS_H

// includes
#ifndef CGAL_QP_ENGINE_BASIC_H
#  include <CGAL/QP_engine/basic.h>
#endif
#ifndef CGAL_FUNCTION_OBJECTS_H
#  include <CGAL/function_objects.h>
#endif

#ifndef CGAL_PROTECT_FUNCTIONAL
#  define CGAL_PROTECT_FUNCTIONAL
#  include <functional>
#endif
#ifndef CGAL_PROTECT_ITERATOR
#  define CGAL_PROTECT_ITERATOR
#  include <iterator>
#endif

CGAL_BEGIN_NAMESPACE

// ==================
// class declarations
// ==================
template < class VectorIt,
           bool check_lower = false, bool check_upper = false >
class QPE_vector_accessor;

template < class MatrixIt,
           bool check_1st_lower = false, bool check_1st_upper = false,
           bool check_2nd_lower = false, bool check_2nd_upper = false >
class QPE_matrix_accessor;

template < class MatrixIt, class IsSymmetric,
           bool check_lower = false, bool check_upper = false >
class QPE_matrix_pairwise_accessor;


template < class RndAccIt >
class Value_by_basic_index;


// =====================
// class implementations
// =====================

// -------------------
// QPE_vector_accessor
// -------------------
template < class VectorIt, bool check_lower, bool check_upper >
class QPE_vector_accessor : public std::unary_function<
    int, typename std::iterator_traits<VectorIt>::value_type > {

  public:
    typedef typename 
        std::unary_function<
           int, 
           typename std::iterator_traits<VectorIt>::value_type >::result_type
    result_type;
    QPE_vector_accessor( VectorIt it, int lower = 0, int upper = 0)
        : z( 0), v( it)
    {
	if ( check_lower) l = lower;
	if ( check_upper) u = upper;
    }

    result_type  operator ( ) ( int i) const
    {
	if ( check_lower && i <  l) return z;
	if ( check_upper && i >= u) return z;
	return v[ i];
    }

  private:
    const result_type  z;
    VectorIt           v;
    int                l;
    int                u;
};

// ----------------------------------------------------------------------------

// -------------------
// QPE_matrix_accessor
// -------------------
template < class MatrixIt, bool check_1st_lower, bool check_1st_upper,
                           bool check_2nd_lower, bool check_2nd_upper >
class QPE_matrix_accessor {

  public:    
    typedef CGAL::Arity_tag<2> Arity;
    typedef int argument1_type;
    typedef int argument2_type;
    typedef typename std::iterator_traits<
            typename std::iterator_traits<MatrixIt>::value_type>::value_type  
    result_type;

    QPE_matrix_accessor( MatrixIt it, int lower_1 = 0, int upper_1 = 0,
			              int lower_2 = 0, int upper_2 = 0)
	: z( 0), m( it)
    {
	if ( check_1st_lower) l1 = lower_1;
	if ( check_1st_upper) u1 = upper_1;
	if ( check_2nd_lower) l2 = lower_2;
	if ( check_2nd_upper) u2 = upper_2;
    }

    result_type  operator () ( int r, int c) const
    {
	if ( check_1st_lower && ( r <  l1)) return z;
	if ( check_1st_upper && ( r >= u1)) return z;
	if ( check_2nd_lower && ( c <  l2)) return z;
	if ( check_2nd_upper && ( c >= u2)) return z;
	return m[ r][ c];
    }

  private:
    const result_type  z;
    MatrixIt           m;
    int   l1, u1, l2, u2;
};		 		 

// ----------------------------------------------------------------------------

// ----------------------------
// QPE_matrix_pairwise_accessor
// ----------------------------
template < class MatrixIt,class IsSymmetric,bool check_lower,bool check_upper >
class QPE_matrix_pairwise_accessor : public std::unary_function<
    int, typename std::iterator_traits<
             typename std::iterator_traits<MatrixIt>::value_type>::value_type>{

    typedef  typename std::iterator_traits<MatrixIt>::value_type  VectorIt;

  public:
    typedef typename
    std::unary_function<
    int, typename std::iterator_traits<
             typename std::iterator_traits
                  <MatrixIt>::value_type>::value_type>::result_type
    result_type;

    // following operator= added by kf. to make Join_input_iterator_1 work
    QPE_matrix_pairwise_accessor() : z( 0)
    {
    }

    QPE_matrix_pairwise_accessor( MatrixIt it, int row, int lower = 0,
				                        int upper = 0)
	: z( 0)
    {
	if ( check_tag( IsSymmetric())) {               // store i-th row
	    v = it[ row];
	} else {                                        // matrix and index
	    m = it;
	    r = row;
	}
	if ( check_lower) l = lower;
	if ( check_upper) u = upper;
    }
    
    // following operator= added by kf. to make Join_input_iterator_1 work
    QPE_matrix_pairwise_accessor&
      operator=(const QPE_matrix_pairwise_accessor& a)
    {
      m = a.m;
      v = a.v;
      r = a.r;
      l = a.l;
      u = a.u;
      return *this;
    }
    
    result_type  operator () ( int c) const
    {
	if ( check_lower && ( ( r <  l) || ( c <  l))) return z;
	if ( check_upper && ( ( r >= u) || ( c >= u))) return z;
	return entry_pair( c, IsSymmetric());
    }

    result_type  entry_pair( int c, Tag_true ) const           // symmetric
        { return v[ c] * result_type( 2); }

    result_type  entry_pair( int c, Tag_false) const           // not symmetric
        { return m[ r][ c] + m[ c][ r];	}

  private:
    const result_type  z;
    MatrixIt           m;
    VectorIt           v;
    int                r;
    int                l;
    int                u;
};


// ----------------------------------------------------------------------------

// --------------------
// Value_by_basic_index
// --------------------
template < class RndAccIt >
class Value_by_basic_index : public std::unary_function<
    int, typename std::iterator_traits<RndAccIt>::value_type > {

  public:
    typedef typename
    std::unary_function<
      int, typename std::iterator_traits
         <RndAccIt>::value_type >::result_type
    result_type;

    Value_by_basic_index( RndAccIt x_B_O_it, int n_original)
	: o( x_B_O_it), s( x_B_O_it),
	  l( n_original), u( n_original+0),
	  z( 0)
	{ }

    Value_by_basic_index( RndAccIt x_B_O_it, int n_original,
			  RndAccIt x_B_S_it, int n_slack = 0)
	: o( x_B_O_it), s( x_B_S_it),
	  l( n_original), u( n_original+n_slack),
	  z( 0)
	{ }

    result_type  operator () ( int i) const
        {
	    if ( i < 0) return z;
	    if ( i >= l && i < u) return s[ i];
	    return o[ i];
	}

  private:
    RndAccIt     o, s;
    int          l, u;
    result_type  z;
};


CGAL_END_NAMESPACE

#endif // CGAL_QP_ENGINE_FUNCTORS_H

// ===== EOF ==================================================================
