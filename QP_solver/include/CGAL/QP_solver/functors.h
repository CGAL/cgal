// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr 
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp 
//                 Kaspar Fischer

#ifndef CGAL_QP_SOLVER_FUNCTORS_H
#define CGAL_QP_SOLVER_FUNCTORS_H

#include <CGAL/license/QP_solver.h>


#include <CGAL/QP_solver/basic.h>
#include <CGAL/function_objects.h>
#include <functional>
#include <iterator>

namespace CGAL {

// ==================
// class declarations
// ==================
template < class VectorIt,
           bool check_lower = false, bool check_upper = false >
class QP_vector_accessor;

template < class MatrixIt,
           bool check_1st_lower = false, bool check_1st_upper = false,
           bool check_2nd_lower = false, bool check_2nd_upper = false >
class QP_matrix_accessor;

template < class MatrixIt, typename ET >
class QP_matrix_pairwise_accessor;


template < class RndAccIt >
class Value_by_basic_index;

template<typename Map>
class Map_with_default;

// =====================
// class implementations
// =====================

// -------------------
// QP_vector_accessor
// -------------------
template < class VectorIt, bool check_lower, bool check_upper >
class QP_vector_accessor : public CGAL::unary_function<
    int, typename std::iterator_traits<VectorIt>::value_type > {

  public:
    typedef typename 
        CGAL::unary_function<
           int, 
           typename std::iterator_traits<VectorIt>::value_type >::result_type
    result_type;
    QP_vector_accessor( VectorIt it, int lower = 0, int upper = 0)
        : z( 0), v( it), l(lower), u(upper)
    {}

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
// QP_matrix_accessor
// -------------------
template < class MatrixIt, bool check_1st_lower, bool check_1st_upper,
                           bool check_2nd_lower, bool check_2nd_upper >
class QP_matrix_accessor {

  public:    
    typedef int argument1_type;
    typedef int argument2_type;
    typedef typename std::iterator_traits<MatrixIt>::value_type VectorIt;
    typedef typename std::iterator_traits<VectorIt>::value_type result_type;

    QP_matrix_accessor( MatrixIt it, int lower_1 = 0, int upper_1 = 0,
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
	return VectorIt(m[ r])[ c];
    }

  private:
    const result_type  z;
    MatrixIt           m;
    int   l1, u1, l2, u2;
};		 		 

// ----------------------------------------------------------------------------

// ----------------------------
// QP_matrix_pairwise_accessor
// ----------------------------
template < class MatrixIt, typename ResultType >
class QP_matrix_pairwise_accessor {
  typedef  typename std::iterator_traits<MatrixIt>::value_type  VectorIt;
  
public:
  typedef int        argument_type;
  typedef ResultType result_type;
  
  // The following default constructor is needed to make it possible
  // to use QP_matrix_pairwise_accessor with CGAL's Join_input_iterator_1
  // (more precisely: once Join_input_iterator_1 should not use an internal
  // mutable variable 'val' anymore, you can remove the following default
  // constructor).
  //QP_matrix_pairwise_accessor() {}
  
  QP_matrix_pairwise_accessor( MatrixIt it, int row)
    : m (it), v (*(it + row)), r (row)                                  
  {}
  
  ResultType operator () ( int c) const
  {
    // make sure that only entries on or below the diagonal are
    // accessed
    if (c <= r)
      return ResultType(v[ c]); 
    else
      return ResultType((*(m + c))[ r]);
  }
  
private:
  MatrixIt           m;
  VectorIt           v;
  int                r;
};


// ----------------------------------------------------------------------------

// --------------------
// Value_by_basic_index
// --------------------
template < class RndAccIt >
class Value_by_basic_index : public CGAL::unary_function<
    int, typename std::iterator_traits<RndAccIt>::value_type > {

  public:
    typedef typename
    CGAL::unary_function<
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

// ----------------------------------------------------------------------------

// --------------------
// Access_by_index
// --------------------
// A functor whose operator(int i) provides access to the i-th element
// of a random access iterator.
template < typename RndAccIt, typename ArgType >
class QP_access_by_index {
public:
  typedef  typename std::iterator_traits<RndAccIt>::value_type result_type;

  QP_access_by_index(RndAccIt it = RndAccIt()) : a(it) {}

  result_type operator () (ArgType i) const { return a[i]; }

private:
  RndAccIt     a;
};


// -------------------
// Map_with_default
// -------------------
template<typename Map>
class Map_with_default {
  // public types
public:
  typedef typename Map::mapped_type       mapped_type;
  //typedef typename Map::difference_type   difference_type;
  typedef typename Map::size_type         size_type;
  typedef typename Map::key_type          key_type;
  typedef mapped_type                     result_type;
  // data members
private:
  const Map* map; // pointer to map
  mapped_type d;  // default value
  
public:
  // construction
  Map_with_default ()
    : map(0), d()
  {}

  Map_with_default (const Map* m, const mapped_type& v = mapped_type())
    : map(m), d(v)
  {}
  
  // operator()
  const mapped_type& operator() (key_type n) const {
    CGAL_qpe_precondition (map != 0);
    typename Map::const_iterator i = map->find (n);
    if (i != map->end())
      return i->second;
    else
      return d;
  }
};

} //namespace CGAL

#endif // CGAL_QP_SOLVER_FUNCTORS_H

// ===== EOF ==================================================================
