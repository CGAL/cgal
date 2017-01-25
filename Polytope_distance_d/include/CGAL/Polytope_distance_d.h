// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
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
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_POLYTOPE_DISTANCE_D_H
#define CGAL_POLYTOPE_DISTANCE_D_H

#include <CGAL/license/Polytope_distance_d.h>


#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244) // conversion warning in Boost iterator_adaptor
#endif

// includes
// --------
#include <CGAL/Optimisation/basic.h>
#include <CGAL/function_objects.h>
#include <boost/functional.hpp>

#include <CGAL/QP_options.h>
#include <CGAL/QP_solver/QP_solver.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_solver/functors.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>

namespace CGAL {

// To determine the polytope distance, we set up the following QP:
// 
// We let v-w be nonnegative vectors such that
// v-w = p-q, where p and q are the points that realize the
// optimal distance. That way, we obtain a nonnegative QP
// 
//
//min (v-w)^T(w-v)
//     v-w - \sum x_i p_i + \sum y_j q_j    = 0
//           \sum x_i                       = 1
//                           \sum y_j       = 1
//     v, w, x, y >= 0

namespace PD_detail {
  // The functors necessary to realize access to A
  // 
  // A has the form
  //    I  -I    -P    Q
  //    0   0     1    0
  //    0   0     0    1
  // where  I and -I are blocks of size d*d, and the 0's and 1's are
  // rows vectors of 0's and 1's 
  //
  // we have one functor for a fixed column, and one functor for the
  // whole matrix

  // functor for a fixed column of A
  template <class NT, class Iterator>
  class A_column : public std::unary_function <int, NT>
  {
  public:
    typedef NT result_type;
    A_column()
    {}
  
    A_column (int j, int d, bool in_p, Iterator it)
      : j_ (j), d_ (d), in_p_ (in_p), it_ (it), nt_0_ (0), nt_1_ (1)
    {}

    result_type operator() (int i) const
    {
      if (j_ < d_) {
	// column for v
	return (i == j_ ? nt_1_ : nt_0_);
      }
      if (j_ < 2*d_) {
	// column for w
	return (i == j_-d_ ? -nt_1_ : nt_0_);
      }
      // point column
      if (in_p_) {
	// column for P
	if (i < d_) return -(*(it_ + i));
	if (i == d_) return *(it_ + d_); // homogenizing coordinate
	return nt_0_;
      } else {
	// column for Q
	if (i < d_) return (*(it_ + i));
	if (i == d_+1) return *(it_ + d_); // homogenizing coordinate
	return nt_0_;
      }
      // never get here
    }
    
  private:
    int j_;                  // column number
    int d_;                  // dimension
    bool in_p_;              // point in P ?
    Iterator it_;            // the iterator through the column's point
    NT nt_0_;
    NT nt_1_; 
  };
  
  // functor for matrix A
  template <class NT, class Access_coordinate_begin_d,
	    class Point_iterator >
  class A_matrix : public std::unary_function
  <int, boost::transform_iterator <A_column
    <NT, typename Access_coordinate_begin_d::Coordinate_iterator>, 
				   boost::counting_iterator<int> > >
  { 
    typedef PD_detail::A_column
    <NT, typename Access_coordinate_begin_d::Coordinate_iterator> A_column;
  public:
    typedef  boost::transform_iterator
    <A_column, boost::counting_iterator<int> > result_type;
    
    A_matrix ()
    {}

    A_matrix (int d, 
	      Access_coordinate_begin_d da_coord,
	      Point_iterator P, 
	      int r, 
	      Point_iterator Q) 
      : d_ (d), da_coord_ (da_coord), P_ (P), r_ (r), Q_ (Q)
    {}

    result_type operator () (int j) const
    { 
      if (j < 2*d_) {
	// column of v or w
	return result_type
	  (0, A_column (j, d_, false  /*dummy*/, 
			da_coord_ (*P_) /*dummy*/));
      }
      if (j < 2*d_+r_) {
	// column of P
	return result_type
	  (0, A_column (j , d_, true, da_coord_ (*(P_+(j-2*d_)))));
      }
      // column of Q
      return result_type
	(0, A_column (j, d_, false, da_coord_ (*(Q_+(j-2*d_-r_))))); 
    }

  private:
    int d_;                  // dimension
    Access_coordinate_begin_d da_coord_; // data accessor
    Point_iterator P_;       // point set P
    int r_;                  // size of point set P
    Point_iterator Q_;       // point set Q
  };

  // The functor necessary to realize access to b
  // b has the form
  // 0
  // 0
  // 1  (row d)
  // 1  (row d+1)

  template <class NT>
  class B_vector : public std::unary_function<int, NT>
  {
  public:
    typedef NT result_type;
    B_vector()
    {}

    B_vector (int d)
      : d_ (d), nt_0_ (0), nt_1_ (1)
    {}

    result_type operator() (int i) const
    {
      return ( (i == d_ || i == d_+1) ? nt_1_ : nt_0_);
    }

  private:
    int d_;
    NT nt_0_;
    NT nt_1_;
  };


  // The functors necessary to realize access to D
  // 
  // D has the form
  //   I  -I    0
  //  -I   I    0
  //   0   0    0
  // where all I and -I are blocks of size d*d
  //
  // we have one functor for a fixed row, and one functor for the
  // whole matrix

  // functor for a fixed row of D; note that we have to return 2D in
  // order to please the QP_solver
  template <class NT> 
  class D_row : public std::unary_function <int, NT> 
  {
  public:
    typedef NT result_type;
    D_row () 
    {}
    D_row (int i, int d) 
      : i_ (i), d_ (d), nt_0_ (0), nt_2_ (2)
    {}
    
    result_type operator () (int j) const
    {
      if (j < d_) {
	//  I ( 1 iff i = j)
	// -I (-1 iff i = j + d)
	//  0
	if (i_ == j) return nt_2_;
	if (i_ == j + d_) return -nt_2_;
	return nt_0_;
      }
      if (j < 2*d_) {
	// -I (-1 iff i = j - d)
	//  I ( 1 iff i = j)
	//  0
	if (i_ == j - d_) return -nt_2_;
	if (i_ == j) return nt_2_;
	return nt_0_;
      }
      // 0
      // 0
      // 0
      return nt_0_;
    }
  private:
    int i_; // row number
    int d_; // dimension
    NT nt_0_;
    NT nt_2_;
  };

  // functor for matrix D
  template <class NT>
  class D_matrix : public std::unary_function
  <int, boost::transform_iterator<D_row<NT>,
				  boost::counting_iterator<int> > >
  { 
  public:
    typedef boost::transform_iterator<D_row<NT>,
	    boost::counting_iterator<int> > result_type; 
    D_matrix ()
    {}
    D_matrix (int d)
      :  d_ (d)
    {}

    result_type operator()(int i) const 
    {
      return result_type (0, D_row<NT>(i, d_));
    }

  private:
    int d_; // dimension
  };
}

// Class interfaces
// ================
template < class Traits_ >
class Polytope_distance_d {
public:
  // self
  typedef  Traits_                    Traits;
  typedef  Polytope_distance_d<Traits>
  Self;

  // types from the traits class
  typedef  typename Traits::Point_d   Point;

  typedef  typename Traits::Rep_tag   Rep_tag;

  typedef  typename Traits::RT        RT;
  typedef  typename Traits::FT        FT;

  typedef  typename Traits::Access_dimension_d
  Access_dimension_d;
  typedef  typename Traits::Access_coordinates_begin_d
  Access_coordinates_begin_d;

  typedef  typename Traits::Construct_point_d
  Construct_point_d;

  typedef  typename Traits::ET        ET;
  typedef  typename Traits::NT        NT; 

private:
  // private types
  typedef  std::vector<Point>         Point_vector;
  typedef  std::vector<ET>            ET_vector;
    
  typedef  QP_access_by_index
  <typename std::vector<Point>::const_iterator, int> Point_by_index;
    
  typedef  std::vector<NT>            NT_vector;
  typedef  std::vector<NT_vector>     NT_matrix;
  
  typedef  std::vector<int>           Index_vector;
public:    
  // public types
  typedef  typename Point_vector::const_iterator
  Point_iterator;
    
  typedef typename Index_vector::const_iterator IVCI;
  typedef CGAL::Join_input_iterator_1< IVCI, Point_by_index >
  Support_point_iterator;

  typedef typename Index_vector::const_iterator
  Support_point_index_iterator;
    
  typedef  typename ET_vector::const_iterator
  Coordinate_iterator;

private:
  
  // QP solver iterator types
  typedef PD_detail::A_matrix <NT, Access_coordinates_begin_d,
			       Point_iterator> A_matrix;
  typedef boost::transform_iterator<
    A_matrix, boost::counting_iterator<int> > A_iterator;

  typedef PD_detail::B_vector <NT> B_vector;
  typedef boost::transform_iterator<
    B_vector, boost::counting_iterator<int> > B_iterator;

  typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_iterator;  
  typedef CGAL::Const_oneset_iterator<NT> C_iterator; 

  typedef PD_detail::D_matrix <NT> D_matrix;
  typedef boost::transform_iterator <
    D_matrix, boost::counting_iterator<int> > D_iterator;

  // Program type
  typedef CGAL::Nonnegative_quadratic_program_from_iterators
      <A_iterator, B_iterator, R_iterator, D_iterator, C_iterator> QP;

  // Tags
  typedef QP_solver_impl::QP_tags <Tag_false, Tag_true> QP_tags;

  // Solver types
  typedef  CGAL::QP_solver <QP, ET, QP_tags > Solver;

  typedef  typename Solver::Pricing_strategy Pricing_strategy;

public:    

  // creation
  Polytope_distance_d( const Traits&  traits  = Traits())
    : nt_0(0), nt_1(1),
      tco( traits), da_coord (tco.access_coordinates_begin_d_object()),
      d( -1), solver(0) 
  {}
    
  template < class InputIterator1, class InputIterator2 >
  Polytope_distance_d( InputIterator1 p_first,
		       InputIterator1 p_last,
		       InputIterator2 q_first,
		       InputIterator2 q_last,
		       const Traits&  traits = Traits())
    : nt_0(0), nt_1(1), 
      tco( traits), da_coord (tco.access_coordinates_begin_d_object()),
      solver(0)
  {
    set( p_first, p_last, q_first, q_last);
  }

  ~Polytope_distance_d() {    
    if (solver)
      delete solver;
  }

  // access to point sets
  int  ambient_dimension( ) const { return d; }
    
  int  number_of_points( ) const { return static_cast<int>(p_points.size()+q_points.size());}
    
  int  number_of_points_p( ) const { return static_cast<int>(p_points.size()); }
  int  number_of_points_q( ) const { return static_cast<int>(q_points.size()); }
    
  Point_iterator  points_p_begin( ) const { return p_points.begin(); }
  Point_iterator  points_p_end  ( ) const { return p_points.end  (); }
    
  Point_iterator  points_q_begin( ) const { return q_points.begin(); }
  Point_iterator  points_q_end  ( ) const { return q_points.end  (); }
    
  // access to support points
  int
  number_of_support_points( ) const
  { return is_finite() ? 
      static_cast<int>(p_support_indices.size()+q_support_indices.size()) : 0; }
    
  int  number_of_support_points_p() const { return static_cast<int>(p_support_indices.size());}
  int  number_of_support_points_q() const { return static_cast<int>(q_support_indices.size());}
    
  Support_point_iterator
  support_points_p_begin() const
  { return Support_point_iterator(
				  p_support_indices.begin(),
				  Point_by_index( p_points.begin())); }
    
  Support_point_iterator
  support_points_p_end() const
  { return Support_point_iterator(
				  is_finite() ? p_support_indices.end()
				  : p_support_indices.begin(),
				  Point_by_index( p_points.begin())); }
    
  Support_point_iterator
  support_points_q_begin() const
  { return Support_point_iterator(
				  q_support_indices.begin(),
				  Point_by_index( q_points.begin())); }
    
  Support_point_iterator
  support_points_q_end() const
  { return Support_point_iterator(
				  is_finite() ? q_support_indices.end()
				  : q_support_indices.begin(),
				  Point_by_index( q_points.begin())); }

  Support_point_index_iterator
  support_points_p_indices_begin() const
  { return p_support_indices.begin(); }
    
  Support_point_index_iterator
  support_points_p_indices_end() const
  { return p_support_indices.end(); }

  Support_point_index_iterator
  support_points_q_indices_begin() const
  { return q_support_indices.begin(); }

  Support_point_index_iterator
  support_points_q_indices_end() const
  { return q_support_indices.end(); }
    
  // access to realizing points (rational representation)
  Coordinate_iterator
  realizing_point_p_coordinates_begin( ) const { return p_coords.begin(); }
    
  Coordinate_iterator
  realizing_point_p_coordinates_end  ( ) const { return p_coords.end  (); }
    
  Coordinate_iterator
  realizing_point_q_coordinates_begin( ) const { return q_coords.begin(); }
    
  Coordinate_iterator
  realizing_point_q_coordinates_end  ( ) const { return q_coords.end  (); }
    
  // access to squared distance (rational representation)
  ET  squared_distance_numerator  ( ) const
  { return solver->solution_numerator(); }
    
  ET  squared_distance_denominator( ) const
  { return solver->solution_denominator(); }
    
  // access to realizing points and squared distance
  // NOTE: an implicit conversion from ET to RT must be available!
  Point
  realizing_point_p( ) const
  { CGAL_optimisation_precondition( is_finite());
  return tco.construct_point_d_object()
    ( ambient_dimension(),
      realizing_point_p_coordinates_begin(),
      realizing_point_p_coordinates_end  ()); }
    
  Point
  realizing_point_q( ) const
  { CGAL_optimisation_precondition( is_finite());
  return tco.construct_point_d_object()
    ( ambient_dimension(),
      realizing_point_q_coordinates_begin(),
      realizing_point_q_coordinates_end  ()); }
    
  FT
  squared_distance( ) const
  { 
    return FT( squared_distance_numerator  ()) /
      FT( squared_distance_denominator()); }
    
  bool  is_finite( ) const
  { return ( number_of_points_p() > 0) && ( number_of_points_q() > 0); }
    
  bool  is_zero( ) const
  { return CGAL_NTS is_zero( squared_distance_numerator()); }
    
  bool  is_degenerate( ) const { return ( ! is_finite()); }
    
  // modifiers
  template < class InputIterator1, class InputIterator2 >
  void
  set( InputIterator1 p_first, InputIterator1 p_last,
       InputIterator2 q_first, InputIterator2 q_last)
  { 
    p_points.clear();
    q_points.clear();
    std::copy( p_first, p_last, std::back_inserter( p_points));
    std::copy( q_first, q_last, std::back_inserter( q_points));
    set_dimension();
    CGAL_optimisation_precondition_msg
      (check_dimension( p_points.begin(), p_points.end())
       && check_dimension( q_points.begin(), q_points.end()),
       "Not all points have the same dimension.");

    compute_distance(); 
  }
    
  template < class InputIterator >
  void
  set_p( InputIterator p_first, InputIterator p_last)
  { 
    p_points.clear();
    std::copy( p_first, p_last, std::back_inserter( p_points));
    set_dimension();
    CGAL_optimisation_precondition_msg
      (check_dimension( p_points.begin(), p_points.end()),
       "Not all points have the same dimension.");

    compute_distance(); 
  }
    
  template < class InputIterator >
  void
  set_q( InputIterator q_first, InputIterator q_last)
  {
    q_points.clear();
    std::copy( q_first, q_last, std::back_inserter( q_points));
    set_dimension();
    CGAL_optimisation_precondition_msg
      (check_dimension( q_points.begin(), q_points.end()),
       "Not all points have the same dimension.");
  
    compute_distance(); 
  }
    
  void
  insert_p( const Point& p)
  { 
    CGAL_optimisation_precondition
      ( ( ! is_finite()) ||
	( tco.access_dimension_d_object()( p) == d));
    p_points.push_back( p);
    set_dimension(); // it might no longer be -1 
    compute_distance(); 
  }
    
  void
  insert_q( const Point& q)
  { 
    CGAL_optimisation_precondition
      ( ( ! is_finite()) ||
	( tco.access_dimension_d_object()( q) == d));
    q_points.push_back( q); 
    set_dimension(); // it might no longer be -1 
    compute_distance(); 
  }
    
  template < class InputIterator1, class InputIterator2 >
  void
  insert( InputIterator1 p_first, InputIterator1 p_last,
	  InputIterator2 q_first, InputIterator2 q_last)
  { 
    CGAL_optimisation_precondition_code(int old_r = static_cast<int>(p_points.size()));
    CGAL_optimisation_precondition_code(int old_s = static_cast<int>(q_points.size()));
    p_points.insert( p_points.end(), p_first, p_last);
    q_points.insert( q_points.end(), q_first, q_last);
    set_dimension();
    CGAL_optimisation_precondition_msg
      (check_dimension( p_points.begin()+old_r, p_points.end())
       && check_dimension( q_points.begin()+old_s, q_points.end()),
       "Not all points have the same dimension.");
    compute_distance(); 
  }
    
  template < class InputIterator >
  void
  insert_p( InputIterator p_first, InputIterator p_last)
  { 
    CGAL_optimisation_precondition_code(int old_r = static_cast<int>(p_points.size()));
    p_points.insert( p_points.end(), p_first, p_last);
    set_dimension();
    CGAL_optimisation_precondition_msg
      (check_dimension( p_points.begin()+old_r, p_points.end()),
       "Not all points have the same dimension.");
    compute_distance(); 
  }
    
  template < class InputIterator >
  void
  insert_q( InputIterator q_first, InputIterator q_last)
  { 
    CGAL_optimisation_precondition_code( int old_s = static_cast<int>(q_points.size()));
    q_points.insert( q_points.end(), q_first, q_last);
    set_dimension();
    CGAL_optimisation_precondition_msg
      (check_dimension( q_points.begin()+old_s, q_points.end()),
       "Not all points have the same dimension.");
    compute_distance(); 
  }
    
  void
  clear( )
  { 
    p_points.clear();
    q_points.clear();
    compute_distance(); 
  }
    
  // validity check
  bool  is_valid( bool verbose = false, int level = 0) const;
    
  // traits class access
  const Traits&  traits( ) const { return tco; }
    

private:
  NT nt_0;
  NT nt_1;

  Traits                     tco;        // traits class object
  Access_coordinates_begin_d da_coord;   // data accessor object
    
  Point_vector             p_points;  // points of P
  Point_vector             q_points;  // points of Q
  int                      d;         // dimension of input points
    
  ET_vector                p_coords;          // realizing point of P
  ET_vector                q_coords;          // realizing point of Q
    
  Solver                   *solver;    // quadratic programming solver
    
  Index_vector             p_support_indices;
  Index_vector             q_support_indices;
    
private:    
  // set dimension of input points
  void
  set_dimension( )
  { 
    d = ( p_points.size() > 0 ?
	  tco.access_dimension_d_object()( p_points[ 0]) :
	  q_points.size() > 0 ?
	  tco.access_dimension_d_object()( q_points[ 0]) :
	  -1); 
  }
    
  // check dimension of input points
  template < class InputIterator >
  bool
  check_dimension( InputIterator first, InputIterator last)
  { return ( std::find_if
	     ( first, last,
	       CGAL::compose1_1
	       ( boost::bind2nd(std::not_equal_to<int>(), d),
		 tco.access_dimension_d_object()))
	     == last); }
    
  // compute (squared) distance
  void
  compute_distance( )
  {
    // clear support points
    p_support_indices.clear();
    q_support_indices.clear();
    if ( ( p_points.size() == 0) || ( q_points.size() == 0)) return;
        
    // construct program
    int n = 2 * d + static_cast<int>(p_points.size() + q_points.size());
    int m = d + 2;
    CGAL_optimisation_precondition (p_points.size() > 0);
    QP qp (n, m, 
	   A_iterator 
	   (boost::counting_iterator<int>(0), 
	    A_matrix (d, da_coord, p_points.begin(), static_cast<int>(p_points.size()), 
			 q_points.begin())),
	   B_iterator (boost::counting_iterator<int>(0), B_vector (d)), 
	   R_iterator (CGAL::EQUAL), 
	   D_iterator (boost::counting_iterator<int>(0), D_matrix (d)),
	   C_iterator (nt_0));

    delete solver;
    Quadratic_program_options options;
    options.set_pricing_strategy(pricing_strategy(NT()));
    solver = new Solver(qp, options);
    CGAL_optimisation_assertion(solver->status() == QP_OPTIMAL);
    // compute support and realizing points
    ET  et_0 = 0;
    int r = static_cast<int>(p_points.size());
    p_coords.resize( ambient_dimension()+1);
    q_coords.resize( ambient_dimension()+1);
    std::fill( p_coords.begin(), p_coords.end(), et_0);
    std::fill( q_coords.begin(), q_coords.end(), et_0);
    for (int i = 0; i < solver->number_of_basic_original_variables(); ++i) {
      ET  value = solver->basic_original_variables_numerator_begin()[ i];
      int index = solver->basic_original_variable_indices_begin()[ i];
      if (index < 2*d) continue; // v or w variable
      if (index < 2*d + r) {
	// a point of p
	for ( int j = 0; j < d; ++j) {
	  p_coords[ j]
	    += value * ET(da_coord (p_points[ index-2*d ])[ j]);
	}
	p_support_indices.push_back( index-2*d);
      } else {
	// a point of q
	for ( int j = 0; j < d; ++j) {
	  q_coords[ j]
	    += value * ET(da_coord (q_points[ index-2*d-r])[ j]);
	}
	q_support_indices.push_back( index-2*d-r);
      }
    }
    p_coords[ d] = q_coords[ d] = solver->variables_common_denominator();
  }
    
  template < class NT >
  Quadratic_program_pricing_strategy pricing_strategy( NT) {
    return QP_PARTIAL_FILTERED_DANTZIG;
  }
  
  Quadratic_program_pricing_strategy pricing_strategy( ET) {
    return QP_PARTIAL_DANTZIG;
  }
    
};

// Function declarations
// =====================
// output operator
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Polytope_distance_d<Traits_>& poly_dist);


// ============================================================================

// Class implementation
// ====================

// validity check
template < class Traits_ >
bool
Polytope_distance_d<Traits_>::
is_valid( bool verbose, int level) const
{
  using namespace std;

  CGAL::Verbose_ostream verr( verbose);
  verr << "CGAL::Polytope_distance_d<Traits>::" << endl;
  verr << "is_valid( true, " << level << "):" << endl;
  verr << "  |P+Q| = " << number_of_points_p()
       <<          '+' << number_of_points_q()
       <<   ", |S| = " << number_of_support_points_p()
       <<          '+' << number_of_support_points_q() << endl;

  if ( is_finite()) {

    // compute normal vector
    ET_vector  normal( d), diff( d);
    ET  et_0 = 0, den = p_coords[d];
    CGAL_optimisation_assertion (den > et_0);
    CGAL_optimisation_assertion (den == q_coords[d]);
    int i, j;
    for ( j = 0; j < d; ++j) normal[ j] = p_coords[ j] - q_coords[ j];

    // check correctness of computed squared distance
    verr << "  checking squared_distance..." << flush;
    ET sqr_dist_num (0);
    for ( j = 0; j < d; ++j)
      sqr_dist_num += normal[ j] * normal[ j];
    ET sqr_dist_den = 
      p_coords[ d] * p_coords[ d];
    if (sqr_dist_num * squared_distance_denominator() !=
	sqr_dist_den * squared_distance_numerator())
      return CGAL::_optimisation_is_valid_fail
	( verr, "realizing points don't have correct squared distance");


    // check P
    // -------
    verr << "  checking P..." << flush;
        
    // check point set
    for ( i = 0; i < number_of_points_p(); ++i) {
      for ( j = 0; j < d; ++j) {
	// compute (a positive multiple of) p^* - p_i
	diff[ j] = 
	  ET(da_coord(p_points[ i])[d]) * p_coords[ j] - 
	  den * ET(da_coord( p_points[ i])[ j]);
      }
      if ( std::inner_product( diff.begin(), diff.end(),
			       normal.begin(), et_0) > et_0)
	return CGAL::_optimisation_is_valid_fail
	  ( verr, "polytope P is not separated by its hyperplane");
    }
        
    verr << "passed." << endl;

    // check Q
    // -------
    verr << "  checking Q..." << flush;
        
    // check point set
    for ( i = 0; i < number_of_points_q(); ++i) {
      for ( j = 0; j < d; ++j) {
	// compute (a positive multiple of) q^* - q_i
	diff[ j] = 
	  ET(da_coord(q_points[ i])[d]) * q_coords[ j] - 
	  den * ET(da_coord( q_points[ i])[ j]);
      }
      if ( std::inner_product( diff.begin(), diff.end(),
			       normal.begin(), et_0) < et_0)
	return CGAL::_optimisation_is_valid_fail
	  ( verr, "polytope Q is not separated by its hyperplane");
    }
        
    verr << "passed." << endl;
  }

  verr << "  object is valid!" << endl;
  return( true);
}

// output operator
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Polytope_distance_d<Traits_>& poly_dist)
{
  using namespace std;

  typedef  typename Polytope_distance_d<Traits_>::Point  Point;
  typedef  ostream_iterator<Point>       Os_it;
  typedef  typename Traits_::ET          ET;
  typedef  ostream_iterator<ET>          Et_it;

  switch ( CGAL::get_mode( os)) {

  case CGAL::IO::PRETTY:
    os << "CGAL::Polytope_distance_d( |P+Q| = "
       << poly_dist.number_of_points_p() << '+'
       << poly_dist.number_of_points_q() << ", |S| = "
       << poly_dist.number_of_support_points_p() << '+'
       << poly_dist.number_of_support_points_q() << endl;
    os << "  P = {" << endl;
    os << "    ";
    copy( poly_dist.points_p_begin(), poly_dist.points_p_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  Q = {" << endl;
    os << "    ";
    copy( poly_dist.points_q_begin(), poly_dist.points_q_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  S_P = {" << endl;
    os << "    ";
    copy( poly_dist.support_points_p_begin(),
	  poly_dist.support_points_p_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  S_Q = {" << endl;
    os << "    ";
    copy( poly_dist.support_points_q_begin(),
	  poly_dist.support_points_q_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  p = ( ";
    copy( poly_dist.realizing_point_p_coordinates_begin(),
	  poly_dist.realizing_point_p_coordinates_end(),
	  Et_it( os, " "));
    os << ")" << endl;
    os << "  q = ( ";
    copy( poly_dist.realizing_point_q_coordinates_begin(),
	  poly_dist.realizing_point_q_coordinates_end(),
	  Et_it( os, " "));
    os << ")" << endl;
    os << "  squared distance = "
       << poly_dist.squared_distance_numerator() << " / "
       << poly_dist.squared_distance_denominator() << endl;
    break;

  case CGAL::IO::ASCII:
    os << poly_dist.number_of_points_p() << endl;
    copy( poly_dist.points_p_begin(),
	  poly_dist.points_p_end(),
	  Os_it( os, "\n"));
    os << poly_dist.number_of_points_q() << endl;
    copy( poly_dist.points_q_begin(),
	  poly_dist.points_q_end(),
	  Os_it( os, "\n"));
    break;

  case CGAL::IO::BINARY:
    os << poly_dist.number_of_points_p() << endl;
    copy( poly_dist.points_p_begin(),
	  poly_dist.points_p_end(),
	  Os_it( os));
    os << poly_dist.number_of_points_q() << endl;
    copy( poly_dist.points_q_begin(),
	  poly_dist.points_q_end(),
	  Os_it( os));
    break;

  default:
    CGAL_optimisation_assertion_msg( false,
				     "CGAL::get_mode( os) invalid!");
    break; }

  return( os);
}


} //namespace CGAL

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#endif // CGAL_POLYTOPE_DISTANCE_D_H

// ===== EOF ==================================================================
