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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_MIN_ANNULUS_D_H
#define CGAL_MIN_ANNULUS_D_H

#include <CGAL/license/Bounding_volumes.h>


#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244) // conversion warning in Boost iterator_adaptor
#endif

// includes
// --------
#include <CGAL/Optimisation/basic.h>
#include <CGAL/function_objects.h>
#include <CGAL/QP_options.h>
#include <CGAL/QP_solver/QP_solver.h>
#include <CGAL/QP_solver/functors.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/functional.hpp>

// here is how it works. We have d+2 variables: 
// R (big radius), r (small radius), c (center). The problem is
//
// min R^2 - r^2 
// s.t. ||p - c||^2 >= r^2   for all p
//      ||p - c||^2 <= R^2   for all p
//
// This looks nonlinear, but we can in fact make the substitutions
// u = R^2 - c^Tc, v = r^2 - c^Tc and get the following equivalent
// linear program:
//
// min u - v
// s.t. p^Tp - 2p^Tc >= v  for all p
//      p^Tp - 2p^Tc <= u  for all p
// 
// or 
//
// max  v - u
// s.t. v     + 2p_1c_1 + 2p_2c_2 + ...  + 2p_dc_d <=  p^Tp for all p
//        - u - 2p_1c_1 - 2p_2c_2 - ...  - 2p_dc_d <= -p^Tp for all p
//
// When we introduce a dual variable x_p for every constraint in the first
// set and a dual variable y_p for every constraint in the second set,
// we obtain the following dual program:
//
// min \sum_p x_p p^Tp -  \sum_p y_p p^Tp
// s.t.
//    2\sum_p x_p p    - 2\sum_p y_p p     =  0
//     \sum_p x_p                          =  1  (constraint for v)
//                     -  \sum_p y_p       = -1  (constraint for u)
//                                               x_p >= 0  for all p        
//                                               y_p >= 0  for all p 
//
// in the following functors, the ordering of the constraints is as above; 
// the indices of the variables are: x_p_j <-> 2 * j, y_p_j <-> 2 * j + 1 
// we also make the substitutions x'_p = x_p / h_p^2, y'_p = y_p / h_p^2
// where h_p is the homogenizing coordinate of p, in order to allow
// homogeneous points. This, however, means that the computed annulus is
// not necessarily correct. If P is a set of homogeneous points, 
//     P = { (p_0,...,p_{d-1}, h_p) }, 
// then we always get a *feasible* annulus for the point set 
//     P' = { (p_0*h_p,...,p_{d-1}*h_p, h_p*h_p) }. 
// If the type NT is inexact, this annulus might not even be optimal, since
// the objective function involves terms p^Tp that might not be exactly
// computed -> document all this!!!

namespace CGAL {

namespace MA_detail {

  // functor for a fixed column of A
  template <class NT, class Iterator>
  class A_column : public CGAL::unary_function <int, NT>
  {
  public:
    typedef NT result_type;
    A_column()
    {}
  
    A_column (int j, int d, Iterator it)
      : j_ (j), d_ (d), it_ (it), h_p (*(it+d)), nt_0_ (0), nt_2_ (2)
    {}

    result_type operator() (int i) const
    {
      if (j_ % 2 == 0) {
	// column for x_p
	if (i < d_) return  *(it_ + i) * h_p * nt_2_;  
	if (i == d_) return h_p * h_p;  
	return nt_0_;
      } else {
	// column for y_p
	if (i < d_) return  -(*(it_ + i)) * h_p * nt_2_;
	if (i == d_+1) return -h_p * h_p;
	return nt_0_;
      }
    }
    
  private:
    int j_;                  // column number
    int d_;                  // dimension
    Iterator it_;            // the iterator through the column's point
    NT h_p;                  // the homogenizing coordinate of p
    NT nt_0_;
    NT nt_2_; 
  };
  
  // functor for matrix A
  template <class NT, class Access_coordinate_begin_d,
	    class Point_iterator >
  class A_matrix : public CGAL::unary_function
  <int, boost::transform_iterator <A_column
    <NT, typename Access_coordinate_begin_d::Coordinate_iterator>, 
				   boost::counting_iterator<int> > >
  { 
    typedef typename MA_detail::A_column
    <NT, typename Access_coordinate_begin_d::Coordinate_iterator> A_column;
  public:
    typedef boost::transform_iterator
    <A_column, boost::counting_iterator<int> > result_type;
    
    A_matrix ()
    {}

    A_matrix (int d, 
	      const Access_coordinate_begin_d& da_coord,
	      Point_iterator P) 
      : d_ (d), da_coord_ (da_coord), P_ (P)
    {}

    result_type operator () (int j) const
    { 
      return result_type
	(0, A_column (j, d_, da_coord_ (*(P_+j/2)))); 
    }

  private:
    int d_;                  // dimension
    Access_coordinate_begin_d da_coord_; // data accessor
    Point_iterator P_;       // point set P
  };

  // The functor necessary to realize access to b
  template <class NT>
  class B_vector : public CGAL::unary_function<int, NT>
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
      if (i == d_) return nt_1_;
      if (i == d_+1) return -nt_1_;
      return nt_0_;
    }

  private:
    int d_;
    NT nt_0_;
    NT nt_1_;
  };

}

// Class interfaces
// ================
template < class Traits_ >
class Min_annulus_d {
public:
  // self
  typedef  Traits_                    Traits;
  typedef  Min_annulus_d<Traits>      Self;

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

  // public types 
  typedef  std::vector<Point>         Point_vector;
  typedef  typename Point_vector::const_iterator
  Point_iterator;

private:  
  // QP solver iterator types
  typedef MA_detail::A_matrix <NT, Access_coordinates_begin_d,
			       Point_iterator> A_matrix;
  typedef boost::transform_iterator
  <A_matrix,  
   boost::counting_iterator<int> >   A_iterator;

  typedef MA_detail::B_vector <NT> B_vector;
  typedef  boost::transform_iterator
  <B_vector,
    boost::counting_iterator<int> >  B_iterator;

  typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> R_iterator;  
 
  typedef  std::vector<NT>                                     C_vector; 
  typedef  typename C_vector::const_iterator                   C_iterator; 

  // Program type
  typedef CGAL::Nonnegative_linear_program_from_iterators
  <A_iterator, B_iterator, R_iterator, C_iterator> LP;

  // Tags
  typedef QP_solver_impl::QP_tags <Tag_true, Tag_true> QP_tags;

  // Solver types
  typedef  CGAL::QP_solver <LP, ET, QP_tags > Solver;

  typedef  typename Solver::Pricing_strategy Pricing_strategy;

  // types from the QP solver
  typedef  typename Solver::Basic_variable_index_iterator
  Basic_variable_index_iterator;
    
  // private types
  typedef  std::vector<ET>            ET_vector;
    
  typedef  QP_access_by_index
  <typename std::vector<Point>::const_iterator, int> Point_by_index;
    
  typedef  boost::binder2nd< std::divides<int> >
  Divide;
    
  typedef  std::vector<int>           Index_vector;
    
  typedef  std::vector<NT>            NT_vector;
  typedef  std::vector<NT_vector>     NT_matrix;
    

public:
  // public types
    
  typedef  CGAL::Join_input_iterator_1<
    Basic_variable_index_iterator,
    CGAL::Unary_compose_1<Point_by_index,Divide> >
  Support_point_iterator;
    
    
  typedef  typename Index_vector::const_iterator IVCI;
  typedef  CGAL::Join_input_iterator_1<
    IVCI, Point_by_index >
  Inner_support_point_iterator;
  typedef  CGAL::Join_input_iterator_1<
    IVCI, Point_by_index >
  Outer_support_point_iterator;

  typedef IVCI Inner_support_point_index_iterator;
  typedef IVCI Outer_support_point_index_iterator;
    
  typedef  typename ET_vector::const_iterator
  Coordinate_iterator;
    

  // creation
  Min_annulus_d( const Traits&  traits  = Traits())
    : tco( traits), da_coord(tco.access_coordinates_begin_d_object()),
      d( -1), solver(0){}
    
  template < class InputIterator >
  Min_annulus_d( InputIterator  first,
		 InputIterator  last,
		 const Traits&  traits = Traits())
    : tco( traits), da_coord(tco.access_coordinates_begin_d_object()), 
      solver(0) {
    set( first, last);
  }

  ~Min_annulus_d() {
    if (solver)
      delete solver;
  }
    
  // access to point set
  int  ambient_dimension( ) const { return d; }
    
  int  number_of_points( ) const { return static_cast<int>(points.size()); }
    
  Point_iterator  points_begin( ) const { return points.begin(); }
  Point_iterator  points_end  ( ) const { return points.end  (); }
    
  // access to support points
  int
  number_of_support_points( ) const
  { return number_of_points() < 2 ?
      number_of_points() :
    solver->number_of_basic_variables(); }
    
  Support_point_iterator
  support_points_begin() const {
    CGAL_optimisation_assertion_msg(number_of_points() >= 2,
				    "support_points_begin: not enough points");
    return Support_point_iterator(
				  solver->basic_original_variable_indices_begin(),
				  CGAL::compose1_1(
						   Point_by_index( points.begin()),
						   boost::bind2nd( std::divides<int>(), 2)));
  }
    
  Support_point_iterator
  support_points_end() const {
    CGAL_optimisation_assertion_msg(number_of_points() >= 2,
				    "support_points_begin: not enough points");
    return Support_point_iterator(
				  solver->basic_original_variable_indices_end(),
				  CGAL::compose1_1(
						   Point_by_index( points.begin()),
						   boost::bind2nd( std::divides<int>(), 2)));
  }
    
  int  number_of_inner_support_points() const { return static_cast<int>(inner_indices.size());}
  int  number_of_outer_support_points() const { return static_cast<int>(outer_indices.size());}
    
  Inner_support_point_iterator
  inner_support_points_begin() const
  { return Inner_support_point_iterator(
					inner_indices.begin(),
					Point_by_index( points.begin())); }
    
  Inner_support_point_iterator
  inner_support_points_end() const
  { return Inner_support_point_iterator(
					inner_indices.end(),
					Point_by_index( points.begin())); }
    
  Outer_support_point_iterator
  outer_support_points_begin() const
  { return Outer_support_point_iterator(
					outer_indices.begin(),
					Point_by_index( points.begin())); }
    
  Outer_support_point_iterator
  outer_support_points_end() const
  { return Outer_support_point_iterator(
					outer_indices.end(),
					Point_by_index( points.begin())); }

  Inner_support_point_index_iterator
  inner_support_points_indices_begin() const
  { return inner_indices.begin(); }
    
  Inner_support_point_index_iterator
  inner_support_points_indices_end() const
  { return inner_indices.end(); }

  Outer_support_point_index_iterator
  outer_support_points_indices_begin() const
  { return outer_indices.begin(); }

  Outer_support_point_index_iterator
  outer_support_points_indices_end() const
  { return outer_indices.end(); }


  // access to center (rational representation)
  Coordinate_iterator
  center_coordinates_begin( ) const { return center_coords.begin(); }
    
  Coordinate_iterator
  center_coordinates_end  ( ) const { return center_coords.end  (); }
    
  // access to squared radii (rational representation)
  ET  squared_inner_radius_numerator( ) const { return sqr_i_rad_numer; }
  ET  squared_outer_radius_numerator( ) const { return sqr_o_rad_numer; }
  ET  squared_radii_denominator     ( ) const { return sqr_rad_denom; }
    
  // access to center and squared radii
  // NOTE: an implicit conversion from ET to RT must be available!
  Point
  center( ) const
  { CGAL_optimisation_precondition( ! is_empty());
  return tco.construct_point_d_object()( ambient_dimension(),
					 center_coordinates_begin(),
					 center_coordinates_end()); }
    
  FT
  squared_inner_radius( ) const
  { CGAL_optimisation_precondition( ! is_empty());
  return FT( squared_inner_radius_numerator()) /
    FT( squared_radii_denominator()); }
    
  FT
  squared_outer_radius( ) const
  { CGAL_optimisation_precondition( ! is_empty());
  return FT( squared_outer_radius_numerator()) /
    FT( squared_radii_denominator()); }
    
  // predicates
  CGAL::Bounded_side
  bounded_side( const Point& p) const
  { CGAL_optimisation_precondition(
				   is_empty() || tco.access_dimension_d_object()( p) == d);
  ET sqr_d = sqr_dist( p);
  ET h_p_sqr = da_coord(p)[d] * da_coord(p)[d];
  return CGAL::Bounded_side
    (CGAL_NTS sign( sqr_d - h_p_sqr * sqr_i_rad_numer)
     * CGAL_NTS sign( h_p_sqr * sqr_o_rad_numer - sqr_d)); }
    
  bool
  has_on_bounded_side( const Point& p) const
  { CGAL_optimisation_precondition(
				   is_empty() || tco.access_dimension_d_object()( p) == d);
  ET sqr_d = sqr_dist( p);
  ET h_p_sqr = da_coord(p)[d] * da_coord(p)[d];
  return ( ( h_p_sqr * sqr_i_rad_numer < sqr_d) && 
	   ( sqr_d < h_p_sqr * sqr_o_rad_numer)); }
    
  bool
  has_on_boundary( const Point& p) const
  { CGAL_optimisation_precondition(
				   is_empty() || tco.access_dimension_d_object()( p) == d);
  ET sqr_d = sqr_dist( p);
  ET h_p_sqr = da_coord(p)[d] * da_coord(p)[d];
  return (( sqr_d == h_p_sqr * sqr_i_rad_numer) || 
	  ( sqr_d == h_p_sqr * sqr_o_rad_numer));}
    
  bool
  has_on_unbounded_side( const Point& p) const
  { CGAL_optimisation_precondition(
				   is_empty() || tco.access_dimension_d_object()( p) == d);
  ET sqr_d = sqr_dist( p);
  ET h_p_sqr = da_coord(p)[d] * da_coord(p)[d];
  return ( ( sqr_d < h_p_sqr * sqr_i_rad_numer) || 
	   ( h_p_sqr * sqr_o_rad_numer < sqr_d)); }
    
  bool  is_empty     ( ) const { return number_of_points() == 0; }
  bool  is_degenerate( ) const
  { return ! CGAL_NTS is_positive( sqr_o_rad_numer); }
    
  // modifiers
  template < class InputIterator >
  void
  set( InputIterator first, InputIterator last)
  { if ( points.size() > 0) points.erase( points.begin(), points.end());
  std::copy( first, last, std::back_inserter( points));
  set_dimension();
  CGAL_optimisation_precondition_msg( check_dimension(),
				      "Not all points have the same dimension.");
  compute_min_annulus(); }
    
  void
  insert( const Point& p)
  { if ( is_empty()) d = tco.access_dimension_d_object()( p);
  CGAL_optimisation_precondition(
				 tco.access_dimension_d_object()( p) == d);
  points.push_back( p);
  compute_min_annulus(); }
    
  template < class InputIterator >
  void
  insert( InputIterator first, InputIterator last)
  { CGAL_optimisation_precondition_code( std::size_t old_n = points.size());
  points.insert( points.end(), first, last);
  set_dimension();
  CGAL_optimisation_precondition_msg( check_dimension( old_n),
				      "Not all points have the same dimension.");
  compute_min_annulus(); }
    
  void
  clear( )
  { points.erase( points.begin(), points.end());
  compute_min_annulus(); }
    
  // validity check
  bool  is_valid( bool verbose = false, int level = 0) const;
    
  // traits class access
  const Traits&  traits( ) const { return tco; }
    

private:
    
  Traits                   tco;       // traits class object
  Access_coordinates_begin_d da_coord; // data accessor
    
  Point_vector             points;    // input points
  int                      d;         // dimension of input points
    
  ET_vector                center_coords;     // center of small.encl.annulus
    
  ET                       sqr_i_rad_numer;   // squared inner radius of
  ET                       sqr_o_rad_numer;   // ---"--- outer ----"----
  ET                       sqr_rad_denom;     // smallest enclosing annulus
    
  Solver                   *solver;    // linear programming solver
     
  Index_vector             inner_indices;
  Index_vector             outer_indices;
    
  NT_matrix                a_matrix;  // matrix `A' of dual LP
  NT_vector                b_vector;  // vector `b' of dual LP
  NT_vector                c_vector;  // vector `c' of dual LP
    
private:
  // squared distance to center * h_p^2 * c_d^2 
  ET sqr_dist_exact( const Point& p) const
  {
    ET result(0);
    typename Access_coordinates_begin_d::Coordinate_iterator 
      p_it (da_coord ( p));     // this is p * h_p
    ET c_d = center_coords[d];
    ET h_p = p_it[d];
    for (int i=0; i<d; ++i) {
      ET x = 
	c_d * ET(p_it[i]) -      /* this is c_d *    p_i * h_p */
	h_p * center_coords[i]   /* this is h_p *    c_i * c_d */ ;
      result += x * x;
    }
    return result;
  }

  // the function above computes sqr_dist as ||p-c||^2 
  // (endowed with a factor of c_d^2 * h_p^2)
  // but we know that c was computed from (possibly slightly wrong)
  // data if NT is inexact; in order to compensate for this, let
  // us instead compute sqr_dist as p^Tp - 2c^Tp + c^Tc, where we use
  // the (potentially wrong) values of p^Tp and p that went into the
  // linear program; this will give us correct containment / on_boundary 
  // checks also in the inexact-NT case.
  ET sqr_dist( const Point& p) const
  {
    ET result(0), two(2);
    NT pTp(0); // computed over input type, possibly slightly wrong
    ET cTc(0); 
    ET two_pTc(0);
    typename Access_coordinates_begin_d::Coordinate_iterator 
      p_it (da_coord ( p));     // this is p * h_p
    NT h_p = p_it[d]; // input type!
    for (int i=0; i<d; ++i) {
      NT p_i (p_it[i]);  // input type!
      pTp += p_i * p_i;  // p_i^2 * h_p^2 
      cTc += center_coords[i] * center_coords[i]; // c_i^2 * c_d^2
      two_pTc += 
	// 2 * c_i * c_d * p_i * h_p^2  
	2 * center_coords[i] * ET(h_p * p_i); 
    } 
    ET c_d = center_coords[d];
    result = 
      ET(pTp) * c_d * c_d +      
      cTc * ET (h_p * h_p) +
      - two_pTc * c_d;
    return result;
  }
    
  // set dimension of input points
  void
  set_dimension( )
  { d = ( points.size() == 0 ? -1 :
	  tco.access_dimension_d_object()( points[ 0])); }
    
  // check dimension of input points
  bool
  check_dimension( std::size_t  offset = 0)
  { return ( std::find_if( points.begin()+offset, points.end(),
			   CGAL::compose1_1( boost::bind2nd(
							  std::not_equal_to<int>(), d),
					     tco.access_dimension_d_object()))
	     == points.end()); }
    
  // compute smallest enclosing annulus
  void
  compute_min_annulus( )
  {
    // clear inner and outer support points
    inner_indices.erase( inner_indices.begin(), inner_indices.end());
    outer_indices.erase( outer_indices.begin(), outer_indices.end());
    
    if ( is_empty()) {
      center_coords.resize( 1);
      sqr_i_rad_numer = -ET( 1);
      sqr_o_rad_numer = -ET( 1);
      return;
    }
    
    if ( number_of_points() == 1) {
      inner_indices.push_back( 0);
      outer_indices.push_back( 0);
      center_coords.resize( d+1);
      std::copy( da_coord( points[ 0]),
		 da_coord( points[ 0])+d+1,
		 center_coords.begin());
      sqr_i_rad_numer = ET( 0);
      sqr_o_rad_numer = ET( 0);
      sqr_rad_denom   = ET( 1);
      return;
    }
    
    // set up vector c and solve dual LP
    // the ordering of the constraints is as above; the ordering
    // of the variables is: z_p_j <-> 2 * j, y_p_j <-> 2 * j + 1 
    c_vector.resize( 2*points.size());
    for ( int j = 0; j < number_of_points(); ++j) {
      typename Traits::Access_coordinates_begin_d::Coordinate_iterator
	coord_it = da_coord( points[j]);
      NT  sum = 0;
      for ( int i = 0; i < d; ++i) {
	sum += NT( coord_it[ i])*NT( coord_it[ i]);
      }
      c_vector[ 2*j  ] =  sum;
      c_vector[ 2*j+1] = -sum;
    }
        
    LP lp (2*static_cast<int>(points.size()), d+2, 
	   A_iterator ( boost::counting_iterator<int>(0), 
			A_matrix (d, da_coord, points.begin())),
	   B_iterator ( boost::counting_iterator<int>(0), 
			B_vector (d)), 
	   R_iterator (CGAL::EQUAL), 
	   c_vector.begin());
    
    Quadratic_program_options options;
    options.set_pricing_strategy(pricing_strategy(NT()));
    delete solver;
    solver = new Solver(lp, options);
    CGAL_optimisation_assertion(solver->status() == QP_OPTIMAL);
 
    // compute center and squared radius
    ET sqr_sum = 0;
    center_coords.resize( ambient_dimension()+1);
    for ( int i = 0; i < d; ++i) {
      center_coords[ i] = -solver->dual_variable( i);
      sqr_sum += center_coords[ i] * center_coords[ i];
    }
    center_coords[ d] = solver->variables_common_denominator();
    sqr_i_rad_numer = sqr_sum
      - solver->dual_variable( d  )*center_coords[ d];
    sqr_o_rad_numer = sqr_sum
      - solver->dual_variable( d+1)*center_coords[ d];
    sqr_rad_denom   = center_coords[ d] * center_coords[ d];
        
    // split up support points
    for ( int i = 0; i < solver->number_of_basic_original_variables(); ++i) {
      int index = solver->basic_original_variable_indices_begin()[ i];
      if ( index % 2 == 0) {
	inner_indices.push_back( index/2);
      } else {
	outer_indices.push_back( index/2);
      }
    }
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
// I/O operators
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os, const Min_annulus_d<Traits_>& min_annulus);

template < class Traits_ >
std::istream&
operator >> ( std::istream& is,       Min_annulus_d<Traits_>& min_annulus);

// ============================================================================

// Class implementation
// ====================

// validity check
template < class Traits_ >
bool
Min_annulus_d<Traits_>::
is_valid( bool verbose, int level) const
{
  using namespace std;

  CGAL::Verbose_ostream verr( verbose);
  verr << "CGAL::Min_annulus_d<Traits>::" << endl;
  verr << "is_valid( true, " << level << "):" << endl;
  verr << "  |P| = " << number_of_points()
       << ", |S| = " << number_of_support_points() << endl;

  // containment check (a)
  // ---------------------
  verr << "  (a) containment check..." << flush;
    
  Point_iterator  point_it = points_begin();
  for ( ; point_it != points_end(); ++point_it) {
    if ( has_on_unbounded_side( *point_it))
      return CGAL::_optimisation_is_valid_fail( verr,
						"annulus does not contain all points");
  }
    
  verr << "passed." << endl;

  // support set check (b)
  // ---------------------
  verr << "  (b) support set check..." << flush;
    
  // all inner support points on inner boundary?
  Inner_support_point_iterator  i_pt_it = inner_support_points_begin();
  for ( ; i_pt_it != inner_support_points_end(); ++i_pt_it) {
    ET h_p_sqr = da_coord (*i_pt_it)[d] * da_coord (*i_pt_it)[d];
    if ( sqr_dist( *i_pt_it) != h_p_sqr * sqr_i_rad_numer)
      return CGAL::_optimisation_is_valid_fail( verr,
						"annulus does not have all inner support points on its inner boundary");
  }
    
  // all outer support points on outer boundary?
  Outer_support_point_iterator  o_pt_it = outer_support_points_begin();
  for ( ; o_pt_it != outer_support_points_end(); ++o_pt_it) {
    ET h_p_sqr = da_coord (*o_pt_it)[d] * da_coord (*o_pt_it)[d];
    if ( sqr_dist( *o_pt_it) != h_p_sqr * sqr_o_rad_numer)
      return CGAL::_optimisation_is_valid_fail( verr,
						"annulus does not have all outer support points on its outer boundary");
  }
  /*
  // center strictly in convex hull of support points?
  typename Solver::Basic_variable_numerator_iterator
  num_it = solver.basic_variables_numerator_begin();
  for ( ; num_it != solver.basic_variables_numerator_end(); ++num_it) {
  if ( ! (    CGAL_NTS is_positive( *num_it)
  && *num_it <= solver.variables_common_denominator()))
  return CGAL::_optimisation_is_valid_fail( verr,
  "center does not lie strictly in convex hull of support points");
  }
  */
    
  verr << "passed." << endl;

  verr << "  object is valid!" << endl;
  return( true);
}

// output operator
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Min_annulus_d<Traits_>& min_annulus)
{
  using namespace std;

  typedef  typename Min_annulus_d<Traits_>::Point  Point;
  typedef  ostream_iterator<Point>       Os_it;
  typedef  typename Traits_::ET          ET;
  typedef  ostream_iterator<ET>          Et_it;

  switch ( CGAL::get_mode( os)) {

  case CGAL::IO::PRETTY:
    os << "CGAL::Min_annulus_d( |P| = "
       << min_annulus.number_of_points() << ", |S| = "
       << min_annulus.number_of_inner_support_points() << '+'
       << min_annulus.number_of_outer_support_points() << endl;
    os << "  P = {" << endl;
    os << "    ";
    copy( min_annulus.points_begin(), min_annulus.points_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  S_i = {" << endl;
    os << "    ";
    copy( min_annulus.inner_support_points_begin(),
	  min_annulus.inner_support_points_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  S_o = {" << endl;
    os << "    ";
    copy( min_annulus.outer_support_points_begin(),
	  min_annulus.outer_support_points_end(),
	  Os_it( os, ",\n    "));
    os << "}" << endl;
    os << "  center = ( ";
    copy( min_annulus.center_coordinates_begin(),
	  min_annulus.center_coordinates_end(),
	  Et_it( os, " "));
    os << ")" << endl;
    os << "  squared inner radius = "
       << min_annulus.squared_inner_radius_numerator() << " / "
       << min_annulus.squared_radii_denominator() << endl;
    os << "  squared outer radius = "
       << min_annulus.squared_outer_radius_numerator() << " / "
       << min_annulus.squared_radii_denominator() << endl;
    break;

  case CGAL::IO::ASCII:
    copy( min_annulus.points_begin(), min_annulus.points_end(),
	  Os_it( os, "\n"));
    break;

  case CGAL::IO::BINARY:
    copy( min_annulus.points_begin(), min_annulus.points_end(),
	  Os_it( os));
    break;

  default:
    CGAL_optimisation_assertion_msg( false,
				     "CGAL::get_mode( os) invalid!");
    break; }

  return( os);
}

// input operator
template < class Traits_ >
std::istream&
operator >> ( std::istream& is, CGAL::Min_annulus_d<Traits_>& min_annulus)
{
  using namespace std;

  switch ( CGAL::get_mode( is)) {

  case CGAL::IO::PRETTY:
    cerr << endl;
    cerr << "Stream must be in ascii or binary mode" << endl;
    break;

  case CGAL::IO::ASCII:
  case CGAL::IO::BINARY:
    typedef  typename CGAL::Min_annulus_d<Traits_>::Point  Point;
    typedef  istream_iterator<Point>             Is_it;
    min_annulus.set( Is_it( is), Is_it());
    break;

  default:
    CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
    break; }

  return( is);
}

} //namespace CGAL

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#endif // CGAL_MIN_ANNULUS_D_H

// ===== EOF ==================================================================
