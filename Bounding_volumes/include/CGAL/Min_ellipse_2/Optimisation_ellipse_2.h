// Copyright (c) 1997-2001  
// ETH Zurich (Switzerland).  All rights reserved.
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

#ifndef CGAL_OPTIMISATION_ELLIPSE_2_H
#define CGAL_OPTIMISATION_ELLIPSE_2_H

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/Conic_2.h>
#include <CGAL/Optimisation/assertions.h>
#include <CGAL/Kernel/global_functions_2.h>


namespace CGAL {

// Class interface
// ===============
template < class K_ >
class Optimisation_ellipse_2 {
    /*
    friend  std::ostream&  operator << <> (
        std::ostream&, const Optimisation_ellipse_2<K_>&);
    friend  std::istream&  operator >> <> (
        std::istream&, Optimisation_ellipse_2<K_> &);
    */
  public:
    // types
    typedef           K_                K;
    typedef  typename K_::RT            RT;
    typedef  typename K_::FT            FT;
    typedef  typename K::Point_2        Point;
    typedef  typename K::Conic_2        Conic;
    
    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.
    
    // creation
    Optimisation_ellipse_2( );
    
    void  set( );
    void  set( const Point& p);
    void  set( const Point& p,  const Point& q);
    void  set( const Point& p1, const Point& p2, const Point& p3);
    void  set( const Point& p1, const Point& p2,
               const Point& p3, const Point& p4);
    void  set( const Point& p1, const Point& p2,
               const Point& p3, const Point& p4, const Point& p5);
    
    // access functions
    int  number_of_boundary_points()
    
    // equality tests
    bool  operator == ( const Optimisation_ellipse_2<K>& e) const;
    bool  operator != ( const Optimisation_ellipse_2<K>& e) const;
    
    // predicates
    CGAL::Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;
    
    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
    **************************************************************************/

  /* private: */
    // private data members
    int    n_boundary_points;                   // number of boundary points
    Point  boundary_point1, 
           boundary_point2,
           boundary_point3,
           boundary_point4,
           boundary_point5;                     // <= 5 support point 
    Conic  conic1, conic2;                      // two conics 

    // this gradient vector has dr=0 and is used in testing the
    // position of a point relative to an ellipse through 4 points
    mutable RT     dr, ds, dt, du, dv, dw;  
    mutable bool   d_values_set; 
    
    // this gradient vector is just conic2 - conic1 and is used in
    // obtaining an explicit conic representing an ellipse through 4 poinnts
    mutable RT     er, es, et, eu, ev, ew;
    mutable bool e_values_set;

    // needed in bounded-side predicate over ellipse with 4 support points
    mutable Conic helper_ellipse; // needed in bounded-side predicate over 
    mutable bool helper_ellipse_set;

    mutable Conic helper_conic; // also needed in bounded-side test

// ============================================================================

// Class implementation
// ====================

  public:
    // Constructor
    // -----------
    inline
    Optimisation_ellipse_2( ) { }

    // Set functions
    // -------------
    void
    set( )
    {
        n_boundary_points = 0;
    }
    
    void
    set( const Point& p)
    {
        n_boundary_points = 1;
        boundary_point1   = p;
    }
    
    void
    set( const Point& p, const Point& q)
    {
        n_boundary_points = 2;
        CGAL_optimisation_precondition(boundary_point1 == p); CGAL_USE(p);
        boundary_point2 = q;
    }
    
    void
    set( const Point& p1, const Point& p2, const Point& p3)
    {       
        n_boundary_points = 3;        
	CGAL_optimisation_precondition(boundary_point1 == p1);
        CGAL_optimisation_precondition(boundary_point2 == p2);
	boundary_point3 = p3;
        helper_conic.set_ellipse( p1, p2, p3);
	CGAL_optimisation_assertion(helper_conic.is_ellipse());
    }
    
    void
    set( const Point& p1, const Point& p2, const Point& p3, const Point& p4)
    {
        n_boundary_points = 4;	
	CGAL_optimisation_precondition(boundary_point1 == p1);
        CGAL_optimisation_precondition(boundary_point2 == p2);
	CGAL_optimisation_precondition(boundary_point3 == p3);
        boundary_point4 = p4;
        Conic::set_two_linepairs( p1, p2, p3, p4, conic1, conic2);

	d_values_set = false;
	e_values_set = false;
	helper_ellipse_set = false;	
    }

    void
    set_d_values() const
    {
      if (!d_values_set) {
        dr = RT( 0);
        ds = conic1.r() * conic2.s() - conic2.r() * conic1.s(),
        dt = conic1.r() * conic2.t() - conic2.r() * conic1.t(),
        du = conic1.r() * conic2.u() - conic2.r() * conic1.u(),
        dv = conic1.r() * conic2.v() - conic2.r() * conic1.v(),
        dw = conic1.r() * conic2.w() - conic2.r() * conic1.w();
	d_values_set = true;
      }
    }

    void
    set_e_values() const
    {
      if (!e_values_set) {
       	er = conic2.r() - conic1.r();
	es = conic2.s() - conic1.s();
	et = conic2.t() - conic1.t();
	eu = conic2.u() - conic1.u();
	ev = conic2.v() - conic1.v();
	ew = conic2.w() - conic1.w();
	e_values_set = true;
      }
    }

    void
    set_helper_ellipse () const
    {
      if (!helper_ellipse_set) {
        helper_ellipse.set_ellipse( conic1, conic2);
        helper_ellipse.analyse();
	CGAL_optimisation_assertion (helper_ellipse.is_ellipse());
	helper_ellipse_set= true;
      }
    }

    void
    set( const Point& p1, const Point& p2,
         const Point& p3, const Point& p4, const Point& p5)
    { 
        helper_conic.set(conic1, conic2, p5);
	helper_conic.analyse();
	// an optimization is possible if this set-call arose from
	// a successful violation test of ME(p1,p2,p3,p4) and p5.
	// In that case, helper_conic is already correct, 
	// but in general, this optimization is NOT valid.
	n_boundary_points = 5;
	CGAL_optimisation_assertion(helper_conic.is_ellipse());	
	// the following assertion is too strict if we run under
	// double (which is sometimes the case, e.g. in demos)
	// CGAL_optimisation_assertion(helper_conic.has_on_boundary(p5));
	CGAL_optimisation_precondition(boundary_point1 == p1);
        CGAL_optimisation_precondition(boundary_point2 == p2);
	CGAL_optimisation_precondition(boundary_point3 == p3);
        CGAL_optimisation_precondition(boundary_point4 == p4);
	CGAL_USE(p1); CGAL_USE(p2); CGAL_USE(p3); CGAL_USE(p4);
	boundary_point5 = p5;
    }

    // Access functions
    // ----------------
    inline
    int
    number_of_boundary_points( ) const
    {
        return( n_boundary_points);
    }
    
    template <typename DoubleConic_2>
    void
    double_conic(DoubleConic_2& e) const
    {
        double r,s,t,u,v,w;
	double_coefficients(r,s,t,u,v,w);
        e.set(r,s,t,u,v,w);
	// NOTE: the set method calls analyze, so the conic is clean
    }

    void 
    double_coefficients (double &r, double &s,double &t,
                         double &u, double &v, double &w) const
    {
      // just like double_conic, but we only get the coefficients
      CGAL_optimisation_precondition( ! is_degenerate());
    
      if ( n_boundary_points == 4) {
        set_e_values();
        double tau = conic1.vol_minimum( er, es, et, eu, ev, ew);
	r = CGAL::to_double( conic1.r()) + tau*CGAL::to_double( er);
	s = CGAL::to_double( conic1.s()) + tau*CGAL::to_double( es);
	t = CGAL::to_double( conic1.t()) + tau*CGAL::to_double( et);
	u = CGAL::to_double( conic1.u()) + tau*CGAL::to_double( eu);
	v = CGAL::to_double( conic1.v()) + tau*CGAL::to_double( ev);
	w = CGAL::to_double( conic1.w()) + tau*CGAL::to_double( ew);
      } else {
	// it's the helper_conic
	r = CGAL::to_double(helper_conic.r());
	s = CGAL::to_double(helper_conic.s());
	t = CGAL::to_double(helper_conic.t());
	u = CGAL::to_double(helper_conic.u());
	v = CGAL::to_double(helper_conic.v());
        w = CGAL::to_double(helper_conic.w());
      }
    }


    // Equality tests
    // --------------
    bool
    operator == ( const Optimisation_ellipse_2<K>& e) const
    {
        if ( n_boundary_points != e.n_boundary_points)
            return( false);
    
        switch ( n_boundary_points) {
          case 0:
            return( true);
          case 1:
            return( boundary_point1 == e.boundary_point1);
          case 2:
            return(    (    ( boundary_point1 == e.boundary_point1)
                         && ( boundary_point2 == e.boundary_point2))
                    || (    ( boundary_point1 == e.boundary_point2)
                         && ( boundary_point2 == e.boundary_point1)));
          case 3:
          case 5:
            return( helper_conic == e.helper_conic);
          case 4:
            return(    (    ( conic1 == e.conic1)
                         && ( conic2 == e.conic2))
                    || (    ( conic1 == e.conic2)
                         && ( conic2 == e.conic1)));
          default:
            CGAL_optimisation_assertion(    ( n_boundary_points >= 0)
                                         && ( n_boundary_points <= 5)); }
        // keeps g++ happy
        return( false);
    }
    
    inline
    bool
    operator != ( const Optimisation_ellipse_2<K>& e) const
    {
        return( ! operator == ( e));
    }

    // Predicates
    // ----------
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        switch ( n_boundary_points) {
          case 0:
            return( CGAL::ON_UNBOUNDED_SIDE);
          case 1:
            return( ( p == boundary_point1) ?
                           CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
          case 2:
            return(    ( p == boundary_point1)
                    || ( p == boundary_point2)
                    || ( CGAL::are_ordered_along_line(
                             boundary_point1, p, boundary_point2)) ?
                         CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
          case 3:
          case 5:
            return(helper_conic.convex_side( p));
          case 4: {
            helper_conic.set( conic1, conic2, p);
            helper_conic.analyse();
            if ( !helper_conic.is_ellipse()) {
	        set_helper_ellipse();
                return( helper_ellipse.convex_side( p)); }
            else {
	        set_d_values();
                int tau_star = 
                  helper_conic.vol_derivative( dr, ds, dt, du, dv, dw);
                return( CGAL::Bounded_side( CGAL_NTS sign( tau_star))); } }
          default:
            CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                         ( n_boundary_points <= 5) ); }
        // keeps g++ happy
        return( CGAL::Bounded_side( 0));
    }
    
    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_BOUNDED_SIDE);
    }
    
    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_BOUNDARY);
    }
    
    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_UNBOUNDED_SIDE);
    }
    
    inline
    bool
    is_empty( ) const
    {
        return( n_boundary_points == 0);
    }
    
    inline
    bool
    is_degenerate( ) const
    {
        return( n_boundary_points < 3);
    }

    bool
    is_circle( ) const
    {
       switch ( n_boundary_points) {
       case 0: 
	 return false; // the empty set is not a circle
       case 1:
	 return true;  
       case 2:
	 return false; // a segment is not a circle
       case 3:
       case 5:
	 return helper_conic.is_circle();
       case 4:
	 // the smallest ellipse through four points is
	 // a circle only if the four points are cocircular;
	 // if so, compute this circle (as a conic) and check
	 // its volume derivative
	 if (CGAL::ON_BOUNDARY !=  CGAL::side_of_bounded_circle
	           (boundary_point1, 
                    boundary_point2,
                    boundary_point3,
                    boundary_point4)) {
	   return false;
	 } else {
	   // ok, they are cocircular, now get the circle and check it
	   Conic c;
	   c.set_circle(boundary_point1, boundary_point2, boundary_point3);
           set_d_values();
	   return (CGAL::ZERO ==  (c.vol_derivative(dr, ds, dt, du, dv, dw)));
	 }
       default:
	 CGAL_optimisation_assertion( ( n_boundary_points >= 0) &&
                                      ( n_boundary_points <= 5) ); 
	 return false;
       }
    }
};

// Function declarations
// =====================
// I/O
// ---
template < class K_ >
std::ostream&
operator << ( std::ostream&, const CGAL::Optimisation_ellipse_2<K_>&);

template < class K_ >
std::istream&
operator >> ( std::istream&, CGAL::Optimisation_ellipse_2<K_>&);

} //namespace CGAL

#include <CGAL/Min_ellipse_2/Optimisation_ellipse_2_impl.h>

#endif // CGAL_OPTIMISATION_ELLIPSE_2_H

// ===== EOF ==================================================================
