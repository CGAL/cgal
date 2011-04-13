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
// file          : include/CGAL/Optimisation_ellipse_2.h
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Min_ellipse_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>, Bernd Gärtner
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: 2D Optimisation Ellipse
// ============================================================================

#ifndef CGAL_OPTIMISATION_ELLIPSE_2_H
#define CGAL_OPTIMISATION_ELLIPSE_2_H

// the following include is needed by `to_double()'
#ifndef CGAL_CARTESIAN_H
#  include <CGAL/Cartesian.h>
#endif

// includes
#ifndef CGAL_POINT_2_H
#  include <CGAL/Point_2.h>
#endif
#ifndef CGAL_CONIC_2_H
#  include <CGAL/Conic_2.h>
#endif
#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#  include <CGAL/Optimisation/assertions.h>
#endif
#ifndef CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
#  include <CGAL/IO/forward_decl_window_stream.h>
#endif

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class K_ >
class Optimisation_ellipse_2;

// Class interface
// ===============
template < class K_ >
class Optimisation_ellipse_2 {
    /*
    friend  std::ostream&  operator << CGAL_NULL_TMPL_ARGS (
        std::ostream&, const Optimisation_ellipse_2<K_>&);
    friend  std::istream&  operator >> CGAL_NULL_TMPL_ARGS (
        std::istream&, Optimisation_ellipse_2<K_> &);
    friend  CGAL::Window_stream& operator << CGAL_NULL_TMPL_ARGS (
        CGAL::Window_stream&, const Optimisation_ellipse_2<K_>&);
    */
  public:
    // types
    typedef           K_                K;
    typedef  typename K_::RT            RT;
    typedef  typename K_::FT            FT;
    typedef           CGAL::Point_2<K>  Point;
    typedef           CGAL::Conic_2<K>  Conic;
    
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
    Point  boundary_point1, boundary_point2;    // two boundary points
    Conic  conic1, conic2;                      // two conics
    RT     dr, ds, dt, du, dv, dw;              // the gradient vector
    

// ============================================================================

// Class implementation
// ====================

  public:
    // Constructor
    // -----------
    inline
    Optimisation_ellipse_2( )
    { }

    // Set functions
    // -------------
    inline
    void
    set( )
    {
        n_boundary_points = 0;
    }
    
    inline
    void
    set( const Point& p)
    {
        n_boundary_points = 1;
        boundary_point1   = p;
    }
    
    inline
    void
    set( const Point& p, const Point& q)
    {
        n_boundary_points = 2;
        boundary_point1   = p;
        boundary_point2   = q;
    }
    
    inline
    void
    set( const Point& p1, const Point& p2, const Point& p3)
    {
        n_boundary_points = 3;
        conic1.set_ellipse( p1, p2, p3);
    }
    
    inline
    void
    set( const Point& p1, const Point& p2, const Point& p3, const Point& p4)
    {
        n_boundary_points = 4;
        Conic::set_two_linepairs( p1, p2, p3, p4, conic1, conic2);
        dr = RT( 0);
        ds = conic1.r() * conic2.s() - conic2.r() * conic1.s(),
        dt = conic1.r() * conic2.t() - conic2.r() * conic1.t(),
        du = conic1.r() * conic2.u() - conic2.r() * conic1.u(),
        dv = conic1.r() * conic2.v() - conic2.r() * conic1.v(),
        dw = conic1.r() * conic2.w() - conic2.r() * conic1.w();
    }
    
    inline
    void
    set( const Point&, const Point&,
         const Point&, const Point&, const Point& p5)
    {
        n_boundary_points = 5;
        conic1.set( conic1, conic2, p5);
        conic1.analyse();
    }

    // Access functions
    // ----------------
    inline
    int
    number_of_boundary_points( ) const
    {
        return( n_boundary_points);
    }
    
    Conic_2< Cartesian< double > >
    to_double( ) const
    {
        CGAL_optimisation_precondition( ! is_degenerate());
    
        double t = 0.0;
    
        if ( n_boundary_points == 4)
            t = conic1.vol_minimum( dr, ds, dt, du, dv, dw);
    
        Conic_2<K> c( conic1);
        Conic_2< Cartesian<double> > e;
        e.set( CGAL::to_double( c.r()) + t*CGAL::to_double( dr),
               CGAL::to_double( c.s()) + t*CGAL::to_double( ds),
               CGAL::to_double( c.t()) + t*CGAL::to_double( dt),
               CGAL::to_double( c.u()) + t*CGAL::to_double( du),
               CGAL::to_double( c.v()) + t*CGAL::to_double( dv),
               CGAL::to_double( c.w()) + t*CGAL::to_double( dw));
    
        return( e);
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
            return( conic1 == e.conic1);
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
            return( conic1.convex_side( p));
          case 4: {
            Conic c;
            c.set( conic1, conic2, p);
            c.analyse();
            if ( ! c.is_ellipse()) {
                c.set_ellipse( conic1, conic2);
                c.analyse();
                return( c.convex_side( p)); }
            else {
                int tau_star = c.vol_derivative( dr, ds, dt, du, dv, dw);
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

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Optimisation_ellipse_2.C>
#endif

#endif // CGAL_OPTIMISATION_ELLIPSE_2_H

// ===== EOF ==================================================================
