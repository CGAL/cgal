// Copyright (c) 1997-2001
// ETH Zurich (Switzerland).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

#ifndef CGAL_MIN_ELLIPSE_2_ADAPTERH2_H
#define CGAL_MIN_ELLIPSE_2_ADAPTERH2_H

#include <CGAL/license/Bounding_volumes.h>


// includes
#  include <CGAL/Homogeneous/ConicHPA2.h>
#  include <CGAL/Optimisation/assertions.h>

namespace CGAL {

// Class declarations
// ==================
template < class Traits_ >
class Min_ellipse_2;

template < class PT_, class DA_ >
class Min_ellipse_2_adapterH2;

template < class PT_, class DA_ >
class _Min_ellipse_2_adapterH2__Ellipse;

// Class interface and implementation
// ==================================
template < class PT, class DA >
bool
are_ordered_along_lineH2( const PT& p, const PT& q, const PT& r,
                          const DA& da)
{
    typedef  typename DA::RT  RT;

    RT  phx;
    RT  phy;
    RT  phw;
    RT  qhx;
    RT  qhy;
    RT  qhw;
    RT  rhx;
    RT  rhy;
    RT  rhw;

    da.get( p, phx, phy, phw);
    da.get( q, qhx, qhy, qhw);
    da.get( r, rhx, rhy, rhw);

    // p,q,r collinear?
    if ( ! CGAL_NTS is_zero(  ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
                            - ( phy*rhw - rhy*phw) * ( qhx*rhw - rhx*qhw)))
        return( false);

    // p,q,r vertical?
    if ( phx*rhw != rhx*phw)
        return(    ( ( phx*qhw < qhx*phw) && ( qhx*rhw < rhx*qhw))
                || ( ( rhx*qhw < qhx*rhw) && ( qhx*phw < phx*qhw)));
    else
        return(    ( ( phy*qhw < qhy*phw) && ( qhy*rhw < rhy*qhw))
                || ( ( rhy*qhw < qhy*rhw) && ( qhy*phw < phy*qhw)));
}

template < class PT_, class DA_ >
class Min_ellipse_2_adapterH2 {
  public:
    // types
    typedef  PT_  PT;
    typedef  DA_  DA;

    // nested types
    typedef  PT                                             Point;
    typedef  _Min_ellipse_2_adapterH2__Ellipse<PT,DA>  Ellipse;

  private:
    DA      dao;                                    // data accessor object
    Ellipse ellipse;                                // current ellipse
    friend class Min_ellipse_2< Min_ellipse_2_adapterH2<PT,DA> >;

  public:
    // creation
    Min_ellipse_2_adapterH2( const DA& da = DA())
        : dao( da), ellipse( da)
    { }

    // operations
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        typedef  typename DA_::RT  RT;

        RT  phx;
        RT  phy;
        RT  phw;
        RT  qhx;
        RT  qhy;
        RT  qhw;
        RT  rhx;
        RT  rhy;
        RT  rhw;

        dao.get( p, phx, phy, phw);
        dao.get( q, qhx, qhy, qhw);
        dao.get( r, rhx, rhy, rhw);

        return( static_cast< CGAL::Orientation>(
                 CGAL_NTS sign( ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
                              - ( phy*rhw - rhy*phw) * ( qhx*rhw - rhx*qhw))));
    }
};

// Nested type `Ellipse'
template < class PT_, class DA_ >
std::ostream&
operator << ( std::ostream&,
              const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

template < class PT_, class DA_ >
std::istream&
operator >> ( std::istream&,
              CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

template < class PT_, class DA_ >
class _Min_ellipse_2_adapterH2__Ellipse {
  public:
    // typedefs
    typedef  PT_  PT;
    typedef  DA_  DA;

    typedef           CGAL::ConicHPA2< PT, DA>  CT;
    typedef  typename DA_::RT                   RT;

  private:
    // data members
    int  n_boundary_points;                 // number of boundary points
    PT   boundary_point1, boundary_point2;  // two boundary points
    CT   conic1, conic2;                    // two conics
    RT   dr, ds, dt, du, dv, dw;            // the gradient vector

    friend  std::ostream&  operator << <> ( std::ostream&,
        const _Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

    friend  std::istream&  operator >> <> ( std::istream&,
        _Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>&);

  public:
    // types
    typedef  PT  Point;

    // creation
    _Min_ellipse_2_adapterH2__Ellipse( const DA& da)
        : conic1( da), conic2( da)
    { }

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
        boundary_point1   = p;
        boundary_point2   = q;
    }

    void
    set( const Point& p1, const Point& p2, const Point& p3)
    {
        n_boundary_points = 3;
        conic1.set_ellipse( p1, p2, p3);
    }

    void
    set( const Point& p1, const Point& p2,
         const Point& p3, const Point& p4)
    {
        n_boundary_points = 4;
        CT::set_two_linepairs( p1, p2, p3, p4, conic1, conic2);
        dr = RT( 0);
        ds = conic1.r() * conic2.s() - conic2.r() * conic1.s(),
        dt = conic1.r() * conic2.t() - conic2.r() * conic1.t(),
        du = conic1.r() * conic2.u() - conic2.r() * conic1.u(),
        dv = conic1.r() * conic2.v() - conic2.r() * conic1.v(),
        dw = conic1.r() * conic2.w() - conic2.r() * conic1.w();
    }

    void
    set( const Point&, const Point&,
         const Point&, const Point&, const Point& p5)
    {
        n_boundary_points = 5;
        conic1.set( conic1, conic2, p5);
        conic1.analyse();
    }

    // predicates
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
                    || ( CGAL::are_ordered_along_lineH2( boundary_point1, p,
                                           boundary_point2, conic1.da())) ?
                                CGAL::ON_BOUNDARY : CGAL::ON_UNBOUNDED_SIDE);
          case 3:
          case 5:
            return( conic1.convex_side( p));
          case 4: {
            CT c( conic1.da());
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

    bool
    has_on_bounded_side( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_BOUNDED_SIDE);
    }

    bool
    has_on_boundary( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_BOUNDARY);
    }

    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( bounded_side( p) == CGAL::ON_UNBOUNDED_SIDE);
    }

    bool
    is_empty( ) const
    {
        return( n_boundary_points == 0);
    }

    bool
    is_degenerate( ) const
    {
        return( n_boundary_points < 3);
    }

    // additional operations for checking
    bool
    operator == (
        const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e) const
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

    bool
    operator != (
        const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e) const
    {
        return( ! ( *this == e));
    }
};

// I/O
template < class PT_, class DA_ >
std::ostream&
operator << ( std::ostream& os,
              const CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e)
{
    const char* const  empty       = "";
    const char* const  pretty_head =
                             "CGAL::Min_ellipse_2_adapterH2::Ellipse( ";
    const char* const  pretty_sep  = ", ";
    const char* const  pretty_tail = ")";
    const char* const  ascii_sep   = " ";

    const char*  head = empty;
    const char*  sep  = empty;
    const char*  tail = empty;

    switch ( CGAL::get_mode( os)) {
      case CGAL::IO::PRETTY:
        head = pretty_head;
        sep  = pretty_sep;
        tail = pretty_tail;
        break;
      case CGAL::IO::ASCII:
        sep  = ascii_sep;
        break;
      case CGAL::IO::BINARY:
        break;
      default:
        CGAL_optimisation_assertion_msg( false,
                                        "CGAL::get_mode( os) invalid!");
        break; }

    os << head << e.n_boundary_points;
    switch ( e.n_boundary_points) {
      case 0:
        break;
      case 1:
        os << sep << e.boundary_point1;
        break;
      case 2:
        os << sep << e.boundary_point1
           << sep << e.boundary_point2;
        break;
      case 3:
      case 5:
        os << sep << e.conic1;
        break;
      case 4:
        os << sep << e.conic1
           << sep << e.conic2;
        break; }
    os << tail;

    return( os);
}

template < class PT_, class DA_ >
std::istream&
operator >> ( std::istream& is,
              CGAL::_Min_ellipse_2_adapterH2__Ellipse<PT_,DA_>& e)
{
    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        CGAL::read( is, e.n_boundary_points);
        switch ( e.n_boundary_points) {
          case 0:
            break;
          case 1:
            is >> e.boundary_point1;
            break;
          case 2:
            is >> e.boundary_point1
               >> e.boundary_point2;
            break;
          case 3:
          case 5:
            is >> e.conic1;
            break;
          case 4:
            is >> e.conic1
               >> e.conic2;
            break; }
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

} //namespace CGAL

#endif // CGAL_MIN_ELLIPSE_2_ADAPTERH2_H

// ===== EOF ==================================================================
