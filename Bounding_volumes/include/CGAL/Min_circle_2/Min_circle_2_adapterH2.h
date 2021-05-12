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

#ifndef CGAL_MIN_CIRCLE_2_ADAPTERH2_H
#define CGAL_MIN_CIRCLE_2_ADAPTERH2_H

#include <CGAL/license/Bounding_volumes.h>


// includes
#  include <CGAL/Optimisation/basic.h>

namespace CGAL {

// Class declarations
// ==================
template < class Traits_ >
class Min_circle_2;

template < class PT_, class DA_ >
class Min_circle_2_adapterH2;

template < class PT_, class DA_ >
class _Min_circle_2_adapterH2__Circle;

// Class interface and implementation
// ==================================
template < class PT_, class DA_ >
class Min_circle_2_adapterH2 {
  public:
    // types
    typedef  PT_  PT;
    typedef  DA_  DA;

    // nested types
    typedef  PT                                            Point;
    typedef  CGAL::_Min_circle_2_adapterH2__Circle<PT,DA>  Circle;

  private:
    DA      dao;                                    // data accessor object
    Circle  circle;                                 // current circle
    friend  class CGAL::Min_circle_2< CGAL::Min_circle_2_adapterH2<PT,DA> >;

  public:
    // creation
    Min_circle_2_adapterH2( const DA& da = DA())
        : dao( da), circle( da)
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

// Nested type `Circle'
template < class PT_, class DA_ >
std::ostream&
operator << ( std::ostream&,
              const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

template < class PT_, class DA_ >
std::istream&
operator >> ( std::istream&,
              CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

template < class PT_, class DA_ >
class _Min_circle_2_adapterH2__Circle {
  public:
    // typedefs
    typedef  PT_  PT;
    typedef  DA_  DA;

    typedef  typename DA_::RT             RT;
    typedef           CGAL::Quotient<RT>  FT;

  private:
    // data members
    DA  dao;                                // data accessor object

    RT  center_hx;                          // center's hx-coordinate
    RT  center_hy;                          // center's hy-coordinate
    RT  center_hw;                          // center's hw-coordinate
    FT  sqr_rad;                            // squared radius

    // private member functions
    FT
    sqr_dist( const RT& phx, const RT& phy, const RT& phw,
              const RT& qhx, const RT& qhy, const RT& qhw) const
    {
        RT  dhx( phx*qhw - qhx*phw);
        RT  dhy( phy*qhw - qhy*phw);
        RT  dhw( phw*qhw);
        return( FT( dhx*dhx + dhy*dhy, dhw*dhw));
    }

    friend  std::ostream&  operator << <> ( std::ostream&,
        const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

    friend  std::istream&  operator >> <> ( std::istream&,
        CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

  public:
    // types
    typedef  PT  Point;
    typedef  FT  Distance;

    // creation
    _Min_circle_2_adapterH2__Circle( const DA& da) : dao( da) { }

    void  set( )
    {
        center_hx =  RT( 0);
        center_hy =  RT( 0);
        center_hw =  RT( 1);
        sqr_rad   = -FT( 1);
    }

    void  set( const Point& p)
    {
        dao.get( p, center_hx, center_hy, center_hw);
        sqr_rad = FT( 0);
    }

    void  set( const Point& p, const Point& q)
    {
        RT  phx;
        RT  phy;
        RT  phw;
        RT  qhx;
        RT  qhy;
        RT  qhw;

        dao.get( p, phx, phy, phw);
        dao.get( q, qhx, qhy, qhw);
        center_hx = ( phx*qhw + qhx*phw);
        center_hy = ( phy*qhw + qhy*phw);
        center_hw = ( phw*qhw * RT( 2));
        sqr_rad   = sqr_dist( phx, phy, phw,
                              center_hx, center_hy, center_hw);
    }

    void  set( const Point& p, const Point& q, const Point& r)
    {
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

        RT  qhx_phx( qhx*phw - phx*qhw);
        RT  qhy_phy( qhy*phw - phy*qhw);    // denominator: qhw*phw

        RT  rhx_phx( rhx*phw - phx*rhw);
        RT  rhy_phy( rhy*phw - phy*rhw);    // denominator: rhw*phw

        RT  phw2( phw*phw);
        RT  qhw2( qhw*qhw);
        RT  rhw2( rhw*rhw);

        RT  p2( phx*phx + phy*phy);         // denominator: phw2

        RT  q2_p2( ( qhx*qhx + qhy*qhy) * phw2 - p2 * qhw2);
                                            // denominator: qhw2*phw2

        RT  r2_p2( ( rhx*rhx + rhy*rhy) * phw2 - p2 * rhw2);
                                            // denominator: rhw2*phw2

        center_hx = q2_p2*rhy_phy * rhw - r2_p2*qhy_phy * qhw;
        center_hy = r2_p2*qhx_phx * qhw - q2_p2*rhx_phx * rhw;
        center_hw = ( qhx_phx*rhy_phy - rhx_phx*qhy_phy)
                      * phw*qhw*rhw * RT( 2);
        sqr_rad   = sqr_dist( phx, phy, phw,
                              center_hx, center_hy, center_hw);
    }

    // predicates
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        RT  phx;
        RT  phy;
        RT  phw;
        dao.get( p, phx, phy, phw);
        return( CGAL::Bounded_side( CGAL_NTS sign(
            sqr_rad - sqr_dist( phx, phy, phw,
                                center_hx, center_hy, center_hw))));
    }

    bool
    has_on_bounded_side( const Point& p) const
    {
        RT  phx;
        RT  phy;
        RT  phw;
        dao.get( p, phx, phy, phw);
        return( sqr_dist( phx, phy, phw,
                          center_hx, center_hy, center_hw) < sqr_rad);
    }

    bool
    has_on_boundary( const Point& p) const
    {
        RT  phx;
        RT  phy;
        RT  phw;
        dao.get( p, phx, phy, phw);
        return( sqr_dist( phx, phy, phw,
                          center_hx, center_hy, center_hw) == sqr_rad);
    }

    bool
    has_on_unbounded_side( const Point& p) const
    {
        RT  phx;
        RT  phy;
        RT  phw;
        dao.get( p, phx, phy, phw);
        return( sqr_rad < sqr_dist( phx, phy, phw,
                                    center_hx, center_hy, center_hw));
    }

    bool
    is_empty( ) const
    {
        return( CGAL::is_negative( sqr_rad));
    }

    bool
    is_degenerate( ) const
    {
        return( ! CGAL::is_positive( sqr_rad));
    }

    // additional operations for checking
    bool
    operator == (
        const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c) const
    {
        return( ( center_hx*c.center_hw == c.center_hx*center_hw) &&
                ( center_hy*c.center_hw == c.center_hy*center_hw) &&
                ( sqr_rad  == c.sqr_rad ) );
    }

    bool
    operator != (
        const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c) const
    {
        return( ! ( *this == c));
    }

    Point
    center( ) const
    {
        Point  p;
        dao.set( p, center_hx, center_hy, center_hw);
        return( p);
    }

    const Distance&
    squared_radius( ) const
    {
        return( sqr_rad);
    }
};

// I/O
template < class PT_, class DA_ >
std::ostream&
operator << ( std::ostream& os,
              const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c)
{
    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << "CGAL::Min_circle_2_adapterH2::Circle( "
           << c.center_hx << ", "
           << c.center_hy << ", "
           << c.center_hw << ", "
           << c.sqr_rad   << ')';
        break;

      case CGAL::IO::ASCII:
        os << c.center_hx << ' '
           << c.center_hy << ' '
           << c.center_hw << ' '
           << c.sqr_rad;
        break;

      case CGAL::IO::BINARY:
        CGAL::write( os, c.center_hx);
        CGAL::write( os, c.center_hy);
        CGAL::write( os, c.center_hw);
        CGAL::write( os, c.sqr_rad);
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                        "CGAL::get_mode( os) invalid!");
        break; }

    return( os);
}

template < class PT_, class DA_ >
std::istream&
operator >> ( std::istream& is,
              CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c)
{
    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII:
        is >> c.center_hx >> c.center_hy >> c.center_hw >> c.sqr_rad;
        break;

      case CGAL::IO::BINARY:
        CGAL::read( is, c.center_hx);
        CGAL::read( is, c.center_hy);
        CGAL::read( is, c.center_hw);
        CGAL::read( is, c.sqr_rad);
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

} //namespace CGAL

#endif // CGAL_MIN_CIRCLE_2_ADAPTERH2_H

// ===== EOF ==================================================================
