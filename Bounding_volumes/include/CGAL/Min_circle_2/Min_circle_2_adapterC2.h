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

#ifndef CGAL_MIN_CIRCLE_2_ADAPTERC2_H
#define CGAL_MIN_CIRCLE_2_ADAPTERC2_H

#include <CGAL/license/Bounding_volumes.h>


// includes
#  include <CGAL/Optimisation/basic.h>

namespace CGAL {

// Class declarations
// ==================
template < class Traits_ >
class Min_circle_2;

template < class PT_, class DA_ >
class Min_circle_2_adapterC2;

template < class PT_, class DA_ >
class _Min_circle_2_adapterC2__Circle;

// Class interface and implementation
// ==================================
template < class PT_, class DA_ >
class Min_circle_2_adapterC2 {
  public:
    // types
    typedef  PT_  PT;
    typedef  DA_  DA;

    // nested types
    typedef  PT                                            Point;
    typedef  CGAL::_Min_circle_2_adapterC2__Circle<PT,DA>  Circle;

  private:
    DA      dao;                                    // data accessor object
    Circle  circle;                                 // current circle
    friend  class CGAL::Min_circle_2< CGAL::Min_circle_2_adapterC2<PT,DA> >;

  public:
    // creation
    Min_circle_2_adapterC2( const DA& da = DA())
        : dao( da), circle( da)
    { }

    // operations
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        typedef  typename DA_::FT  FT;
    
        FT  px;
        FT  py;
        FT  qx;
        FT  qy;
        FT  rx;
        FT  ry;
    
        dao.get( p, px, py);
        dao.get( q, qx, qy);
        dao.get( r, rx, ry);
    
        return( static_cast< CGAL::Orientation>(
                  CGAL_NTS sign( ( px-rx) * ( qy-ry) - ( py-ry) * ( qx-rx))));
    }
};

// Nested type `Circle'
template < class PT_, class DA_ >
std::ostream&
operator << ( std::ostream&,
              const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

template < class PT_, class DA_ >
std::istream&
operator >> ( std::istream&,
              CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

template < class PT_, class DA_ >
class _Min_circle_2_adapterC2__Circle {
  public:
    // typedefs
    typedef  PT_  PT;
    typedef  DA_  DA;

    typedef  typename DA_::FT  FT;

  private:
    // data members
    DA  dao;                                // data accessor object

    FT  center_x;                           // center's x-coordinate
    FT  center_y;                           // center's y-coordinate
    FT  sqr_rad;                            // squared radius

    // private member functions
    FT
    sqr_dist( const FT& px, const FT& py, const FT& qx, const FT& qy) const
    {
        FT  dx( px - qx);
        FT  dy( py - qy);
        return( dx*dx + dy*dy);
    }

    friend  std::ostream&  operator << <> ( std::ostream&,
        const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

    friend  std::istream&  operator >> <> ( std::istream&,
        CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

  public:
    // types
    typedef  PT  Point;
    typedef  FT  Distance;

    // creation
    _Min_circle_2_adapterC2__Circle( const DA& da) : dao( da) { }

    void  set( )
    {
        center_x =  FT( 0);
        center_y =  FT( 0);
        sqr_rad  = -FT( 1);
    }

    void  set( const Point& p)
    {
        dao.get( p, center_x, center_y);
        sqr_rad = FT( 0);
    }

    void  set( const Point& p, const Point& q)
    {
        FT  px;
        FT  py;
        FT  qx;
        FT  qy;

        dao.get( p, px, py);
        dao.get( q, qx, qy);

        center_x = ( px+qx) / FT( 2);
        center_y = ( py+qy) / FT( 2);
        sqr_rad  = sqr_dist( px, py, center_x, center_y);
    }

    void  set( const Point& p, const Point& q, const Point& r)
    {
        FT  px;
        FT  py;
        FT  qx;
        FT  qy;
        FT  rx;
        FT  ry;

        dao.get( p, px, py);
        dao.get( q, qx, qy);
        dao.get( r, rx, ry);

        FT  qx_px( qx - px);
        FT  qy_py( qy - py);
        FT  rx_px( rx - px);
        FT  ry_py( ry - py);

        FT  p2   ( px*px + py*py);
        FT  q2_p2( qx*qx + qy*qy - p2);
        FT  r2_p2( rx*rx + ry*ry - p2);
        FT  denom( ( qx_px*ry_py - rx_px*qy_py) * FT( 2));

        center_x = ( q2_p2*ry_py - r2_p2*qy_py) / denom;
        center_y = ( r2_p2*qx_px - q2_p2*rx_px) / denom;
        sqr_rad  = sqr_dist( px, py, center_x, center_y);
    }

    // predicates
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        FT  px;
        FT  py;
        dao.get( p, px, py);
        return( CGAL::Bounded_side(
         CGAL_NTS sign( sqr_rad - sqr_dist( px, py, center_x, center_y))));
    }

    bool
    has_on_bounded_side( const Point& p) const
    {
        FT  px;
        FT  py;
        dao.get( p, px, py);
        return( sqr_dist( px, py, center_x, center_y) < sqr_rad);
    }

    bool
    has_on_boundary( const Point& p) const
    {
        FT  px;
        FT  py;
        dao.get( p, px, py);
        return( sqr_dist( px, py, center_x, center_y) == sqr_rad);
    }

    bool
    has_on_unbounded_side( const Point& p) const
    {
        FT  px;
        FT  py;
        dao.get( p, px, py);
        return( sqr_rad < sqr_dist( px, py, center_x, center_y));
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
        const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c) const
    {
        return( ( center_x == c.center_x) &&
                ( center_y == c.center_y) &&
                ( sqr_rad  == c.sqr_rad ) );
    }

    bool
    operator != (
        const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c) const
    {
        return( ! ( *this == c));
    }

    Point
    center( ) const
    {
        Point  p;
        dao.set( p, center_x, center_y);
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
              const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c)
{
    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << "CGAL::Min_circle_2_adapterC2::Circle( "
           << c.center_x << ", "
           << c.center_y << ", "
           << c.sqr_rad  << ')';
        break;

      case CGAL::IO::ASCII:
        os << c.center_x << ' ' << c.center_y << ' ' << c.sqr_rad;
        break;

      case CGAL::IO::BINARY:
        CGAL::write( os, c.center_x);
        CGAL::write( os, c.center_y);
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
              CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c)
{
    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
	std::cerr << std::endl;
	std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII:
        is >> c.center_x >> c.center_y >> c.sqr_rad;
        break;

      case CGAL::IO::BINARY:
        CGAL::read( is, c.center_x);
        CGAL::read( is, c.center_y);
        CGAL::read( is, c.sqr_rad);
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

} //namespace CGAL

#endif // CGAL_MIN_CIRCLE_2_ADAPTERC2_H

// ===== EOF ==================================================================
