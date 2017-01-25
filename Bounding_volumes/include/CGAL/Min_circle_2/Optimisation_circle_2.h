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
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

#ifndef CGAL_OPTIMISATION_CIRCLE_2_H
#define CGAL_OPTIMISATION_CIRCLE_2_H

#include <CGAL/license/Bounding_volumes.h>


// includes

#  include <CGAL/basic_constructions_2.h>
#  include <CGAL/squared_distance_2.h>

namespace CGAL {

// Class declaration
// =================
template < class K_ >
class Optimisation_circle_2;

// Class interface
// ===============
template < class K_ >
class Optimisation_circle_2 {
  public:
    // types
    typedef           K_                K;
    typedef typename  K::Point_2        Point;
    typedef typename  K_::FT            Distance;
    
    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.
    
    // creation
    Optimisation_circle_2( );
    
    void  set( );
    void  set( const Point& p);
    void  set( const Point& p, const Point& q);
    void  set( const Point& p, const Point& q, const Point& r);
    void  set( const Point& center, const Distance& squared_radius);
    
    // access functions
    const Point&     center        ( ) const;
    const Distance&  squared_radius( ) const
    
    // equality tests
    bool  operator == ( const Optimisation_circle_2<K>& c) const;
    bool  operator != ( const Optimisation_circle_2<K>& c) const;
    
    // predicates
    CGAL::Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;
    
    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
    **************************************************************************/

  private:
    // private data members
    Point     _center;
    Distance  _squared_radius;

// ============================================================================

// Class implementation
// ====================

  public:
    // Constructor
    // -----------
    inline
    Optimisation_circle_2( )
    { }

    // Set functions
    // -------------
    inline
    void
    set( )
    {
        _center         =  Point( CGAL::ORIGIN);
        _squared_radius = -Distance( 1);
    }
    
    inline
    void
    set( const Point& p)
    {
        _center         = p;
        _squared_radius = Distance( 0);
    }
    
    inline
    void
    set( const Point& p, const Point& q)
    {
        _center         = CGAL::midpoint( p, q);
        _squared_radius = CGAL::squared_distance( p, _center);
    }
    
    inline
    void
    set( const Point& p, const Point& q, const Point& r)
    {
        _center         = CGAL::circumcenter( p, q, r);
        _squared_radius = CGAL::squared_distance( p, _center);
    }
    
    inline
    void
    set( const Point& center, const Distance& squared_radius)
    {
        _center         = center;
        _squared_radius = squared_radius;
    }

    // Access functions
    // ----------------
    inline
    const Point&
    center( ) const
    {
        return( _center);
    }
    
    inline
    const Distance&
    squared_radius( ) const
    {
        return( _squared_radius);
    }

    // Equality tests
    // --------------
    bool
    operator == ( const Optimisation_circle_2<K>& c) const
    {
        return( ( _center          == c._center        ) &&
                ( _squared_radius  == c._squared_radius) );
    }
    
    bool
    operator != ( const Optimisation_circle_2<K>& c) const
    {
        return( ! operator==( c));
    }

    // Predicates
    // ----------
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( CGAL::Bounded_side( CGAL_NTS sign(
            _squared_radius - CGAL::squared_distance( p, _center))));
    }
    
    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( CGAL::squared_distance( p, _center) < _squared_radius);
    }
    
    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( CGAL::squared_distance( p, _center) == _squared_radius);
    }
    
    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( _squared_radius < CGAL::squared_distance( p, _center));
    }
    
    inline
    bool
    is_empty( ) const
    {
        return( CGAL::is_negative( _squared_radius));
    }
    
    inline
    bool
    is_degenerate( ) const
    {
        return( ! CGAL::is_positive( _squared_radius));
    }
};

// Function declarations
// =====================
// I/O
// ---
template < class K_ >
std::ostream&
operator << ( std::ostream&, const CGAL::Optimisation_circle_2<K_>&);

template < class K_ >
std::istream&
operator >> ( std::istream&, CGAL::Optimisation_circle_2<K_>&);

} //namespace CGAL

#include <CGAL/Min_circle_2/Optimisation_circle_2_impl.h>

#endif // CGAL_OPTIMISATION_CIRCLE_2_H

// ===== EOF ==================================================================
