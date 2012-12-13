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

#ifndef CGAL_MIN_CIRCLE_2_H
#define CGAL_MIN_CIRCLE_2_H

// includes
#  include <CGAL/Optimisation/basic.h>
#  include <CGAL/Random.h>
#  include <list>
#  include <vector>
#  include <algorithm>
#  include <iostream>

namespace CGAL {

// Class declaration
// =================
template < class Traits_ >
class Min_circle_2;

// Class interface
// ===============
template < class Traits_ >
class Min_circle_2 {
  public:
    // types
    typedef           Traits_                           Traits;
    typedef typename  Traits_::Point                    Point;
    typedef typename  Traits_::Circle                   Circle;
    typedef typename  std::list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                     Support_point_iterator;
    
    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.
    
    // creation
    template < class InputIterator >
    Min_circle_2( InputIterator first,
                  InputIterator last,
                  bool          randomize = false,
                  Random&       random    = default_random,
                  const Traits& traits    = Traits());
    
    Min_circle_2( const Traits& traits = Traits());
    Min_circle_2( const Point&  p,
                  const Traits& traits = Traits());
    Min_circle_2( Point  p, Point  q,
                  const Traits& traits = Traits());
    Min_circle_2( const Point&  p1, const Point&  p2, const Point&  p3,
                  const Traits& traits = Traits());
    ~Min_circle_2( );
    
    // access functions
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;
    
    Point_iterator  points_begin( ) const;
    Point_iterator  points_end  ( ) const;
    
    Support_point_iterator  support_points_begin( ) const;
    Support_point_iterator  support_points_end  ( ) const;
    
    const Point&  support_point( int i) const;
    
    const Circle&  circle( ) const;
    
    // predicates
    Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;
    
    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
    
    // modifiers
    void  insert( const Point& p);
    void  insert( const Point* first,
                  const Point* last );
    void  insert( std::list<Point>::const_iterator first,
                  std::list<Point>::const_iterator last );
    void  insert( std::istream_iterator<Point,std::ptrdiff_t> first,
                  std::istream_iterator<Point,std::ptrdiff_t> last );
    void  clear( );
    
    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;
    
    // miscellaneous
    const Traits&  traits( ) const;
    **************************************************************************/

  private:
    // private data members
    Traits       tco;                           // traits class object
    std::list<Point>  points;                   // doubly linked list of points
  std::size_t         n_support_points;              // number of support points
    Point*       support_points;                // array of support points
    

    // copying and assignment not allowed!
    Min_circle_2( const Min_circle_2<Traits_>&);
    Min_circle_2<Traits_>&
        operator = ( const Min_circle_2<Traits_>&);

// ============================================================================

// Class implementation
// ====================

  public:
    // Access functions and predicates
    // -------------------------------
    // #points and #support points
    inline
    std::size_t
    number_of_points( ) const
    {
        return( points.size());
    }
    
    inline
    std::size_t
    number_of_support_points( ) const
    {
        return( n_support_points);
    }

    // is_... predicates
    inline
    bool
    is_empty( ) const
    {
        return( number_of_support_points() == 0);
    }
    
    inline
    bool
    is_degenerate( ) const
    {
        return( number_of_support_points() <  2);
    }

    // access to points and support points
    inline
    Point_iterator
    points_begin( ) const
    {
        return( points.begin());
    }
    
    inline
    Point_iterator
    points_end( ) const
    {
        return( points.end());
    }
    
    inline
    Support_point_iterator
    support_points_begin( ) const
    {
        return( support_points);
    }
    
    inline
    Support_point_iterator
    support_points_end( ) const
    {
        return( support_points+n_support_points);
    }
    
    // random access for support points
    inline
    const Point&
    support_point( std::size_t i) const
    {
        CGAL_optimisation_precondition(i <  number_of_support_points());
        return( support_points[ i]);
    }
    // circle
    inline
    const Circle&
    circle( ) const
    {
        return( tco.circle);
    }
    

    // in-circle test predicates
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( tco.circle.bounded_side( p));
    }
    
    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( tco.circle.has_on_bounded_side( p));
    }
    
    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( tco.circle.has_on_boundary( p));
    }
    
    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( tco.circle.has_on_unbounded_side( p));
    }

  private:
    // Private member functions
    // ------------------------
    // compute_circle
    inline
    void
    compute_circle( )
    {
        switch ( n_support_points) {
          case 3:
            tco.circle.set( support_points[ 0],
                            support_points[ 1],
                            support_points[ 2]);
            break;
          case 2:
            tco.circle.set( support_points[ 0], support_points[ 1]);
            break;
          case 1:
            tco.circle.set( support_points[ 0]);
            break;
          case 0:
            tco.circle.set( );
            break;
          default:
            CGAL_optimisation_assertion( n_support_points <= 3 ); }
    }

    void
    mc( const Point_iterator& last, std::size_t n_sp)
    {
        // compute circle through support points
        n_support_points = n_sp;
        compute_circle();
        if ( n_sp == 3) return;
    
        // test first n points
        typename std::list<Point>::iterator  point_iter = points.begin();
        for ( ; last != point_iter; ) {
            const Point& p = *point_iter;
    
            // p not in current circle?
            if ( has_on_unbounded_side( p)) {
    
                // recursive call with p as additional support point
                support_points[ n_sp] = p;
                mc( point_iter, n_sp+1);
    
                // move current point to front
                points.splice( points.begin(), points, point_iter++); }
    
            else
                ++point_iter; }
    }

  public:
    // Constructors
    // ------------
    // STL-like constructor (member template)
    template < class InputIterator >
    Min_circle_2( InputIterator first,
                  InputIterator last,
                  bool          randomize
    #if !defined(_MSC_VER) || _MSC_VER > 1300
                                              = false
    #endif
                                                     ,
                      Random&       random    = default_random,
                      const Traits& traits    = Traits())
            : tco( traits)
        {
            // allocate support points' array
            support_points = new Point[ 3];
    
            // range of points not empty?
            if ( first != last) {
    
                // store points
                if ( randomize) {
    
                    // shuffle points at random
                    std::vector<Point> v( first, last);
                    std::random_shuffle( v.begin(), v.end(), random);
                    std::copy( v.begin(), v.end(),
                               std::back_inserter( points)); }
                else
                    std::copy( first, last, std::back_inserter( points)); }
    
            // compute mc
            mc( points.end(), 0);
        }
    
    // default constructor
    inline
    Min_circle_2( const Traits& traits = Traits())
        : tco( traits), n_support_points( 0)
    {
        // allocate support points' array
        support_points = new Point[ 3];
    
        // initialize circle
        tco.circle.set();
    
        CGAL_optimisation_postcondition( is_empty());
    }
    
    // constructor for one point
    inline
    Min_circle_2( const Point& p, const Traits& traits = Traits())
        : tco( traits), points( 1, p), n_support_points( 1)
    {
        // allocate support points' array
        support_points = new Point[ 3];
    
        // initialize circle
        support_points[ 0] = p;
        tco.circle.set( p);
    
        CGAL_optimisation_postcondition( is_degenerate());
    }
    
    // constructor for two points
    // This was const Point& but then Intel 7.0/.net2003 messes it up 
    // with the constructor taking an iterator range
    inline
    Min_circle_2( Point p1, Point p2,
                  const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 3];
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
    
        // compute mc
        mc( points.end(), 0);
    }
    
    // constructor for three points
    inline
    Min_circle_2( const Point& p1, const Point& p2, const Point& p3,
                  const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 3];
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
    
        // compute mc
        mc( points.end(), 0);
    }
    

    // Destructor
    // ----------
    inline
    ~Min_circle_2( )
    {
        // free support points' array
        delete[] support_points;
    }

    // Modifiers
    // ---------
    void
    insert( const Point& p)
    {
        // p not in current circle?
        if ( has_on_unbounded_side( p)) {
    
            // p new support point
            support_points[ 0] = p;
    
            // recompute mc
            mc( points.end(), 1);
    
            // store p as the first point in list
            points.push_front( p); }
        else
    
            // append p to the end of the list
            points.push_back( p);
    }
    
        template < class InputIterator >
        void
        insert( InputIterator first, InputIterator last)
        {
            for ( ; first != last; ++first)
                insert( *first);
        }
    
    void
    clear( )
    {
        points.erase( points.begin(), points.end());
        n_support_points = 0;
        tco.circle.set();
    }
    

    // Validity check
    // --------------
    bool
    is_valid( bool verbose = false, int level = 0) const
    {
        using namespace std;
    
        CGAL::Verbose_ostream verr( verbose);
        verr << endl;
        verr << "CGAL::Min_circle_2<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P| = " << number_of_points()
             << ", |S| = " << number_of_support_points() << endl;
    
        // containment check (a)
        verr << "  a) containment check..." << flush;
        Point_iterator point_iter;
        for ( point_iter  = points_begin();
              point_iter != points_end();
              ++point_iter)
            if ( has_on_unbounded_side( *point_iter))
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "circle does not contain all points"));
        verr << "passed." << endl;
    
        // support set checks (b)+(c)
        verr << "  b)+c) support set checks..." << flush;
        switch( number_of_support_points()) {
        
          case 0:
            if ( ! is_empty())
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "P is nonempty, \
                             but there are no support points."));
            break;
        
          case 1:
            if ( ( circle().center() != support_point( 0)    ) ||
                 ( ! CGAL_NTS is_zero( circle().squared_radius())) )
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "circle differs from the circle \
                             spanned by its single support point."));
            break;
        
          case 2: {
            const Point& p = support_point( 0),
                         q = support_point( 1);
            
            // p equals q?
            if ( p == q)
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "the two support points are equal."));
            
            // segment(p,q) is not diameter?
            if ( ( ! has_on_boundary( p)                                ) ||
                 ( ! has_on_boundary( q)                                ) ||
                 ( tco.orientation( p, q,
                                    circle().center()) != CGAL::COLLINEAR) )
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "circle does not have its \
                             two support points as diameter.")); }
            break;
        
          case 3: {
            const Point& p = support_point( 0),
                         q = support_point( 1),
                         r = support_point( 2);
            
            // p, q, r not pairwise distinct?
            if ( ( p == q) || ( q == r) || ( r == p))
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "at least two of the three \
                             support points are equal."));
            
            // p, q, r collinear?
            if ( tco.orientation( p, q, r) == CGAL::COLLINEAR)
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "the three support points are collinear."));
            
            // current circle not equal the unique circle through p,q,r ?
            Circle c = circle();
            c.set( p, q, r);
            if ( circle() != c)
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "circle is not the unique circle \
                             through its three support points."));
            
            // circle's center on boundary of triangle(p,q,r)?
            const Point& center = circle().center();
            CGAL::Orientation o_pqz = tco.orientation( p, q, center);
            CGAL::Orientation o_qrz = tco.orientation( q, r, center);
            CGAL::Orientation o_rpz = tco.orientation( r, p, center);
            if ( ( o_pqz == CGAL::COLLINEAR) ||
                 ( o_qrz == CGAL::COLLINEAR) ||
                 ( o_rpz == CGAL::COLLINEAR) )
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "one of the three support points is redundant."));
            
            // circle's center not inside triangle(p,q,r)?
            if ( ( o_pqz != o_qrz) || ( o_qrz != o_rpz) || ( o_rpz != o_pqz))
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "circle's center is not in the \
                             convex hull of its three support points.")); }
            break;
        
          default:
            return( CGAL::_optimisation_is_valid_fail( verr,
                        "illegal number of support points, \
                         not between 0 and 3."));
        };
        verr << "passed." << endl;
    
        verr << "  object is valid!" << endl;
        return( true);
    }

    // Miscellaneous
    // -------------
    inline
    const Traits&
    traits( ) const
    {
        return( tco);
    }
};

// Function declarations
// =====================
// I/O
// ---
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os, const Min_circle_2<Traits_>& mc);

template < class Traits_ >
std::istream&
operator >> ( std::istream& is,       Min_circle_2<Traits_>& mc);

} //namespace CGAL

#include <CGAL/Min_circle_2/Min_circle_2_impl.h>

#endif // CGAL_MIN_CIRCLE_2_H

// ===== EOF ==================================================================
