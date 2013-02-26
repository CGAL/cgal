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

#ifndef CGAL_MIN_ELLIPSE_2_H
#define CGAL_MIN_ELLIPSE_2_H

#include <CGAL/Optimisation/basic.h>
#include <CGAL/Random.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>

namespace CGAL {

// Class declaration
// =================
template < class Traits_ >
class Min_ellipse_2;

// Class interface
// ===============
template < class Traits_ >
class Min_ellipse_2 {
  public:
    // types
    typedef           Traits_                           Traits;
    typedef typename  Traits_::Point                    Point;
    typedef typename  Traits_::Ellipse                  Ellipse;
    typedef typename  std::list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                     Support_point_iterator;
    
    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.
    
    // creation
    template < class InputIterator >
    Min_ellipse_2( InputIterator first,
                   InputIterator last,
                   bool          randomize = false,
                   Random&       random    = default_random,
                   const Traits& traits    = Traits());
    
    Min_ellipse_2( const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p,
                   const Traits& traits = Traits());
    Min_ellipse_2( Point  p,
                   Point  q,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p1,
                   const Point&  p2,
                   const Point&  p3,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p1,
                   const Point&  p2,
                   const Point&  p3,
                   const Point&  p4,
                   const Traits& traits = Traits());
    Min_ellipse_2( const Point&  p1,
                   const Point&  p2,
                   const Point&  p3,
                   const Point&  p4,
                   const Point&  p5,
                   const Traits& traits = Traits());
    ~Min_ellipse_2( );
    
    // access functions
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;
    
    Point_iterator  points_begin( ) const;
    Point_iterator  points_end  ( ) const;
    
    Support_point_iterator  support_points_begin( ) const;
    Support_point_iterator  support_points_end  ( ) const;
    
    const Point&  support_point( int i) const;
    
    const Ellipse&  ellipse( ) const;
    
    // predicates
    CGAL::Bounded_side  bounded_side( const Point& p) const;
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
    int          n_support_points;              // number of support points
    Point*       support_points;                // array of support points
    

    // copying and assignment not allowed!
    Min_ellipse_2( const Min_ellipse_2<Traits_>&);
    Min_ellipse_2<Traits_>& operator = ( const Min_ellipse_2<Traits_>&);

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
        return( number_of_support_points() <  3);
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
    // ellipse
    inline
    const Ellipse&
    ellipse( ) const
    {
        return( tco.ellipse);
    }
    

    // in-ellipse test predicates
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( tco.ellipse.bounded_side( p));
    }
    
    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( tco.ellipse.has_on_bounded_side( p));
    }
    
    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( tco.ellipse.has_on_boundary( p));
    }
    
    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( tco.ellipse.has_on_unbounded_side( p));
    }

  private:
    // Private member functions
    // ------------------------
    // compute_ellipse
    inline
    void
    compute_ellipse( )
    {
        switch ( n_support_points) {
          case 5:
            tco.ellipse.set( support_points[ 0],
                             support_points[ 1],
                             support_points[ 2],
                             support_points[ 3],
                             support_points[ 4]);
            break;
          case 4:
            tco.ellipse.set( support_points[ 0],
                             support_points[ 1],
                             support_points[ 2],
                             support_points[ 3]);
            break;
          case 3:
            tco.ellipse.set( support_points[ 0],
                             support_points[ 1],
                             support_points[ 2]);
            break;
          case 2:
            tco.ellipse.set( support_points[ 0], support_points[ 1]);
            break;
          case 1:
            tco.ellipse.set( support_points[ 0]);
            break;
          case 0:
            tco.ellipse.set( );
            break;
          default:
            CGAL_optimisation_assertion( ( n_support_points >= 0) &&
                                         ( n_support_points <= 5) ); }
    }

    void
    me( const Point_iterator& last, int n_sp)
    {
        // compute ellipse through support points
        n_support_points = n_sp;
        compute_ellipse();
        if ( n_sp == 5) return;
    
        // test first n points
        typename std::list<Point>::iterator  point_iter = points.begin();
        for ( ; last != point_iter; ) {
            const Point& p = *point_iter;
    
            // p not in current ellipse?
            if ( has_on_unbounded_side( p)) {
    
                // recursive call with p as additional support point
                support_points[ n_sp] = p;
                me( point_iter, n_sp+1);
    
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
    Min_ellipse_2( InputIterator first,
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
            support_points = new Point[ 5];
    
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
    
            // compute me
            me( points.end(), 0);
        }
    
    // default constructor
    inline
    Min_ellipse_2( const Traits& traits = Traits())
        : tco( traits), n_support_points( 0)
    {
        // allocate support points' array
        support_points = new Point[ 5];
    
        // initialize ellipse
        tco.ellipse.set();
    
        CGAL_optimisation_postcondition( is_empty());
    }
    
    // constructor for one point
    inline
    Min_ellipse_2( const Point& p, const Traits& traits = Traits())
        : tco( traits), points( 1, p), n_support_points( 1)
    {
        // allocate support points' array
        support_points = new Point[ 5];
    
        // initialize ellipse
        support_points[ 0] = p;
        tco.ellipse.set( p);
    
        CGAL_optimisation_postcondition( is_degenerate());
    }
    
    // constructor for two points
    // This was const Point& but then Intel 7.0/.net2003 messes it up 
    // with the constructor taking an iterator range
    inline
    Min_ellipse_2( Point p1, Point p2,
                   const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 5];
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
    
        // compute me
        me( points.end(), 0);
    
        CGAL_optimisation_postcondition( is_degenerate());
    }
    
    // constructor for three points
    inline
    Min_ellipse_2( const Point& p1, const Point& p2, const Point& p3,
                   const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 5];
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
    
        // compute me
        me( points.end(), 0);
    }
    
    // constructor for four points
    inline
    Min_ellipse_2( const Point& p1, const Point& p2,
                   const Point& p3, const Point& p4,
                   const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 5];
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
        points.push_back( p4);
    
        // compute me
        me( points.end(), 0);
    }
    
    // constructor for five points
    inline
    Min_ellipse_2( const Point& p1, const Point& p2, const Point& p3,
                   const Point& p4, const Point& p5,
                   const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 5];
    
        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);
        points.push_back( p4);
        points.push_back( p5);
    
        // compute me
        me( points.end(), 0);
    }
    

    // Destructor
    // ----------
    inline
    ~Min_ellipse_2( )
    {
        // free support points' array
        delete[] support_points;
    }

    // Modifiers
    // ---------
    void
    insert( const Point& p)
    {
        // p not in current ellipse?
        if ( has_on_unbounded_side( p)) {
    
            // p new support point
            support_points[ 0] = p;
    
            // recompute me
            me( points.end(), 1);
    
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
        tco.ellipse.set();
    }
    

    // Validity check
    // --------------
    bool
    is_valid( bool verbose = false, int level = 0) const
    {
        using namespace std;
    
        CGAL::Verbose_ostream verr( verbose);
        verr << endl;
        verr << "CGAL::Min_ellipse_2<Traits>::" << endl;
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
                            "ellipse does not contain all points"));
        verr << "passed." << endl;
    
        // support set checks (b)+(c) (not yet implemented)
        
        // alternative support set check
        verr << "  +) support set check..." << flush;
        Support_point_iterator support_point_iter;
        for ( support_point_iter  = support_points_begin();
              support_point_iter != support_points_end();
              ++support_point_iter)
            if ( ! has_on_boundary( *support_point_iter))
                return( CGAL::_optimisation_is_valid_fail( verr,
                            "ellipse does not have all \
                             support points on the boundary"));
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
operator << ( std::ostream& os, const Min_ellipse_2<Traits_>& me);

template < class Traits_ >
std::istream&
operator >> ( std::istream& is,       Min_ellipse_2<Traits_>& me);

} //namespace CGAL

#include <CGAL/Min_ellipse_2/Min_ellipse_2_impl.h>

#endif // CGAL_MIN_ELLIPSE_2_H

// ===== EOF ==================================================================
