// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : extremal_polygon_2_demo.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Demo program: Compute extremal polygons of a convex polygon
// ============================================================================


#include <CGAL/basic.h>

#if (defined(CGAL_USE_LEDA) || defined(CGAL_USE_CGAL_WINDOW))

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/circulator.h>
#include <CGAL/extremal_polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Window_stream.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

using CGAL::cgalize;
using CGAL::Timer;
using CGAL::has_smaller_dist_to_point;
using CGAL::RED;
using CGAL::BLUE;
using std::back_inserter;
using std::copy;
using std::max;
typedef CGAL::Window_stream Window;

using std::vector;
using CGAL::Cartesian;
using CGAL::Random_access_circulator_from_container;
using CGAL::convex_hull_points_2;
using CGAL::squared_distance;
using CGAL::maximum_area_inscribed_k_gon_2;
using CGAL::maximum_perimeter_inscribed_k_gon_2;
using CGAL::Random_points_in_square_2;
using CGAL::Creator_uniform_2;
using CGAL::random_convex_set_2;

typedef double                                    FT;
typedef Cartesian< FT >                           K;
typedef K::Point_2                                Point;
typedef K::Segment_2                              Segment;
typedef vector< Point >                           PointCont;
typedef PointCont::iterator                       PointIter;
typedef vector< PointIter >                       PointIterCont;
typedef PointIterCont::iterator                   PointIterIter;
typedef PointCont                                 Polygon;
typedef PointCont::iterator                       Vertex_iterator;
typedef Random_access_circulator_from_container< Polygon >
  Vertex_circulator;
typedef PointCont::const_iterator                 Vertex_const_iterator;
typedef vector< Vertex_iterator >                 Vertex_iteratorCont;
typedef Vertex_iteratorCont::iterator             Vertex_iteratorIter;


void
wait_for_button_release( Window& W)
{
  // wait until mouse button is released
  double x, y;
  int v;
  do {}
#ifdef CGAL_USE_LEDA
  while ( W.read_event( v, x, y) != button_release_event);
#else
  while ( W.read_event( v, x, y) != CGAL::button_release_event);
#endif // CGAL_USE_LEDA
}

int
main()
{
  Window W( 650, 650);
  int k( 3);
  W.int_item( "k", k, 2, 12, "#vertices of polygon to inscribe");
  int n( 3);
  W.int_item( "n", n, "#vertices of polygon");
  int compute_mode( 0);
  W.choice_item( "criterion",
                 compute_mode,
                 "Determine the criterion to maximize.",
                 3,
                 "Area",
                 "Perimeter",
                 "Both");
  int compute_button(
    W.button( "Compute",
              "Compute largest inscribed k-gon"));
  int generate_button(
    W.button( "Generate",
              "Generate random convex polygon"));
  int help_button(
    W.button( "Help",
              "Explain the program and its mouse interaction."));
  int quit_button(
    W.button( "Quit",
              "Leave the Program"));
  cgalize( W);
  W.display();
  W.init( -1.5, 1.5, -1.5);

  PointCont points;
  Polygon p;
  bool polygon_changed(false);
  bool done(false);

  while(!done) {
    if (polygon_changed) {
      // compute convex hull:
      PointCont ch_points;
      convex_hull_points_2(points.begin(),
                           points.end(),
                           back_inserter(ch_points));
    
      // replace points by ch_points:
      points.erase(points.begin(), points.end());
      points.insert(points.begin(), ch_points.begin(), ch_points.end());
    
      // construct polygon and data structure for point location:
      p.erase(p.begin(), p.end());
      copy(points.begin(), points.end(), back_inserter(p));
    
      polygon_changed = false;
    } // if ( polygon_changed)
    
    // show polygon:
    W.clear();
    if ( !p.empty()) {
      Vertex_const_iterator i( p.begin());
      for (;;) {
        W << RED << *i << BLUE;
    #ifndef _MSC_VER
        if ( (i+1) == p.end())
    #else
        if ( (i+1) == Vertex_const_iterator(p.end()))
    #endif // _MSC_VER
          {
            W << Segment( *i, *(p.begin()));
            break;
          }
        else {
          W << Segment( *i, *(i+1));
          ++i;
        }
      } // for (;;)
    } // if ( !p.empty())
    
    #ifndef _MSC_VER
    char vertices_message[80];
    int num( points.size());
    std::sprintf( vertices_message,
                  "Polygon has %d vertices.",
                  num);
    W.message( vertices_message);
    #endif // ! _MSC_VER
    
    

    if (compute_mode != 1) {
      // compute maximum area inscribed k-gon:
      k = max(3, k);
      PointCont k_gon;
      if (p.size() >= 3) {
#ifdef EXTREMAL_POLYGON_MEASURE
        Timer t;
        t.start();
#endif // EXTREMAL_POLYGON_MEASURE
        maximum_area_inscribed_k_gon_2(
          p.begin(),
          p.end(),
          k,
          back_inserter(k_gon));
#ifdef EXTREMAL_POLYGON_MEASURE
        t.stop();
        cout << "[time: " << t.time() << " msec]" << endl;
#endif // EXTREMAL_POLYGON_MEASURE
        W << CGAL::GREEN;
        if ( !p.empty()) {
          PointIter i( k_gon.begin());
          while ( ++i != k_gon.end())
            W << Segment( *(i-1), *i);
          W << Segment( *--i, *(k_gon.begin()));
        } // if ( !p.empty())
      } // if ( p.size() >= 3)
    } // if ( compute_mode != 1)
    if (compute_mode != 0) {
      // compute maximum perimeter inscribed k-gon:
      PointCont k_gon;
#ifdef EXTREMAL_POLYGON_MEASURE
      Timer t;
      t.start();
#endif // EXTREMAL_POLYGON_MEASURE
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
      maximum_perimeter_inscribed_k_gon_2(
        Vertex_circulator( &p),
        Vertex_circulator( &p),
        k,
        back_inserter( k_gon));
#else
      maximum_perimeter_inscribed_k_gon_2(
        p.begin(),
        p.end(),
        k,
        back_inserter( k_gon));
#endif
#ifdef EXTREMAL_POLYGON_MEASURE
      t.stop();
      cout << "[time: " << t.time() << " msec]" << endl;
#endif // EXTREMAL_POLYGON_MEASURE
      W << CGAL::ORANGE;
      if ( !p.empty()) {
        PointIter i( k_gon.begin());
        while ( ++i != k_gon.end())
          W << Segment( *(i-1), *i);
        W << Segment( *--i, *(k_gon.begin()));
      } // if ( !p.empty())
    } // if ( compute_mode != 0)

    // wait for input:
    int input;
    for (;;) {
      double x( 0), y( 0);
      input = W.get_mouse( x, y);
      if ( input == MOUSE_BUTTON( 2)) {
        // move point
        
        // find nearest vertex
        Point sp(x, y);
        PointIter nearest = points.begin();
        PointIter k = nearest;
        while (++k != points.end())
          if (has_smaller_dist_to_point(sp, *k, *nearest))
            nearest = k;
        
        // test for snapping:
        if ( squared_distance( *nearest, Point( x, y)) < FT( 0.05)) {
          int v( 0);
        #ifdef CGAL_USE_LEDA
          leda_drawing_mode old_drawing_mode( W.set_mode( leda_xor_mode));
          leda_color old_color( W.set_color( leda_grey1));
        #else
          CGAL::drawing_mode old_drawing_mode = W.set_mode( CGAL::xor_mode);
          CGAL::color old_color = W.set_color( CGAL::grey1);
        #endif // CGAL_USE_LEDA
          x = (*nearest).x();
          y = (*nearest).y();
          do {
            if ( nearest == points.begin()) {
              W.draw_segment( (*(points.end() - 1)).x(),
                              (*(points.end() - 1)).y(),
                              (*nearest).x(),
                              (*nearest).y());
              if ( points.size() > 2)
                W.draw_segment( (*(points.end() - 1)).x(),
                                (*(points.end() - 1)).y(),
                                x,
                                y);
            }
            else {
              W.draw_segment( (*(nearest - 1)).x(),
                              (*(nearest - 1)).y(),
                              (*nearest).x(),
                              (*nearest).y());
              if ( points.size() > 2)
                W.draw_segment( (*(nearest - 1)).x(),
                                (*(nearest - 1)).y(),
                                x,
                                y);
            }
            if ( nearest == (points.end() - 1)) {
              W.draw_segment( (*(points.begin())).x(),
                              (*(points.begin())).y(),
                              (*nearest).x(),
                              (*nearest).y());
              if ( points.size() > 2)
                W.draw_segment( (*(points.begin())).x(),
                                (*(points.begin())).y(),
                                x,
                                y);
            }
            else {
              W.draw_segment( (*(nearest + 1)).x(),
                              (*(nearest + 1)).y(),
                              (*nearest).x(),
                              (*nearest).y());
              if ( points.size() > 2)
                W.draw_segment( (*(nearest + 1)).x(),
                                (*(nearest + 1)).y(),
                                x,
                                y);
            }
        
            // set new point:
            *nearest = Point( x, y);
        
          }
        #ifdef CGAL_USE_LEDA
          while ( W.read_event( v, x, y) != button_release_event);
        #else
          while ( W.read_event( v, x, y) != CGAL::button_release_event);
        #endif // CGAL_USE_LEDA
        
          // restore parameters of W:
          W.set_mode( old_drawing_mode);
          W.set_color( old_color);
        
          polygon_changed = true;
          break;
        }
      }
      else if ( input == MOUSE_BUTTON( 3)) {
        // delete point
        
        // find nearest vertex
        Point sp(x, y);
        PointIter nearest = points.begin();
        PointIter k = nearest;
        while (++k != points.end())
          if (has_smaller_dist_to_point(sp, *k, *nearest))
            nearest = k;
        
        // test for snapping:
        if (squared_distance(*nearest, Point(x, y)) < FT(0.05)) {
          points.erase(nearest);
          polygon_changed = true;
          break;
        }    
      }
      else if ( input == help_button) {
        // display help text
        W.del_messages();
        W.message( "CGAL MONOTONE MATRIX SEARCH TEST");
        W.message( "");
        W.message( "compute maximal inscribed k-gon of a convex polygon");
        W.message( "(maximal with respect to either area or perimeter)");
    
        W.message( "");
        W.message( "");
        W.message( "Mouse Input:");
        W.message( "");
        W.message( " <left mouse button>      -  insert vertex");
        W.message( " <middle mouse button> -  move vertex");
        W.message( " <right mouse button>     -  delete vertex");
      }
      else if ( input == quit_button) {
        // quit program
        done = true;
        break;
      } else if ( input == generate_button) {
        // generate random convex polygon with n vertices
        typedef Random_points_in_square_2<
          Point,
          Creator_uniform_2< FT, Point >
        >
        Point_generator;
        
        points.erase( points.begin(), points.end());
        n = max( n, 3);
        random_convex_set_2( n, back_inserter( points), Point_generator( 1));
        
        polygon_changed = true;
        break;
      }
      else if ( input == compute_button)
        // recompute largest inscribing k-gon
        break;
      else if ( input == MOUSE_BUTTON( 1)) {
        // insert point:
        points.push_back( Point( x, y));
        polygon_changed = true;
        break;
      }
    } // for (;;)
  } // while (!done)

  return 0;
} // int main()



#else

#include <iostream>

int main()
{
  std::cerr << "This demo requires LEDA." << std::endl;
  return 0;
}

#endif

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

