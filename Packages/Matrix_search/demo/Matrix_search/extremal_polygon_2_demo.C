// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Demo program: Compute extremal polygons of a convex polygon
// ============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/circulator.h>
#include <CGAL/copy_n.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <CGAL/IO/leda_window.h>
#include <LEDA/point_set.h>
#include <iostream.h>

typedef double                                    FT;
typedef CGAL_Cartesian< FT >                      RepClass;
typedef CGAL_Point_2< RepClass >                  Point_2;
typedef CGAL_Segment_2< RepClass >                Segment;
typedef vector< Point_2 >                         PointCont;
typedef PointCont::iterator                       PointIter;
typedef vector< PointIter >                       PointIterCont;
typedef vector< PointIter >::iterator             PointIterIter;
typedef RepClass::FT                              FT;
typedef PointCont                                 Polygon;
typedef PointCont::iterator                       Vertex_iterator;
typedef CGAL_Random_access_circulator_from_container<
  Polygon >                                       Vertex_circulator;
typedef PointCont::const_iterator                 Vertex_const_iterator;
typedef vector< Vertex_iterator >                 Vertex_iteratorCont;
typedef Vertex_iteratorCont::iterator             Vertex_iteratorIter;
#include <LEDA/REDEFINE_NAMES.h>
typedef point                                     LEDA_Point;
typedef point_set< PointIter >                    LEDA_Point_set_PointIter;
#include <LEDA/UNDEFINE_NAMES.h>

#ifdef CGAL_EXTREMAL_POLYGON_MEASURE
#ifndef CGAL_PROTECT_TIME_H
#include <time.h>
#define CGAL_PROTECT_TIME_H
#endif // CGAL_PROTECT_TIME_H
static time_t Measure;
static long long int measure;
#define MEASURE(comm) \
  Measure = clock(); \
  comm; \
  measure = \
    (long long int)((float)(clock() - Measure) * 1000 / CLOCKS_PER_SEC); \
    cout << "[time: " << measure << " msec]" << endl;
#else
#define MEASURE(comm) comm
#endif
void
wait_for_button_release( leda_window& W)
{
  // wait until mouse button is released
  double x, y;
  int v;
  do {}
  while ( W.read_event( v, x, y) != button_release_event);
}

int
main()
{
  leda_window W( 650, 650);
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
    W.button( "End",
              "Quit Program"));
  W.display();
  W.init( -1.5, 1.5, -1.5);

  PointCont points;
  Polygon p;
  bool polygon_changed( false);
  LEDA_Point_set_PointIter point_location;

  for (;;) {
    if ( polygon_changed) {
      // compute convex hull:
      PointCont ch_points;
      CGAL_convex_hull_points_2( points.begin(),
                                 points.end(),
                                 back_inserter( ch_points));
    
      // replace points by ch_points:
      points.erase( points.begin(), points.end());
      points.insert( points.begin(),
                     ch_points.begin(),
                     ch_points.end());
    
      // construct polygon and data structure for point location:
      p.erase( p.begin(), p.end());
      point_location.clear();
      for ( PointIter p1( points.begin());
            p1 != points.end();
            ++p1) {
        p.push_back( *p1);
        point_location.insert(
          LEDA_Point( (*p1).x(), (*p1).y()),
          p1);
      }
    
      polygon_changed = false;
    } // if ( polygon_changed)
    
    // show polygon:
    W.clear();
    if ( !p.empty()) {
      Vertex_const_iterator i( p.begin());
      for (;;) {
        W << CGAL_RED << *i;
        W.set_fg_color( leda_blue);
        if ( (i+1) == p.end()) {
          W << Segment( *i, *(p.begin()));
          break;
        }
        else {
          W << Segment( *i, *(i+1));
          ++i;
        }
      } // for (;;)
    } // if ( !p.empty())
    
    char vertices_message[80];
    int num( points.size());
    sprintf( vertices_message,
             "Polygon has %d vertices.",
             num);
    W.message( vertices_message);
    
    

    if ( compute_mode != 1) {
      // compute maximum area inscribed k-gon:
      k = max( 3, k);
      PointCont k_gon;
      if ( p.size() >= 3) {
        MEASURE(CGAL_maximum_area_inscribed_k_gon(
          p.begin(),
          p.end(),
          k,
          back_inserter( k_gon));)
        W.set_fg_color( leda_green);
        if ( !p.empty()) {
          PointIter i( k_gon.begin());
          while ( ++i != k_gon.end())
            W << Segment( *(i-1), *i);
          W << Segment( *--i, *(k_gon.begin()));
        } // if ( !p.empty())
      } // if ( p.size() >= 3)
    } // if ( compute_mode != 1)
    if ( compute_mode != 0) {
      // compute maximum perimeter inscribed k-gon:
      PointCont k_gon;
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
      MEASURE(CGAL_maximum_perimeter_inscribed_k_gon(
        Vertex_circulator( &p),
        Vertex_circulator( &p),
        k,
        back_inserter( k_gon));)
#else
      MEASURE(CGAL_maximum_perimeter_inscribed_k_gon(
        p.begin(),
        p.end(),
        k,
        back_inserter( k_gon));)
#endif
      W.set_fg_color( leda_orange);
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
        PointIter nearest(
          point_location.inf(
            point_location.nearest_neighbor(
              LEDA_Point( x, y))));
        // test for snapping:
        /*
        #ifdef __GNUG__
        FT di( CGAL_squared_distance( *nearest, Point_2( x, y)));
        if ( di < FT( 0.05)) {
        #else
        */
        if ( CGAL_squared_distance( *nearest, Point_2( x, y)) < FT( 0.05)) {
        // #endif
          int v( 0);
          drawing_mode old_drawing_mode( W.set_mode( leda_xor_mode));
          leda_color old_color( W.set_color( leda_grey1));
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
            *nearest = Point_2( x, y);
        
          } while ( W.read_event( v, x, y) != button_release_event);
        
        
          // restore parameters of W:
          W.set_mode( old_drawing_mode);
          W.set_color( old_color);
        
          polygon_changed = true;
          break;
        }
      }
      else if ( input == MOUSE_BUTTON( 3)) {
        // delete point
        PointIter nearest(
          point_location.inf(
            point_location.nearest_neighbor(
              LEDA_Point( x, y))));
        // test for snapping:
        /*
        #ifdef __GNUG__
        FT di( CGAL_squared_distance( *nearest, Point_2( x, y)));
        if ( di < FT( 0.05)) {
        #else
        */
        if ( CGAL_squared_distance( *nearest, Point_2( x, y)) < FT( 0.05)) {
        //#endif
          points.erase( nearest);
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
      else if ( input == quit_button)
        // quit program
        return 0;
      else if ( input == generate_button) {
        // generate random convex polygon with n vertices
        typedef CGAL_Random_points_in_square_2<
          Point_2,
          CGAL_Creator_uniform_2< FT, Point_2 >
        >
        Point_generator;
        
        points.erase( points.begin(), points.end());
        n = max( n, 3);
        CGAL_random_convex_set_2( n,
                                  back_inserter( points),
                                  Point_generator( 1));
        
        polygon_changed = true;
        break;
      }
      else if ( input == compute_button)
        // recompute largest inscribing k-gon
        break;
      else if ( input == MOUSE_BUTTON( 1)) {
        // insert point:
        points.push_back( Point_2( x, y));
        polygon_changed = true;
        break;
      }
    } // for (;;)
  } // for (;;)

} // int main()


// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

