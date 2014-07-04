// Copyright (c) 2013 Technical University Braunschweig (Germany).
// All rights reserved.
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
// Author(s):  Ning Xu <longyin0904@gmail.com>
//

#ifndef CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H
#define CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H

#include <CGAL/Arrangement_2.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/Visibility_2/visibility_utils.h>
#include <vector>
#include <map>

//#define MYDEBUG
#ifdef MYDEBUG
  #include <iostream>
  using namespace std;
#endif

namespace CGAL {

template < class Arrangement_2_, class RegularizationTag = CGAL::Tag_true >
class Rotational_sweep_visibility_2
{
public:
  typedef Arrangement_2_                              Arrangement_2;
  typedef RegularizationTag                           Regularization_tag;
  typedef typename Arrangement_2::Geometry_traits_2   Geometry_traits_2;
  typedef typename Arrangement_2::Traits_2            Traits_2;
  typedef CGAL::Tag_true                              Supports_general_polygon_tag;
  typedef CGAL::Tag_true                              Supports_simple_polygon_tag;

  typedef typename Geometry_traits_2::FT              Number_type;
  typedef typename Geometry_traits_2::Point_2         Point_2;
  typedef typename Geometry_traits_2::Ray_2           Ray_2;
  typedef typename Geometry_traits_2::Segment_2       Segment_2;
  typedef typename Geometry_traits_2::Direction_2     Direction_2;

  typedef typename Arrangement_2::Face_const_handle   Face_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle
                                                      Halfedge_const_handle;
  typedef typename Arrangement_2::Hole_const_iterator Hole_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                      Circulator;


private:
  template < class Arr_2_ED >
  class Edge
  {
  public:
    typedef Arr_2_ED                                    Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2   Geometry_traits_2;

    typedef typename Geometry_traits_2::Point_2         Point_2;
    typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                        Circulator;
  private:
    Circulator circ;
    CGAL::Orientation orient;
    enum { NOT_APPLIED, INWARD, OUTWARD, IN_EDGE, AT_SOURCE, AT_TARGET } mode;
    int idx;
  public:
    Edge ( const Point_2& query_pt, const Circulator& he, bool on_edge, bool on_vertex, int _idx )
      : circ( he ), idx( _idx )
    {
      if ( on_vertex ) {
        orient = CGAL::COLLINEAR;
        if ( query_pt == he->source()->point() )
          mode = AT_SOURCE;
        else if ( query_pt == he->target()->point() )
          mode = AT_TARGET;
        else
          mode = IN_EDGE;
      } else if ( on_edge ) {
        orient = CGAL::COLLINEAR;
        mode = IN_EDGE;
      } else {
        orient = CGAL::orientation( query_pt,
                                    he->source()->point(),
                                    he->target()->point() );
        if ( orient == CGAL::COLLINEAR ) {
          if ( CGAL::collinear_are_ordered_along_line( query_pt,
                                                       he->source()->point(),
                                                       he->target()->point() ) )
            mode = OUTWARD;
          else
            mode = INWARD;
        } else {
            mode = NOT_APPLIED;
        }
      }
    }

    const Point_2& source () const
      { return circ->source()->point(); }
    const Point_2& target () const
      { return circ->target()->point(); }
    const Point_2& prev_source() const
      { return circ->prev()->source()->point(); }
    const Point_2& next_target() const
      { return circ->next()->target()->point(); }
    Circulator circulator () const
      { return circ; }
    CGAL::Orientation orientation () const
      { return orient; }
    bool inward () const
      { return ( mode == INWARD ); }
    bool outward () const
      { return ( mode == OUTWARD ); }
    bool in_edge () const
      { return ( mode == IN_EDGE ); }
    bool at_source () const
      { return ( mode == AT_SOURCE ); }
    bool at_target () const
      { return ( mode == AT_TARGET ); }
    int index () const
      { return idx; }

#ifdef MYDEBUG
    void trace ( ostream& os )
    {
      os << "source=[" << source() << "],target=[" << target() << "],";
      os << "orientation=[";
      switch( orient ) {
        case CGAL::LEFT_TURN:
          os << "left";
          break;
        case CGAL::RIGHT_TURN:
          os << "right";
          break;
        case CGAL::COLLINEAR:
          os << "collinear";
          break;
      }
      os << "],mode =[";
      switch( mode ) {
        case INWARD:
          os << "inward";
          break;
        case OUTWARD:
          os << "outward";
          break;
        case IN_EDGE:
          os << "in_edge";
          break;
        case AT_SOURCE:
          os << "at_source";
          break;
        case AT_TARGET:
          os << "at_target";
          break;
        case NOT_APPLIED:
          os << "not_applied";
          break;
      }
      os << "]" << endl;
    }
#endif
  };

  /*
    class Less_Source
    Compare the halfedges with their source point by their polar angle.
    The class rovides a comparator:
        Less_Source( const Edge_type& he1, const Edge_type& he2 )
    where he1 and he2 are two half edges.

    Precondition: he1 != he2

    Special cases:
      1) If two source points have the same angle, the halfedge goes outward
         or makes a right turn is less than the halfedge goes inward 
         or makes a left turn.
         If two halfedges both go outward or make a right turn, the one with
         closer source point is less.
         If two halfedges both go inward or make a left turn, the one with
         farther source point is less.
      2) If the source of an halfedge is the query point, consider the case
         as moving the source point slightly along the line through the
         halfedge, either toward the target or away the target, so that
         the query point is still in the face.
   */
  template < class E >
  class Less_Source
  {
  private:
    typedef E                                           Edge_type;
    typedef typename E::Geometry_traits_2               Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2         Point_2;
    typedef typename Geometry_traits_2::FT              Number_type;
  private:
    // angle_less_than pi
    // Verify whether the polar angle of p is less than pi
    bool angle_less_than_pi( const Point_2& p )
    {
      CGAL::Orientation orient = CGAL::orientation( query_pt, aux, p );
      return ( orient == CGAL::LEFT_TURN ||
               ( orient == CGAL::COLLINEAR && p.x() >= query_pt.x() ) );
    }
    // less_general
    // precondition: query_pt, p1 and p2 are not collinear
    // return true if p1 has smaller polar angle
    bool less_general ( const Point_2& p1, const Point_2& p2, CGAL::Orientation orient )
    {
      if ( angle_less_than_pi( p1 ) ) {
        if ( orient == CGAL::LEFT_TURN )
          return true;
        else
          return ( CGAL::orientation( query_pt, aux, p2 ) == CGAL::RIGHT_TURN );
      } else {
        if ( orient == CGAL::RIGHT_TURN )
          return false;
        else
          return ( CGAL::orientation( query_pt, aux, p2 ) == CGAL::RIGHT_TURN );
      }
    }
    // less_source_at_query_pt
    // Precondition: he1.source() == query_pt
    bool less_source_at_query_pt( const Edge_type& he1, const Edge_type& he2 )
    {
      CGAL::Orientation orient1 = CGAL::orientation( he1.prev_source(), he1.source(), he1.target() );
      CGAL::Orientation orient2 = CGAL::orientation( he1.source(), he1.target(), he2.source() );
      if ( orient1 == CGAL::LEFT_TURN ) {
        // The boundary makes a left turn at query_pt
        // Consider the case as moving he1.source() slightly
        // along the line through he1, away he1.target()
        if ( orient2 == CGAL::COLLINEAR ) {
          // he1.source(), he1.target(), he2.source() are collinear
          if ( CGAL::collinear_are_ordered_along_line( he2.source(), he1.source(), he1.target() ) ) {
            // he1.source() is between he1.target() and he2.source()
            // he1 will be considered going inward
            return false;
          } else {
            // he1.source(), he1.target(), he2.source() are ordered along ray
            return !angle_less_than_pi( he2.source() );
          }
        } else {
          Point_2 image( query_pt.x()+query_pt.x()-he1.target().x(), query_pt.y()+query_pt.y()-he1.target().y() );
          if ( orient2 == CGAL::LEFT_TURN )
            return less_general( image, he2.source(), CGAL::RIGHT_TURN );
          else
            return less_general( image, he2.source(), CGAL::LEFT_TURN );
        }
      } else {
        // The boundary makes a right turn at query_pt,
        // or does not make a turn at query_pt.
        // Consider the case as moving he1.source() slightly
        // along the line through he1, toward he1.target()
        if ( orient2 == CGAL::COLLINEAR ) {
          // he1.source(), he1.target(), he2.source() are collinear
          if ( CGAL::collinear_are_ordered_along_line( he2.source(), he1.source(), he1.target() ) ) {
            // he1.source() is between he1.target() and he2.source()
            return angle_less_than_pi( he1.target() );
          } else {
            // he1.source(), he1.target(), he2.source() are ordered along ray
            return true;
          }
        } else {
          return less_general( he1.target(), he2.source(), orient2 );
        }
      }
    }
    // less_on_ray
    // Precondition: he1.source() and he2.source() are collinear on a ray
    // shooting from the query_pt
    bool less_on_ray ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( CGAL::collinear_are_ordered_along_line( query_pt, he1.source(), he2.source() ) )
        // he1.source() is closer
        return ( he1.orientation() == CGAL::RIGHT_TURN ||
                 ( he1.at_source() || he1.outward() ) );
      else
        return !( he2.orientation() == CGAL::RIGHT_TURN ||
                  ( he2.at_source() || he2.outward() ) );
    }
  public:
    // Constructor
    Less_Source( const Point_2& q, const Point_2& a )
      : query_pt( q ), aux( a )
      {}
    // Comparator
    // Precondition: he1 and he2 are not the same halfedge
    bool operator () ( const Edge_type& he1, const Edge_type& he2 )
    {
      CGAL::Orientation orient = CGAL::orientation( query_pt, he1.source(), he2.source() );
      if ( orient != CGAL::COLLINEAR )
        return less_general( he1.source(), he2.source(), orient );
      // query_pt, he1.source(), he2.source() are collinear
      if ( he2.at_source() )
        return !less_source_at_query_pt( he2, he1 );
      if ( he1.at_source() )
        return less_source_at_query_pt( he1, he2 );
      if ( CGAL::collinear_are_ordered_along_line( he1.source(), query_pt, he2.source() ) )
        // querypt is between he1.source() and he2.source()
        return angle_less_than_pi( he1.source() );
      return less_on_ray ( he1, he2 );
    }
  private:
    Point_2 query_pt;
    Point_2 aux;
  };

  template < class E >
  class Less_Edge
  {
  private:
    typedef E                                           Edge_type;
    typedef typename E::Geometry_traits_2               Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2         Point_2;
    typedef typename Geometry_traits_2::FT              Number_type;
  private:
    // less_vertex, the halfedge he1 passes throught query_pt
    bool less_vertex ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he1.at_source() ) {
        if ( he2.at_target() ) {
          // he2 is the previous halfedge of he1
          return ( CGAL::orientation( query_pt, *p_aux, he1.target() ) != CGAL::COLLINEAR );
        } else {
          return true;
        }
      } else if ( he1.at_target() ) {
        if ( he2.at_source() ) {
          // he2 is the next halfege of he2
          return ( CGAL::orientation( query_pt, *p_aux, he2.target() ) == CGAL::COLLINEAR );
        } else {
          return true;
        }
      } else {
        // The query point lies in the interior of he1
        return true;
      }
    }
    // less_collinear, the halfedge he1 is through the ray (query_pt, *p_aux)
    bool less_collinear ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he1.inward() ) {
        if ( he2.inward() || he2.outward() )
          return ( CGAL::collinear_are_ordered_along_line( query_pt, he1.target(), he2.target() ) );
        else
          return ( CGAL::orientation( he2.source(), he2.target(), query_pt ) ==
                   CGAL::orientation( he2.source(), he2.target(), he1.target() ) );
      } else {
        // he1.outward()
        if ( he2.inward() || he2.outward() )
          return ( CGAL::collinear_are_ordered_along_line( query_pt, he1.source(), he2.source() ) );
        else
          return ( CGAL::orientation( he2.source(), he2.target(), query_pt ) ==
                   CGAL::orientation( he2.source(), he2.target(), he1.source() ) );
      }
    }
    // General setting, he2 is the next halfedge of he1
    bool less_consecutive ( const Edge_type& he1, const Edge_type& he2 )
    {
      CGAL::Orientation orient1 = CGAL::orientation( query_pt, *p_aux, he1.target() );
      if ( orient1 == CGAL::LEFT_TURN ) {
        return ( CGAL::orientation( he1.source(), he1.target(), he2.target() ) == CGAL::RIGHT_TURN );
      } else if ( orient1 == CGAL::RIGHT_TURN ) {
        return ( CGAL::orientation( he1.source(), he1.target(), he2.target() ) == CGAL::LEFT_TURN );
      } else {
        // orient1 == CGAL::COLLINEAR
        if ( he1.orientation() == CGAL::LEFT_TURN ) {
          if ( he2.orientation() == CGAL::LEFT_TURN )
            return false;
          else
            return ( CGAL::orientation( he1.source(), he1.target(), he2.target() ) == CGAL::RIGHT_TURN );
        } else {
          // he1.orientation() == CGAL::RIGHT_TURN
          if ( he2.orientation() == CGAL::RIGHT_TURN )
            return true;
          else
            return ( CGAL::orientation( he1.source(), he1.target(), he2.target() ) == CGAL::LEFT_TURN );
        }
      }
    }
    // less_general, the general case
    bool less_general ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he1.source() == he2.target() ) 
        return !less_consecutive( he2, he1 );
      if ( he1.target() == he2.source() )
        return less_consecutive( he1, he2 );
      CGAL::Orientation orient1 = CGAL::orientation( he1.source(), he1.target(), he2.source() );
      CGAL::Orientation orient2 = CGAL::orientation( he1.source(), he1.target(), he2.target() );
      if ( orient1 == orient2 ) {
        return ( CGAL::orientation( he1.source(), he1.target(), query_pt ) != orient1 );
      } else {
        return ( CGAL::orientation( he2.source(), he2.target(), query_pt ) ==
                 CGAL::orientation( he2.source(), he2.target(), he1.source() ) );
      }
    }
  public:
    Less_Edge ( const Point_2& q, const Point_2 * a )
      : query_pt( q ), p_aux( a )
      {}
    bool operator () ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he1.circulator() == he2.circulator() )
        return false;
      if ( he1.at_source() || he1.at_target() || he1.in_edge() )
        return less_vertex( he1, he2 );
      if ( he2.at_source() || he2.at_target() || he2.in_edge() )
        return !less_vertex( he2, he1 );
      if ( he1.inward() || he1.outward() )
        return less_collinear( he1, he2 );
      if ( he2.inward() || he2.outward() )
        return !less_collinear( he2, he1 );
      return less_general( he1, he2 );
    }
  private:
    Point_2 query_pt;
    const Point_2 * p_aux;
  };

private:
  typedef Edge<Arrangement_2>                         Edge_type;
  typedef std::vector<Edge_type>                      Edge_vector;
  typedef typename Edge_vector::const_iterator        Edge_iterator;
  typedef bool                                        Active_edge_property;
  typedef Less_Edge<Edge_type>                        Active_edge_sorter;
  typedef std::map<Edge_type,Active_edge_property,Active_edge_sorter>
                                                      Active_edge_container;
  typedef typename Active_edge_container::const_iterator
                                                      Active_const_iterator;

private:
  // do_intersect_ray
  // Verify whether the halfedge he will intersect the ray obtained by
  // by slightly rotating the ray (query_pt, aux ) clockwisely
  bool do_intersect_ray ( const Point_2& query_pt, const Point_2& aux, const Edge_type& he )
  {
    if ( he.orientation() != CGAL::COLLINEAR ) {
      Ray_2 ray( query_pt, aux );
      Segment_2 seg( he.source(), he.target() );
      if ( !CGAL::do_intersect( ray, seg ) )
        return false;
      if ( he.orientation() == CGAL::LEFT_TURN )
        return ( CGAL::orientation( query_pt, aux, he.source() ) == CGAL::RIGHT_TURN );
      else
        return ( CGAL::orientation( query_pt, aux, he.target() ) == CGAL::RIGHT_TURN );
    } else if ( he.inward() || he.outward() ) {
      return false;
    } else if ( he.at_source() ) {
      CGAL::Orientation orient1 = CGAL::orientation( he.prev_source(), he.source(), he.target() );
      if ( orient1 == CGAL::LEFT_TURN ) {
        // The boundary makes a left turn at the query_pt
        // Consider the case as moving he.source() slightly along the line
        // through he, away he.target()
        CGAL::Orientation orient2 = CGAL::orientation( he.source(), he.target(), aux );
        if ( orient2 == CGAL::LEFT_TURN ) {
          return false;
        } else if ( orient2 == CGAL::RIGHT_TURN ) {
          return true;
        } else {
          // he.source(), he.target(), aux ) are collinear 
          return !CGAL::collinear_are_ordered_along_line( aux, he.source(), he.target() );
        }
      } else {
        // The boundary makes a right turn or does not turn at the query_pt
        // Consider the case as moving he.source() slightly along the line
        // through he, toward he.target()
        return false;
      }
    } else if ( he.at_target() ) {
      CGAL::Orientation orient1 = CGAL::orientation( he.source(), he.target(), he.next_target() );
      if ( orient1 == CGAL::LEFT_TURN ) {
        // The boundary makes a left turn at the query_pt
        CGAL::Orientation orient2 = CGAL::orientation( he.target(), he.next_target(), aux );
        if ( orient2 == CGAL::LEFT_TURN ) {
          return ( CGAL::orientation( he.source(), he.target(), aux ) == CGAL::RIGHT_TURN );
        } else if ( orient2 == CGAL::RIGHT_TURN ) {
          return false;
        } else {
          return CGAL::collinear_are_ordered_along_line( aux, he.target(), he.next_target() );
        }
      } else if ( orient1 == CGAL::RIGHT_TURN ) {
        // The boundary makes a right turn at the query_pt
        CGAL::Orientation orient2 = CGAL::orientation( he.target(), he.next_target(), aux );
        if ( orient2 == CGAL::LEFT_TURN ) {
          return false;
        } else if ( orient2 == CGAL::RIGHT_TURN ) {
          return ( CGAL::orientation( he.source(), he.target(), aux ) == CGAL::RIGHT_TURN );
        } else {
          return !CGAL::collinear_are_ordered_along_line( aux, he.target(), he.next_target() );
        }
      }
    } else {
      // he.in_edge() == true
      CGAL::Orientation orient1 = CGAL::orientation( query_pt, he.target(), aux );
      if ( orient1 == CGAL::LEFT_TURN ) {
        return false;
      } else if ( orient1 == CGAL::RIGHT_TURN ) {
        return true;
      } else {
        return !CGAL::collinear_are_ordered_along_line( aux, query_pt, he.target() );
      }
    }
  }
  Point_2 calculate_intersection( const Point_2& query_pt, const Point_2& aux, const Circulator& he )
  {
     Ray_2 ray ( query_pt, aux );
     Segment_2 seg ( he->source()->point(), he->target()->point() );
     CGAL::Object res = CGAL::intersection( ray, seg );
     const Point_2 * ipoint = CGAL::object_cast<Point_2>(&res);
     if ( ipoint ) {
       return *ipoint;
     } else {
       assert( seg.has_on( query_pt ) );
       return query_pt;
     }
  }
  void compute_visibility_partition ( const Point_2& query_pt,
                                      Edge_iterator first,
                                      Edge_iterator last,
                                      std::vector<Point_2>& out_points )
  {
    Point_2 aux;
    Active_edge_sorter closer( query_pt, &aux );
    Active_edge_container active_edges( closer );
    std::vector<bool> active;
    Active_const_iterator ait;

    // Initialize the edges intersecting the ray
    aux = first->source();
    if ( query_on_vertex && query_pt == aux ) {
      CGAL::Orientation orient1 = CGAL::orientation( first->prev_source(), first->source(), first->target() );
      if ( orient1 == CGAL::LEFT_TURN ) {
        aux = Point_2( query_pt.x()+query_pt.x()-first->target().x(), query_pt.y()+query_pt.y()-first->target().y() );
      } else {
        aux = first->target();
      }
    }

    active.assign( unsorted_edges.size(), false );
    for ( Edge_iterator eit = edges.begin(); eit != edges.end(); eit++ ) {
      if ( do_intersect_ray( query_pt, aux, *eit ) ) { 
        active_edges.insert( std::make_pair( *eit, true ) );
        active[eit->index()] = true;
      }
    }

#ifdef MYDEBUG
cout << "query_pt = [" << query_pt << "]" <<endl;
cout << "Unsorted edges" << endl;
cout << "================================" << endl;
for ( int i = 0; i < unsorted_edges.size(); i++ ) {
  cout << "Edge: " << unsorted_edges[i].source() << " --> " << unsorted_edges[i].target() << "  idx = " << unsorted_edges[i].index() << endl;
}
cout << endl;
cout << "Sorted edges" << endl;
cout << "================================" << endl;
for ( int i = 0; i < edges.size(); i++ ) {
  cout << "Edge: " << edges[i].source() << " --> " << edges[i].target() << "  idx = " << edges[i].index() << endl;
}
cout << endl;
cout << "After Initialization" << endl;
cout << "================================" << endl;
for ( ait = active_edges.begin(); ait != active_edges.end(); ait++ ) {
  cout << "Edge: " << ait->first.source() << " ---> " << ait->first.target() << endl;
}
cout << endl;
#endif 

    // Rotational sweep the ray, until reach the end of the cone
    Point_2 first_intersection, last_intersection;
    first_intersection = last_intersection = calculate_intersection( query_pt, aux, active_edges.begin()->first.circulator() );
    bool first_intersect_at_vertex = ( first_intersection == active_edges.begin()->first.source() ||  first_intersection == active_edges.begin()->first.target() );
    for ( Edge_iterator eit = first; eit != last; eit++ ) {
      aux = eit->source();
      int idx = eit->index();
      int prev = ( idx-1 );
      if ( prev < 0 )
        prev = unsorted_edges.size() - 1;
      Circulator top = active_edges.begin()->first.circulator();
      assert( unsorted_edges[idx].circulator() == eit->circulator() );

#ifdef MYDEBUG
cout << "idx=" << idx << ",prev=" << prev << ",size=" << unsorted_edges.size() << endl;
cout << "Current edge: " << eit->source() << " --> " << eit->target() << "    Previous edge: " << unsorted_edges[prev].source() << " --> " << unsorted_edges[prev].target() << endl;
cout << "aux = " << aux << "   top: " << top->source()->point() << " --> " << top->target()->point() << endl;
#endif

      if ( active[idx] && active[prev] ) {
        // Both edges incident to the current vertex are active
#ifdef MYDEBUG
cout << "Both Active!" << endl;
#endif

        active_edges.erase( *eit );
        active[idx] = false;
        active_edges.erase( unsorted_edges[prev] );
        active[prev] = false;
        if ( top != active_edges.begin()->first.circulator() ) {
          Point_2 u;
          u = calculate_intersection( query_pt, aux, active_edges.begin()->first.circulator() );
          if ( last_intersection != eit->source() )
            out_points.push_back( eit->source() );
          out_points.push_back( u );
          last_intersection = u;
#ifdef MYDEBUG
cout << "New Top! Intersection = " << u << endl;
#endif
        }
      } else if ( active[idx] ) {
#ifdef MYDEBUG
cout << "Current Active!" << endl;
#endif
        // Only one edge whose source is the current vertex is active.
        active_edges.erase( *eit );
        active[idx] = false;
        active_edges.insert( std::make_pair( unsorted_edges[prev], true ) );
        active[prev] = true;
        if ( top != active_edges.begin()->first.circulator() ) {
          if ( last_intersection != eit->source() ) {
            out_points.push_back( eit->source() );
            last_intersection = eit->source();
          }
#ifdef MYDEBUG
cout << "New Top! Intersection = " << eit->source() << endl;
#endif
        }
      } else if ( active[prev] ) {
#ifdef MYDEBUG
cout << "Previous Active!" << endl;
#endif
        // Only one edge whose target is the current vertex is active.
        active_edges.erase( unsorted_edges[prev] );
        active[prev] = false;
        active_edges.insert( std::make_pair( *eit, true ) );
        active[idx] = true;
        if ( top != active_edges.begin()->first.circulator() ) {
          if ( last_intersection != eit->source() ) {
            out_points.push_back( eit->source() );
            last_intersection = eit->source();
          }
#ifdef MYDEBUG
cout << "New Top! Intersection = " << eit->source() << endl;
#endif
        }
      } else {
        // Both edges incident to the current vertex are not active
#ifdef MYDEBUG
cout << "Both Inctive!" << endl;
#endif
        active_edges.insert( std::make_pair( *eit, true ) );
        active[idx] = true;
        active_edges.insert( std::make_pair( unsorted_edges[prev], true ) );
        active[prev] = true;
        if ( top != active_edges.begin()->first.circulator() ) {
          Point_2 u;
          if ( query_on_vertex && query_pt == aux )
            u = aux;
          else
            u = calculate_intersection( query_pt, aux, top );
          if ( last_intersection != u )
            out_points.push_back( u );
          out_points.push_back( eit->source() );
          last_intersection = eit->source();
#ifdef MYDEBUG
cout << "New Top! Intersection = " << u << endl;
#endif
        }
      }
#ifdef MYDEBUG
cout << "*** After Iteration ***" << endl;
for ( ait = active_edges.begin(); ait != active_edges.end(); ait++ ) {
  cout << "Edge: " << ait->first.source() << " ---> " << ait->first.target() << endl;
}
cout << endl;
#endif
    }
    if ( first_intersection != last_intersection ) {
      if ( first_intersect_at_vertex ||
        ( CGAL::orientation( out_points.front(), out_points.back(), first_intersection ) != CGAL::COLLINEAR ) )
      out_points.push_back( first_intersection );
    }
#ifdef MYDEBUG
cout << "Visibility segments:" << endl;
for ( int i = 0; i < out_points.size(); i++ ) {
  cout << out_points[i] << endl;
}
cout << endl;
#endif
  }


  template < class VARR >
  typename VARR::Face_handle
  compute_visibility_impl ( const Point_2& query_pt, VARR& arr_out )
  {
    Point_2 aux( query_pt.x()+Number_type(1), query_pt.y() );
    Less_Source<Edge_type> comp ( query_pt, aux );
    // Sort halfedges with their source point by polar angle
#ifdef MYDEBUG
for ( int i = 0; i < edges.size(); i++ ) {
  for ( int j = 0; j < edges.size(); j++ ) {
    if ( i== j )
      continue;
    bool res = comp( edges[i], edges[j] );
    cout << "Edges1: " << edges[i].source() << " --> " << edges[i].target() << "    ";
    cout << "Edges2: " << edges[j].source() << " --> " << edges[j].target() << "    ";
    if ( res )
      cout << "Result:  1<2" <<endl;
    else
      cout << "Result:  2<1" <<endl;
  }
}
#endif
    std::sort( edges.begin(), edges.end(), comp );

    std::vector<Point_2> vec;
    compute_visibility_partition( query_pt, edges.begin(), edges.end(), vec );

    // Construct arrangement
    CGAL::Visibility_2::report_while_handling_needles
                        < Rotational_sweep_visibility_2 >
                        ( geom_traits, query_pt, vec, arr_out );

    conditional_regularize( arr_out, Regularization_tag() );

    edges.clear();
    unsorted_edges.clear();
    return arr_out.faces_begin();
  }

  /*! Regularize output if flag is set to true*/
  template <typename VARR> 
  void conditional_regularize(VARR& out_arr, CGAL::Tag_true) {
    regularize_output(out_arr);
  }
  /*! No need to regularize output if flag is set to false*/
  template <typename VARR> 
  void conditional_regularize(VARR& out_arr, CGAL::Tag_false) {
    //do nothing
  }

  /*! Regularizes the output - removes edges that have the same face on both
      sides */
  template <typename VARR> 
  void regularize_output(VARR& out_arr) {
    typename VARR::Edge_iterator e_itr;
    for (e_itr = out_arr.edges_begin() ; 
         e_itr != out_arr.edges_end() ; e_itr++) {

      typename VARR::Halfedge_handle he = e_itr;
      typename VARR::Halfedge_handle he_twin = he->twin();
      if (he->face() == he_twin->face()) {
        out_arr.remove_edge(he);
      }
    }
  }

public:
  // Constructor
  Rotational_sweep_visibility_2 ()
    : p_arr( NULL ), geom_traits( NULL )
    {}
  Rotational_sweep_visibility_2 ( const Arrangement_2& arr )
    : p_arr( &arr )
    { geom_traits = p_arr->geometry_traits(); }

  const std::string name ()
    { return std::string( "R_visibility_2" ); }
  bool is_attached () const
    { return (p_arr != NULL); }
  void attach ( const Arrangement_2& arr )
    { p_arr = &arr; geom_traits = p_arr->geometry_traits(); }
  void detach ()
    { p_arr = NULL; geom_traits = NULL; }
  const Arrangement_2& arr () const
    { return *p_arr; }

  template < typename VARR >
  typename VARR::Face_handle
  compute_visibility ( const Point_2& query_pt,
                       const Halfedge_const_handle he,
                       VARR& arr_out )
  {
    Face_const_handle f = he->face();
    assert( !f->is_unbounded() );

    Halfedge_const_handle he2 = he;
    query_on_edge = true;
    query_on_vertex = false;
    if ( query_pt == he2->target()->point() )
      he2 = he2->next();
    if ( query_pt == he2->source()->point() )
      query_on_vertex = true;

    arr_out.clear();
    edges.clear();
    Circulator circ, curr;
    Circulator qedge( he2 );
    Circulator qprev( he2->prev() );

    int idx = 0;

    if ( f->has_outer_ccb() ) {
      curr = circ = f->outer_ccb();
      do {
        bool on_edge = false, on_vertex = false;
        if ( curr == qedge ) {
          on_edge = true;
          on_vertex = query_on_vertex;
        } else if ( curr == qprev ) {
          on_edge = on_vertex = query_on_vertex;
        }
        edges.push_back( Edge_type( query_pt, curr, on_edge, on_vertex, idx ) );
        idx++;
      } while ( ++curr != circ );
    }
    for ( Hole_const_iterator hi = f->holes_begin();
          hi != f->holes_end(); hi++ ) {
      curr = circ = *hi;
      do {
        bool on_edge = false, on_vertex = false;
        if ( curr == qedge ) {
          on_edge = true;
          on_vertex = query_on_vertex;
        } else if ( curr == qprev ) {
          on_edge = on_vertex = query_on_vertex;
        }
        edges.push_back( Edge_type( query_pt, curr, on_edge, on_vertex, idx ) );
        idx++;
      } while ( ++curr != circ );
    }

    unsorted_edges.assign( edges.begin(), edges.end() );

    return compute_visibility_impl( query_pt, arr_out );
  }

  template < typename VARR >
  typename VARR::Face_handle
  compute_visibility ( const Point_2& query_pt,
                       const Face_const_handle f,
                       VARR& arr_out )
  {
    assert( !f->is_unbounded() );

    query_on_edge = query_on_vertex = false;
    arr_out.clear();
    edges.clear();

    Circulator circ, curr;

    int idx = 0;

    if ( f->has_outer_ccb() ) {
      curr = circ = f->outer_ccb();
      do {
        edges.push_back( Edge_type( query_pt, curr, false, false, idx ) );
        idx++;
      } while ( ++curr != circ );
    }
    for ( Hole_const_iterator hi = f->holes_begin();
          hi != f->holes_end(); hi++ ) {
      curr = circ = *hi;
      do {
        edges.push_back( Edge_type( query_pt, curr, false, false, idx ) );
        idx++;
      } while ( ++curr != circ );
    }

    unsorted_edges.assign( edges.begin(), edges.end() );

    return compute_visibility_impl( query_pt, arr_out );
  }

private:
  const Arrangement_2 * p_arr;
  const Geometry_traits_2 * geom_traits;

  bool query_on_edge;
  bool query_on_vertex;
  Edge_vector edges;
  Edge_vector unsorted_edges;

}; // End of class Rotational_sweep_visibility_2

} // End namespace CGAL

#endif
