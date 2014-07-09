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
#include <set>

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
    typedef typename Geometry_traits_2::FT              Number_type;
    typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                        Circulator;
  private:
    Circulator circ;
    CGAL::Orientation orient;
    CGAL::Orientation boundary_orient;
    enum { NOT_APPLIED, INWARD, OUTWARD, IN_EDGE, AT_SOURCE, AT_TARGET } mode;
    int idx;
    int prev_idx;
    int quad;
  private:
    int compute_quad( const Point_2& query_pt, const Point_2& p )
    {
      CGAL::Comparison_result cx = CGAL::compare_x( query_pt, p );
      CGAL::Comparison_result cy = CGAL::compare_y( query_pt, p );
      if ( cy == CGAL::SMALLER ) {
        if ( cx == CGAL::SMALLER )
          return 0;
        else
          return 1;
      } else if ( cy == CGAL::LARGER ) {
        if ( cx == CGAL::LARGER )
          return 2;
        else
          return 3;
      } else {
        if ( cx != CGAL::LARGER )
          return 0;
        else
          return 2;
      }
    }
  public:
    Edge ( const Point_2& query_pt, const Circulator& he, bool on_edge, bool on_vertex )
      : circ( he )
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
      quad = compute_quad( query_pt, he->source()->point()  );
      boundary_orient = CGAL::orientation( prev_source(), source(), target() );
    }
    void set_index ( int i, int p )
      { idx = i; prev_idx = p; }

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
    CGAL::Orientation boundary_orientation () const
      { return boundary_orient; }
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
    int prev_index () const
      { return prev_idx; }
    int quadrant () const
      { return quad; }
    bool angle_less_than_pi () const
      { return ( quadrant() < 2 ); }

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
      os << "\t\tquadrant=[" << quadrant() << "], idx=[" << index() << "], prev_idx=[" << prev_index() << "]" << endl;
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
    // less_source_at_query_pt
    // Precondition: he1.source() == query_pt
    bool less_source_at_query_pt( const Edge_type& he1, const Edge_type& he2 )
    {
      CGAL::Orientation orient1 = he1.boundary_orientation();
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
            return !he2.angle_less_than_pi();
          }
        } else if ( orient2 == CGAL::LEFT_TURN ) {
          // he1.source(), he1.target(), he2.source() make a left turn
          if ( he2.angle_less_than_pi() ) {
            return false;
          } else {
            return ( he1.target().y() < query_pt.y() ||
                     ( he1.target().y() == query_pt.y() &&
                       he1.target().x() < query_pt.x() ) );
          }
        } else {
          // he1.source(), he1.target(), he2.source() make a right turn
          if ( he2.angle_less_than_pi() ) {
            return ( he1.target().y() < query_pt.y() ||
                     ( he1.target().y() == query_pt.y() &&
                       he1.target().x() < query_pt.x() ) );
          } else {
            return true;
          }
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
            return !he2.angle_less_than_pi();
          } else {
            // he1.source(), he1.target(), he2.source() are ordered along ray
            return true;
          }
        } else if ( orient2 == CGAL::LEFT_TURN ) {
          // he1.source(), he1.target(), he2.source() make a left turn
          if ( he2.angle_less_than_pi() ) {
            return ( he1.target().y() > query_pt.y() ||
                     ( he1.target().y() == query_pt.y() &&
                       he1.target().x() > query_pt.x() ) );
          } else {
            return true;
          }
        } else {
          // he1.source(), he1.target(), he2.source() make a right turn
          if ( he2.angle_less_than_pi() ) {
            return false;
          } else {
            return ( he1.target().y() > query_pt.y() ||
                     ( he1.target().y() == query_pt.y() &&
                       he1.target().x() > query_pt.x() ) );
          }
        }
      }
    }
    // less
    bool less ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he2.at_source() )
        return !less_source_at_query_pt( he2, he1 );
      if ( he1.at_source() )
        return less_source_at_query_pt( he1, he2 );

      if ( he1.quadrant() != he2.quadrant() )
        return ( he1.quadrant() < he2.quadrant() );
      // In the same quadrant
      CGAL::Orientation orient = CGAL::orientation( query_pt, he1.source(), he2.source() );
      if ( orient != CGAL::COLLINEAR ) {
        // General case, in the same quadrant
        return ( orient == CGAL::LEFT_TURN );
      } else {
        // query_pt, he1.source(), he2.source() are collinear on a ray
        if ( CGAL::collinear_are_ordered_along_line( query_pt, he1.source(), he2.source() ) ) {
          // he1.source() is closer
          return ( he1.orientation() == CGAL::RIGHT_TURN || he1.outward() );
        } else {
          // he2.source() is closer
          return !( he2.orientation() == CGAL::RIGHT_TURN || he2.outward() );
        }
      }
    }
  public:
    // Constructor
    Less_Source( const Point_2& q )
      : query_pt( q )
      {}
    // Comparator
    // Precondition: he1 and he2 are not the same halfedge
    bool operator () ( const Edge_type& he1, const Edge_type& he2 )
      { return less( he1, he2 ); }
  private:
    Point_2 query_pt;
    Point_2 aux;
  };

  /*
    class Less_Edge
    Compare the halfedges intersecting a ray, find whose intersection point
    is closer to the query point.
    The class rovides a comparator:
        Less_Edge( const Edge_type& he1, const Edge_type& he2 )
    where he1 and he2 are two half edges.

    Precondition: he1 != he2
  */
  template < class E >
  class Less_Edge
  {
  private:
    typedef E                                           Edge_type;
    typedef typename E::Geometry_traits_2               Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2         Point_2;
  private:
    Point_2 query_pt;
  private:
    // less_vertex
    // Precondition:  (1) he1 contains query_pt
    //                (2) he2 does not contain query_pt
    bool less_vertex ( const Edge_type& he1, const Edge_type& he2 )
    {
      assert( !(he2.at_source() || he2.at_target() || he2.in_edge() ) );
      return true;
    }
    // less_consecutive
    // Precondition:  (1) Both he1 and he2 does not contain query_pt
    //                (2) he1 is the previous halfedge of he2
    bool less_consecutive ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he1.outward() )
        return true;
      else if ( he1.inward() )
        return false;
      else if ( he1.orientation() == CGAL::LEFT_TURN ) {
        return ( he2.boundary_orientation() == CGAL::RIGHT_TURN );
      } else {
        // he1.orientation() == CGAL::RIGHT_TURN
        return ( he2.boundary_orientation() != CGAL::RIGHT_TURN );
      }
    }
    // less_collinear
    // Precondition:  (1) Both he1 and he2 does not contain query_pt
    //                (2) he1 and he2 are not incident
    //                (3) query_pt, he1.source(), he1.target() are collinear
    bool less_collinear ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he2.inward() || he2.outward() ) {
        return CGAL::collinear_are_ordered_along_line( query_pt,
                                                       he1.source(),
                                                       he2.source() );
      } else {
        return ( CGAL::orientation( he2.source(), he2.target(), he1.source() )
                                   == he2.orientation() );
      }
    }
    // less_general
    // Precondition:  (1) Both he1 and he2 does not contain query_pt
    //                (2) he1 and he2 are not incident
    //                (3) he1.orientation() == CGAL::LEFT_TURN || RIGHT_TURN
    //                (4) he2.orientation() == CGAL::LEFT_TURN || RIGHT_TURN
    bool less_general ( const Edge_type& he1, const Edge_type& he2 )
    {
      CGAL::Orientation orient1 = CGAL::orientation( he1.source(),
                                                     he1.target(), 
                                                     he2.source() );
      CGAL::Orientation orient2 = CGAL::orientation( he1.source(),
                                                     he1.target(), 
                                                     he2.target() );
      if ( orient1 == orient2 ) {
        // he2.source() and he2.target() lies on the same side of he1
        return ( CGAL::orientation( he1.source(), he1.target(), query_pt )
                                  != orient1 );
      } else {
        // he2.source() and he2.target() lies on the different side of he1
        return ( CGAL::orientation( he2.source(), he2.target(), he1.source() )
                                 == he2.orientation() );
      }
    }
    // less
    bool less ( const Edge_type& he1, const Edge_type& he2 )
    {
      if ( he1.circulator() == he2.circulator() )
        return false;
      if ( he1.at_source() || he1.at_target() || he1.in_edge() )
        return less_vertex( he1, he2 );
      if ( he2.at_source() || he2.at_target() || he2.in_edge() )
        return !less_vertex( he2, he1 );
      if ( he1.index() == he2.prev_index() )
        return less_consecutive( he1, he2 );
      if ( he1.prev_index() == he2.index() )
        return !less_consecutive( he2, he1 );
      if ( he1.inward() || he1.outward() )
        return less_collinear( he1, he2 );
      if ( he2.inward() || he2.outward() )
        return !less_collinear( he2, he1 );
      return less_general( he1, he2 );
    }
  public:
    Less_Edge ( const Point_2& q )
      : query_pt( q )
      {}
    bool operator () ( const Edge_type& he1, const Edge_type& he2 )
      { return less( he1, he2 ); }
  };

  template < class E, class Comp >
  class Active_Edge
  {
  private:
    typedef E                                           Edge_type;
    typedef Comp                                        Sorter;
    typedef typename E::Geometry_traits_2               Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2         Point_2;
    typedef typename std::set<Edge_type,Sorter>         Active_edge_container;
    typedef typename std::vector<Edge_type>             Edge_vector;
  private:
    Active_edge_container active_edges;
    const Edge_vector * p_edges;
    std::vector<bool> active;
  public:
    Active_Edge ( const Edge_vector * pe, const Sorter& s )
      : p_edges( pe ), active_edges( s )
    { active.assign( pe->size(), false ); }

    void insert ( const Edge_type& he )
    {
      active_edges.insert( he );
      active[he.index()] = true;
    }
    void erase ( const Edge_type& he )
    {
      active_edges.erase( he );
      active[he.index()] = false;
    }
    void replace ( const Edge_type& he1, const Edge_type& he2 )
    {
      typename Active_edge_container::const_iterator ait;
      ait = active_edges.find( he1 );
      assert( ait != active_edges.end() );
      const Edge_type& const_he = *ait;
      Edge_type& tmp_he = const_cast<Edge_type&>(const_he);
      tmp_he = he2;
      active[he1.index()] = false;
      active[he2.index()] = true;
    }
    bool is_active( int idx ) const
      { return active[idx]; }
    const Edge_type& closest () const
      { return *(active_edges.begin()); }
#ifdef MYDEBUG
    void trace ( ostream& os )
    {
      typename Active_edge_container::const_iterator ait;
      for( ait = active_edges.begin(); ait != active_edges.end(); ait++ ) {
        cout << "Edge: " << ait->source() << " --> " << ait->target() << endl;
      }
        cout << "Active:[ ";
      for ( int i = 0; i < active.size(); i++ ) {
        if ( active[i] )
          cout << i << " ";
      }
      cout << "]" << endl;
    }
#endif
  };

private:
  typedef Edge<Arrangement_2>                         Edge_type;
  typedef std::vector<Edge_type>                      Edge_vector;
  typedef typename Edge_vector::const_iterator        Edge_iterator;
  typedef Less_Edge<Edge_type>                        Active_edge_sorter;

private:
  // do_intersect_ray
  // Verify whether the halfedge he will intersect the ray obtained by
  // by slightly rotating the ray (query_pt, aux ) clockwisely
  bool do_intersect_ray ( const Point_2& query_pt, const Point_2& aux, const Edge_type& he )
  {
    if ( he.orientation() != CGAL::COLLINEAR ) {
      CGAL::Orientation orient1 = CGAL::orientation( query_pt, aux, he.source() );
      CGAL::Orientation orient2 = CGAL::orientation( query_pt, aux, he.target() );
      if ( orient1 == orient2 )
        return false;
      if ( orient1 == CGAL::COLLINEAR ) {
        if ( CGAL::collinear_are_ordered_along_line( aux, query_pt, he.source() ) )
          return false;
        // Ray intersects he at he.source()
        return ( he.orientation() == CGAL::RIGHT_TURN );
      } else if ( orient2 == CGAL::COLLINEAR ) {
        if ( CGAL::collinear_are_ordered_along_line( aux, query_pt, he.target() ) )
          return false;
        // Ray intersects he at he.target()
        return ( he.orientation() == CGAL::LEFT_TURN );
      } else {
        return ( CGAL::orientation( query_pt, aux, he.target() ) == he.orientation() );
      }
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
    Active_edge_sorter closer( query_pt );
    Active_Edge< Edge_type, Active_edge_sorter > active( &unsorted_edges, closer );

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

    for ( Edge_iterator eit = edges.begin(); eit != edges.end(); eit++ ) {
      if ( do_intersect_ray( query_pt, aux, *eit ) ) { 
        active.insert( *eit );
      }
    }

#ifdef MYDEBUG
cout << "query_pt = [" << query_pt << "]" <<endl;
cout << "Unsorted edges" << endl;
cout << "================================" << endl;
for ( int i = 0; i < unsorted_edges.size(); i++ ) {
  unsorted_edges[i].trace( cout );
}
cout << endl;
cout << "Sorted edges" << endl;
cout << "================================" << endl;
for ( int i = 0; i < edges.size(); i++ ) {
  edges[i].trace( cout );
}
cout << endl;
cout << "After Initialization" << endl;
cout << "================================" << endl;
active.trace( cout );
cout << endl;
#endif 

    // Rotational sweep the ray, until reach the end of the cone
    Point_2 first_pt, last_pt;
    first_pt = last_pt = calculate_intersection( query_pt, aux, active.closest().circulator() );
    bool first_at_vertex = ( first_pt == active.closest().source() ||  first_pt == active.closest().target() );
    for ( Edge_iterator eit = first; eit != last; eit++ ) {
      aux = eit->source();
      int idx = eit->index();
      int prev = eit->prev_index();
      Circulator top = active.closest().circulator();
      assert( unsorted_edges[idx].circulator() == eit->circulator() );

#ifdef MYDEBUG
cout << "idx = " << idx << " , prev = " << prev << endl;
cout << "Current edge: " << eit->source() << " --> " << eit->target() << "    Previous edge: " << unsorted_edges[prev].source() << " --> " << unsorted_edges[prev].target() << endl;
cout << "top: " << top->source()->point() << " --> " << top->target()->point() << endl;
#endif

      if ( active.is_active( idx ) && active.is_active( prev ) ) {
        // Both edges incident to the current vertex are active
#ifdef MYDEBUG
cout << "Both Active!" << endl;
#endif

        active.erase( *eit );
        active.erase( unsorted_edges[prev] );
        if ( top != active.closest().circulator() ) {
          Point_2 u;
          u = calculate_intersection( query_pt, aux, active.closest().circulator() );
          if ( last_pt != eit->source() )
            out_points.push_back( eit->source() );
          out_points.push_back( u );
          last_pt = u;
#ifdef MYDEBUG
cout << "New Top! Intersection = " << u << endl;
#endif
        }
      } else if ( active.is_active( idx ) ) {
#ifdef MYDEBUG
cout << "Current Active!" << endl;
#endif
        // Only one edge whose source is the current vertex is active.
        active.replace( *eit, unsorted_edges[prev] );
        if ( top != active.closest().circulator() ) {
          if ( last_pt != eit->source() ) {
            out_points.push_back( eit->source() );
            last_pt = eit->source();
          }
#ifdef MYDEBUG
cout << "New Top! Intersection = " << eit->source() << endl;
#endif
        }
      } else if ( active.is_active( prev ) ) {
#ifdef MYDEBUG
cout << "Previous Active!" << endl;
#endif
        // Only one edge whose target is the current vertex is active.
        active.replace( unsorted_edges[prev], *eit );
        if ( top != active.closest().circulator() ) {
          if ( last_pt != eit->source() ) {
            out_points.push_back( eit->source() );
            last_pt = eit->source();
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
        active.insert( *eit );
        active.insert( unsorted_edges[prev] );
        if ( top != active.closest().circulator() ) {
          Point_2 u;
          if ( query_on_vertex && query_pt == aux )
            u = aux;
          else
            u = calculate_intersection( query_pt, aux, top );
          if ( last_pt != u )
            out_points.push_back( u );
          out_points.push_back( eit->source() );
          last_pt = eit->source();
#ifdef MYDEBUG
cout << "New Top! Intersection = " << u << endl;
#endif
        }
      }
#ifdef MYDEBUG
cout << "*** After Iteration ***" << endl;
active.trace( cout );
cout << endl;
#endif
    }

    if ( first_pt != last_pt ) {
      if ( first_at_vertex ||
        ( CGAL::orientation( out_points.front(), out_points.back(), first_pt ) != CGAL::COLLINEAR ) )
      out_points.push_back( first_pt );
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
    Less_Source<Edge_type> comp ( query_pt );
    // Sort halfedges with their source point by polar angle
    std::sort( edges.begin(), edges.end(), comp );
#ifdef MYDEBUG
for ( int i = 0; i < edges.size(); i++ ) {
  for ( int j = 0; j < edges.size(); j++ ) {
    if ( i == j )
      continue;
    bool res = comp( edges[i], edges[j] );
    if ( res != ( i < j ) ) {
      cout << "WRONG: Edges1: " << edges[i].source() << " --> " << edges[i].target() << " ";
      cout << "Edges2: " << edges[j].source() << " --> " << edges[j].target() << "    ";
      if ( res )
        cout << "Result:  1<2" <<endl;
      else
        cout << "Result:  2<1" <<endl;
    }
  }
}
#endif

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

    int hole_base = 0, hole_edge = 0;

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
        edges.push_back( Edge_type( query_pt, curr, on_edge, on_vertex ) );
        hole_edge++;
      } while ( ++curr != circ );
      for ( int i = 0; i < hole_edge; i++ ) {
        edges[hole_base+i].set_index( hole_base+i, hole_base+(i+hole_edge-1)%hole_edge );
      }
      hole_base += hole_edge;
    }
    for ( Hole_const_iterator hi = f->holes_begin();
          hi != f->holes_end(); hi++ ) {
      curr = circ = *hi;
      hole_edge = 0;
      do {
        bool on_edge = false, on_vertex = false;
        if ( curr == qedge ) {
          on_edge = true;
          on_vertex = query_on_vertex;
        } else if ( curr == qprev ) {
          on_edge = on_vertex = query_on_vertex;
        }
        edges.push_back( Edge_type( query_pt, curr, on_edge, on_vertex ) );
        hole_edge++;
      } while ( ++curr != circ );
      for ( int i = 0; i < hole_edge; i++ ) {
        edges[hole_base+i].set_index( hole_base+i, hole_base+(i+hole_edge-1)%hole_edge );
      }
      hole_base += hole_edge;
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

    int hole_base = 0, hole_edge = 0;

    if ( f->has_outer_ccb() ) {
      curr = circ = f->outer_ccb();
      do {
        edges.push_back( Edge_type( query_pt, curr, false, false ) );
        hole_edge++;
      } while ( ++curr != circ );
      for ( int i = 0; i < hole_edge; i++ ) {
        edges[hole_base+i].set_index( hole_base+i, hole_base+(i+hole_edge-1)%hole_edge );
      }
      hole_base += hole_edge;
    }
    for ( Hole_const_iterator hi = f->holes_begin();
          hi != f->holes_end(); hi++ ) {
      curr = circ = *hi;
      hole_edge = 0;
      do {
        edges.push_back( Edge_type( query_pt, curr, false, false ) );
        hole_edge++;
      } while ( ++curr != circ );
      for ( int i = 0; i < hole_edge; i++ ) {
        edges[hole_base+i].set_index( hole_base+i, hole_base+(i+hole_edge-1)%hole_edge );
      }
      hole_base += hole_edge;
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
