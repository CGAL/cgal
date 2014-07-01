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

  template < class E >
  class Less_Source
  {
  private:
    typedef E                                           Edge_type;
    typedef typename E::Geometry_traits_2               Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2         Point_2;
    typedef typename Geometry_traits_2::FT              Number_type;
  private:
    bool angle_less_than_pi( const Point_2& p )
    {
      CGAL::Orientation orient = CGAL::orientation( center,  aux, p );
      return ( orient == CGAL::LEFT_TURN ||
               ( orient == CGAL::COLLINEAR && p.x() >= center.x() ) );
    }
    bool less_general( const Edge_type& he1, const Edge_type& he2, CGAL::Orientation orient2 )
    {
      if ( angle_less_than_pi( he1.source() ) ) {
        if ( orient2 == CGAL::LEFT_TURN )
          return true;
        else 
          return ( CGAL::orientation( center, aux, he2.source() ) == CGAL::RIGHT_TURN );
      } else {
        if ( orient2 == CGAL::RIGHT_TURN )
          return false;
        else
          return ( CGAL::orientation( center, aux, he2.source() ) == CGAL::RIGHT_TURN );
      }  
    }
  public:
    Less_Source ( const Point_2& c )
      : center( c ), aux( c.x()+Number_type(1), c.y() )
      {}
    bool operator () ( const Edge_type& he1, const Edge_type& he2 )
    {
      CGAL::Orientation orient = CGAL::orientation( center, he1.source(), he2.source() );
      if ( orient == CGAL::COLLINEAR ) {
        if ( he2.at_source() )
          return false;
        if ( he1.at_source() )
          return true;
        if ( CGAL::collinear_are_ordered_along_line( he1.source(), center, he2.source() ) ) {
          // center is between he1.source() and he2.source()
          return angle_less_than_pi( he1.source() );
        } else {
          // center, he1.source(), he2.source() are ordered along the line
          if ( CGAL::collinear_are_ordered_along_line( center, he1.source(), he2.source() ) ) {
            // he1.source() is closer than he2.source()
            return ( he1.orientation() == CGAL::RIGHT_TURN ||
                     he1.at_source() || he1.outward() );
          } else {
            // he2.source() is closer than he1.source()
            return !( he2.orientation() == CGAL::RIGHT_TURN ||
                      he2.at_source() || he2.outward() );
          }
        }
      } else {
        return less_general( he1, he2, orient );
      }
    }
  private:
    Point_2 center;
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
  bool do_intersect_ray ( const Point_2& query_pt, const Point_2& aux, const Edge_type& he )
  {
    if ( he.orientation() != CGAL::COLLINEAR ) {
      Ray_2 ray( query_pt, aux );
      Segment_2 seg( he.source(), he.target() );
      return CGAL::do_intersect( ray, seg );
    } else if ( he.inward() || he.outward() ) {
      if ( CGAL::orientation( query_pt, aux, he.source() ) != CGAL::COLLINEAR )
        return false;
      return !CGAL::collinear_are_ordered_along_line( aux, query_pt, he.source() );
    } else if ( he.at_source() ) {
      Direction_2 d1( Segment_2( he.source(), he.target() ) );
      Direction_2 d2( Segment_2( he.source(), he.prev_source() ) );
      Direction_2 d0( Segment_2( query_pt, aux ) );
      return ( d0 == d1 || d0.counterclockwise_in_between( d2, d1 ) );
    } else if ( he.at_target() ) {
      Direction_2 d1( Segment_2( he.target(), he.next_target() ) );
      Direction_2 d2( Segment_2( he.target(), he.source() ) );
      Direction_2 d0( Segment_2( query_pt, aux ) );
      return ( d0 == d2 || d0.counterclockwise_in_between( d2, d1 ) );
    } else {
      // he.in_edge() == true
      return ( CGAL::orientation( aux, query_pt, he.target() ) != CGAL::LEFT_TURN );
    }
  }
  Point_2 calculate_intersection( const Point_2& query_pt, const Point_2& aux, const Circulator& he )
  {
     Ray_2 ray ( query_pt, aux );
     Segment_2 seg ( he->source()->point(), he->target()->point() );
     CGAL::Object res = CGAL::intersection( ray, seg );
     const Point_2 * ipoint = CGAL::object_cast<Point_2>(&res);
     assert( ipoint );
     return *ipoint;
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
      aux = (first+1)->source();
    }

    active.assign( unsorted_edges.size(), false );
    for ( Edge_iterator eit = edges.begin(); eit != edges.end(); eit++ ) {
      if ( do_intersect_ray( query_pt, aux, *eit ) ) { 
        if ( eit->orientation() == CGAL::LEFT_TURN  &&
          CGAL::orientation( query_pt, aux, eit->source() ) == CGAL::RIGHT_TURN ) {
          active_edges.insert( std::make_pair( *eit, true ) );
          active[eit->index()] = true;
        }
        if ( eit->orientation() == CGAL::RIGHT_TURN  &&
          CGAL::orientation( query_pt, aux, eit->target() ) == CGAL::RIGHT_TURN ) {
          active_edges.insert( std::make_pair( *eit, true ) );
          active[eit->index()] = true;
        }
      }
    }

#ifdef MYDEBUG
cout << "Unsorted edges" << endl;
cout << "================================" << endl;
for ( int i = 0; i < unsorted_edges.size(); i++ ) {
  cout << "Edge: " << unsorted_edges[i].source() << " --> " << unsorted_edges[i].target() << "  idx = " << unsorted_edges[i].index() << endl;
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
          Point_2 u = calculate_intersection( query_pt, aux, active_edges.begin()->first.circulator() );
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
          Point_2 u = calculate_intersection( query_pt, aux, top );
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
    if ( first_intersection != last_intersection )
      out_points.push_back( first_intersection );
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
    // Sort halfedges with their source point by polar angle
    Less_Source<Edge_type> comp ( query_pt );
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
