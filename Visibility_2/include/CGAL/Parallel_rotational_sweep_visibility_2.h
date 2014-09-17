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
// Author(s):  Kan Huang <huangkandiy@gmail.com>
//             Ning Xu <longyin0904@gmail.com>
//

#ifndef CGAL_PARALLEL_ROTATIONAL_SWEEP_VISIBILITY_2_H
#define CGAL_PARALLEL_ROTATIONAL_SWEEP_VISIBILITY_2_H

#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/result_of.h>
#include <boost/unordered_set.hpp> 
#include <boost/unordered_map.hpp> 
#include <iostream>
using namespace std;

#ifdef CGAL_LINKED_WITH_TBB
  #include "tbb/parallel_sort.h"
  #include "tbb/parallel_for.h"
#endif

#ifndef NDEBUG
  #include <iostream>
  using std::cout;
  using std::ostream;
  using std::endl;
#endif

namespace CGAL {

template < class Arrangement_2_,
           class RegularizationCategory = CGAL::Tag_true,
           class ConcurrencyCategory = CGAL::Parallel_tag >
class Parallel_rotational_sweep_visibility_2
{

// Public types declaration
public:
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                        Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

  typedef typename Geometry_traits_2::Kernel            K;
  typedef typename Geometry_traits_2::Point_2           Point_2;
  typedef typename Geometry_traits_2::Ray_2             Ray_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Line_2            Line_2;
  typedef typename Geometry_traits_2::Vector_2          Vector_2;
  typedef typename Geometry_traits_2::Direction_2       Direction_2;
  typedef typename Geometry_traits_2::FT                Number_type;
  typedef typename Geometry_traits_2::Object_2          Object_2;

  typedef RegularizationCategory                             Regularization_category;
  typedef CGAL::Tag_true                                Supports_general_polygon_category;
  typedef CGAL::Tag_true                                Supports_simple_polygon_category;

//private data types declaration
  typedef Ccb_halfedge_const_circulator                 Circulator;
  typedef typename Arrangement_2::Hole_const_iterator   Hole_const_iterator;
  typedef std::vector<int>                              Index_vector;

// private nested classes
private:
  class Vertex
  {
  private:
    Vertex_const_handle vh;
    int quad;
    int first_ic;
    int last_ic;
    int alias_idx;
    int sorted_idx;
  private:
  public:
    Vertex ( const Vertex_const_handle& v )
      : vh( v ), first_ic( 0 ), last_ic( 0 ), alias_idx( -1 ), sorted_idx( -1 )
      {}
    void init_quadrant ( const Point_2& query_pt )
    {
      CGAL::Comparison_result cx = CGAL::compare_x( vh->point(), query_pt );
      CGAL::Comparison_result cy = CGAL::compare_y( vh->point(), query_pt );
      if ( cx == CGAL::LARGER )
        quad = ( cy != CGAL::SMALLER ) ? 1 : 4;
      else if ( cx == CGAL::SMALLER )
        quad = ( cy != CGAL::LARGER ) ? 3 : 2;
      else {
        if ( cy == CGAL::LARGER )
          quad = 2;
        else if ( cy == CGAL::SMALLER )
          quad = 4;
        else
          quad = 0;
      }
    }
    const Vertex_const_handle& handle () const
      { return vh; }
    int quadrant () const
      { return quad; }
    int first_incident () const
      { return first_ic; }
    int last_incident () const
      { return last_ic; }
    int alias_index () const
      { return alias_idx; }
    int sorted_index () const
      { return sorted_idx; }

    void set_incident_index( int f, int l )
      { first_ic = f; last_ic = l; }
    void set_alias_index( int idx )
      { alias_idx = idx; }
    void set_sorted_index( int idx )
      { sorted_idx = idx; }

    const Point_2& point () const
      { return vh->point(); }
#ifndef NDEBUG
    void trace ( ostream& os, int level ) const
    {
      for ( int i = 0; i < level; i++ )
        os << "  ";
      os << "point=[" << point()
         << "], quadrant=[" << quadrant()
         << "], first ic=[" << first_incident()
         << "], last ic=[" << last_incident()
         << "]" << endl;
      for ( int i = 0; i < level; i++ )
        os << "  ";
      os << "               alias=[" << alias_index()
         << "], sorted=[" << sorted_index()
         << "]" << endl;
    }
#endif
  };

  class Edge
  {
  private:
    int source_idx;
    int target_idx;
    enum { LTURN, RTURN, OUTWARD, INWARD, AT_SOURCE, AT_TARGET, IN_EDGE } mode;
  public:
    Edge ( int s, int t )
      : source_idx( s ), target_idx( t )
      {}
    void init_orientation ( const Point_2& query_pt,
                            const Point_2& source,
                            const Point_2& target )
    {
      CGAL::Orientation orient = CGAL::orientation( query_pt, source, target );
      if ( orient == CGAL::LEFT_TURN )
        mode = LTURN;
      else if ( orient == CGAL::RIGHT_TURN )
        mode = RTURN;
      else {
        if ( query_pt == source )
          mode = AT_SOURCE;
        else if ( query_pt == target )
          mode = AT_TARGET;
        else if ( CGAL::collinear_are_ordered_along_line
                        ( query_pt, source, target ) )
          mode = OUTWARD;
        else if ( CGAL::collinear_are_ordered_along_line
                        ( query_pt, target, source ) )
          mode = INWARD;
        else
          mode = IN_EDGE;
      }
    }

    int source_index () const
      { return source_idx; }
    int target_index () const
      { return target_idx; }
    void set_index ( int s, int t )
      { source_idx = s; target_idx = t; }
    bool left_turn () const
      { return ( mode == LTURN ); }
    bool right_turn () const
      { return ( mode == RTURN ); }
    bool collinear () const
      { return ( !left_turn() && !right_turn() ); }
    bool inward () const
      { return ( mode == INWARD ); }
    bool outward () const
      { return ( mode == OUTWARD ); }
    bool at_source () const
      { return ( mode == AT_SOURCE ); }
    bool at_target () const
      { return ( mode == AT_TARGET ); }
    bool in_edge () const
      { return ( mode == IN_EDGE ); }
    bool pass_query_pt () const
      { return ( at_source() || at_target() || in_edge() ); }
    CGAL::Orientation orientation() const
    {
      if ( left_turn() ) return CGAL::LEFT_TURN;
      if ( right_turn() ) return CGAL::RIGHT_TURN;
      return CGAL::COLLINEAR;
    }
    CGAL::Orientation reverse_orientation() const
    {
      if ( left_turn() ) return CGAL::RIGHT_TURN;
      if ( right_turn() ) return CGAL::LEFT_TURN;
      return CGAL::COLLINEAR;
    }

#ifndef NDEBUG
    void trace ( ostream& os, int level ) const
    {
      for ( int i = 0; i < level; i++ )
        os << "  ";
      os << "source_idx=[" << source_index()
         << "],target_idx=[" << target_index()
         << "],mode =[";
      switch( mode ) {
      case LTURN:
        os << "left turn";
        break;
      case RTURN:
        os << "right turn";
        break;
      case OUTWARD:
        os << "outward";
        break;
      case INWARD:
        os << "inward";
        break;
      case AT_SOURCE:
        os << "at source";
        break;
      case AT_TARGET:
        os << "at target";
        break;
      case IN_EDGE:
        os << "in edge";
        break;
      }
      os << "]" << endl;
    }
#endif
  };

  class Closer_edge : public std::binary_function<int, int, bool>
  {
  private:
    const Point_2& query_pt;
    const std::vector<Edge>& edges;
    const std::vector<Vertex>& vertices;
  private:
    int vtype( const Edge& e, bool incident_at_source ) const
    {
      if ( incident_at_source ) {
        if ( e.left_turn() ) return 1;
        if ( e.right_turn() ) return 2;
        if ( e.outward()|| e.at_source() || e.in_edge() ) return 3;
        return 0;
      } else {
        if ( e.left_turn() ) return 2;
        if ( e.right_turn() ) return 1;
        if ( e.inward() || e.at_target() || e.in_edge() ) return 3;
        return 0;
      }
    }
    const Point_2& source_point( const Edge& e ) const
      { return vertices[e.source_index()].point(); }
    const Point_2& target_point( const Edge& e ) const
      { return vertices[e.target_index()].point(); }
    bool less ( const Edge& e1, const Edge& e2 ) const
    {
      CGAL::Orientation orient1, orient2, orient3;
      if ( e1.source_index() == e2.source_index() ) {
        int vt1 = vtype( e1, true );
        int vt2 = vtype( e2, true );
        if ( vt1 != vt2 ) return ( vt1 < vt2 );
        orient1 = CGAL::orientation( source_point( e2 ),
                                     target_point( e2 ),
                                     target_point( e1 ) );
        return ( orient1 == e2.orientation() );
      } else if ( e1.target_index() == e2.source_index() ) {
        int vt1 = vtype( e1, false );
        int vt2 = vtype( e2, true );
        if ( vt1 != vt2 ) return ( vt1 < vt2 );
        orient1 = CGAL::orientation( source_point( e2 ),
                                     target_point( e2 ),
                                     source_point( e1 ) );
        return ( orient1 == e2.orientation() );
      } else if ( e1.source_index() == e2.target_index() ) {
        int vt1 = vtype( e1, true );
        int vt2 = vtype( e2, false );
        if ( vt1 != vt2 ) return ( vt1 < vt2 );
        orient1 = CGAL::orientation( target_point( e2 ),
                                     source_point( e2 ),
                                     target_point( e1 ) );
        return ( orient1 == e2.reverse_orientation() );
      } else if ( e1.target_index() == e2.target_index() ) {
        int vt1 = vtype( e1, false );
        int vt2 = vtype( e2, false );
        if ( vt1 != vt2 ) return ( vt1 < vt2 );
        orient1 = CGAL::orientation( target_point( e2 ),
                                     source_point( e2 ),
                                     source_point( e1 ) );
        return ( orient1 == e2.reverse_orientation() );
      } else {
        // General case
        const Point_2& s1 = source_point( e1 );
        const Point_2& t1 = target_point( e1 );
        const Point_2& s2 = source_point( e2 );
        const Point_2& t2 = target_point( e2 );

        if ( e1.collinear() ) {
          if ( e2.collinear() )
            return CGAL::collinear_are_ordered_along_line( query_pt, s1, s2 );
          else
            return ( CGAL::orientation( s2, t2, s1 ) == e2.orientation() );
        } else {
          orient1 = CGAL::orientation( s1, t1, s2 );
          orient2 = CGAL::orientation( s1, t1, t2 );
          if ( orient1 == CGAL::COLLINEAR )
            return ( orient2 != e1.orientation() );
          if ( orient1 == e1.orientation() ) {
            if ( orient2 != e1.orientation() )
              return ( CGAL::orientation( s2, t2, s1 ) == e2.orientation() );
            else
              return false;
          } else {
            if ( orient2 != e1.orientation() )
              return true;
            else
              return ( CGAL::orientation( s2, t2, s1 ) == e2.orientation() );
          }
        }
      }
    }
  public:
    Closer_edge( const Point_2& q, const std::vector<Edge>& e,
                 const std::vector<Vertex>& v )
      : query_pt( q ), edges( e ), vertices( v )
      {}
    bool operator () ( int idx1, int idx2 ) const
    {
      if ( idx1 == idx2 ) return false;
      if ( edges[idx1].pass_query_pt() && edges[idx2].pass_query_pt() )
        return ( idx1 < idx2 );
      if ( edges[idx1].pass_query_pt() )
        return true;
      if ( edges[idx2].pass_query_pt() )
        return false;
      return less( edges[idx1], edges[idx2] );
    }
  };
  
  class Intersection_edges
  {
  private:
    struct IE_node {
      bool exist;
      int prev;
      int next;
      IE_node ()
        : exist( false ), prev( -1 ), next( -1 )
        {}
      IE_node ( bool f, int p, int n )
        : exist( f ), prev( p ), next( n )
        {}
    };
  private:
    std::vector< IE_node > _data;
    int _head;
  public:
    Intersection_edges( int s )
    {
      _data.reserve( s+1 );
      for ( int i = 0; i < s+1; i++ )
        _data.push_back( IE_node( false, i, i ) );
      _head = s;
    }
    void insert( int idx )
    {
      if ( _data[idx].exist )
        return;
      _data[idx].exist = true;
      _data[idx].next = _data[_head].next;
      _data[idx].prev = _head;
      _data[_data[_head].next].prev = idx;
      _data[_head].next = idx;
    }
    void erase ( int idx )
    {
      if ( !_data[idx].exist )
        return;
      _data[_data[idx].prev].next = _data[idx].next;
      _data[_data[idx].next].prev = _data[idx].prev;
      _data[idx].exist = false;
      _data[idx].prev = _data[idx].next = idx;
    }
    std::vector<int> data ()
    {
      std::vector<int> res;
      int curr = _data[_head].next;
      while ( curr != _head ) {
        res.push_back( curr );
        curr = _data[curr].next;
      }
      return res;
    }
    int begin () const
      { return _data[_head].next; }
    int end () const
      { return _head; }
    int next ( int curr ) const
      { return _data[curr].next; }
    bool has ( int idx ) const
    { return _data[idx].exist; }
  };

  class Cone
  {
  private:
    int _begin;
    int _end;
    std::vector<int> _edx;
  public:
    Cone ( int b, int e )
      : _begin( b ), _end( e )
      { _edx.reserve( e - b ); }
    int begin () const
      { return _begin; }
    int end () const
      { return _end; }
    int size () const
      { return _edx.size(); }
    int operator [] ( int pos ) const
      { return _edx[pos]; }
    void insert ( int idx )
      { _edx.push_back( idx ); }
#ifndef NDEBUG
    void trace ( ostream& os, int level ) const
    {
      for ( int i = 0; i < level; i++ ) os << "  ";
      os << "begin = " << begin() << " , end = " << end() << " , size = " << size() << endl;
      for ( int i = 0; i < level; i++ ) os << "  ";
      os << "edx = [ ";
      for ( int i = 0; i < _edx.size(); i++ )
        os << _edx[i] << " ";
      os << "]" << endl;
    }
#endif
  };

  class Sub_region
  {
  private:
    struct Intersection_point
    {
      bool by_edge;
      int v_idx;
      int e_idx;
      Intersection_point( bool b, int v, int e )
        : by_edge( b ), v_idx( v ), e_idx( e )
        {}
    };
  private:
    std::vector<Intersection_point> _pts;
  private:
    Point_2 compute_intersection( const Point_2& q,
                                  const Point_2& dp,
                                  const Point_2& s,
                                  const Point_2& t ) const
   {
     if ( CGAL::collinear( q, dp, s ) ) {
       if ( CGAL::collinear( q, dp, t ) ) {
         if ( CGAL::collinear_are_ordered_along_line( q, s, t ) )
           return s;
         else
           return t;
       } else {
         return s;
       }
     }
     Ray_2 ray( q, dp );
     Segment_2 seg( s, t );
     CGAL::Object obj = CGAL::intersection( ray, seg );
     const Point_2 * ipoint = CGAL::object_cast<Point_2> (&obj);
     assert( ipoint != NULL );
     return *ipoint;
   }
  public:
    Sub_region ( int est_pts_num )
      { _pts.reserve( est_pts_num*2+2 ); }
    int size () const
      { return _pts.size(); }
    Point_2 get_point ( const Point_2& q,
                        const std::vector<Vertex>& vertices,
                        const std::vector<Edge>& edges,
                        int pos ) const
    {
      if ( _pts[pos].by_edge ) {
        return compute_intersection( q, vertices[_pts[pos].v_idx].point(),
                      vertices[edges[_pts[pos].e_idx].source_index()].point(),
                      vertices[edges[_pts[pos].e_idx].target_index()].point() );
      } else {
        return vertices[_pts[pos].v_idx].point();
      }
    }
    void insert ( int v_idx )
      { _pts.push_back( Intersection_point( false, v_idx, -1 ) ); }
    void insert ( int v_idx, int e_idx )
      { _pts.push_back( Intersection_point( true, v_idx, e_idx ) ); }
#ifndef NDEBUG
    void trace ( ostream& os, int level ) const
    {
      for ( int i = 0; i < level; i++ ) os << "  ";
      os << "size = " << size() << endl;
      for ( int i = 0; i < level; i++ ) os << "  ";
      os << "pts = [ ";
      for ( int i = 0; i < _pts.size(); i++ )
      os << _pts[i].by_edge << " " << _pts[i].v_idx << " " <<_pts[i].e_idx << " , ";
      os << "]" << endl;
    }
#endif
  };

//private data types declaration
private:
  typedef Vertex                                        Vertex_type;
  typedef std::vector<Vertex_type>                      Vertex_vector;
  typedef Edge                                          Edge_type;
  typedef std::vector<Edge_type>                        Edge_vector;
  typedef std::vector<Point_2>                          Point_vector;
  
private:
  class Is_swept_earlier
  {
  private:
    const Point_2& query_pt;
    const Vertex_vector * vertices;
    const Point_2& shifted_source;
    int shifted_quadrant;

    bool is_vertex_query;
  private:
    // less_vertex
    // Precondition: v1 == query_pt, v2 != query_pt
    bool less_vertex ( const Vertex_type& v1, const Vertex_type& v2 ) const
    {
      CGAL::Orientation orient = CGAL::orientation( query_pt,
                                                    shifted_source,
                                                    v2.point() );
      if ( orient == CGAL::COLLINEAR ) {
        if ( shifted_quadrant != v2.quadrant() )
          return ( shifted_quadrant < v2.quadrant() );
        else
          return true;
      } else {
        if ( shifted_quadrant != v2.quadrant() )
          return ( shifted_quadrant < v2.quadrant() );
        else
          return ( orient == CGAL::LEFT_TURN );
      }
    }
    bool less_general ( const Vertex_type& v1, const Vertex_type& v2 ) const
    {
      if ( v1.quadrant() != v2.quadrant() )
        return ( v1.quadrant() < v2.quadrant() );
      CGAL::Orientation orient = CGAL::orientation( query_pt,
                                 v1.point(), v2.point() );
      if ( orient != CGAL::COLLINEAR )
        return ( orient == CGAL::LEFT_TURN );
      else
        return CGAL::collinear_are_ordered_along_line( query_pt,
                                 v1.point(), v2.point() );
    }
    // less
    // compare two points by their polar angle
    bool less ( const Vertex_type& v1, const Vertex_type& v2 ) const
    {
      if ( v1.handle() == v2.handle() )
        return false;
      if ( is_vertex_query ) {
        if ( v1.quadrant() == 0 )
          return less_vertex( v1, v2 );
        else if ( v2.quadrant() == 0 )
          return !less_vertex( v2, v1 );
        else
          return less_general( v1, v2 );
      } else
        return less_general( v1, v2 );
    }
  public:
    Is_swept_earlier ( const Point_2& q, const Vertex_vector * pv,
                       const Point_2& fs, int fq )
      : query_pt( q ), vertices( pv ), shifted_source( fs ),
        shifted_quadrant( fq )
      { is_vertex_query = ( fs != query_pt ); }
    bool operator () ( int vdx1, int vdx2 ) const
      { return less( vertices->at( vdx1 ), vertices->at( vdx2 ) ); }
  };

#ifdef CGAL_LINKED_WITH_TBB
  class Parallel_sweep
  {
  private:
    Parallel_rotational_sweep_visibility_2 * algo;
  public:
    Parallel_sweep ( Parallel_rotational_sweep_visibility_2 * a  )
      : algo( a )
      {}
    void operator () ( const tbb::blocked_range<int>& range ) const
    {
      for ( int i = range.begin(); i != range.end(); i++ )
        algo->compute_visibility_partition( i );
    }
  };

  class Parallel_init_quadrant
  {
  private:
    Parallel_rotational_sweep_visibility_2 * algo;
  public:
    Parallel_init_quadrant ( Parallel_rotational_sweep_visibility_2 * a  )
      : algo( a )
      {}
    void operator () ( const tbb::blocked_range<int>& range ) const
    {
      for ( int i = range.begin(); i != range.end(); i++ )
        algo->init_quadrant_each_vertex( i );
    }
  };

  class Parallel_init_orientation
  {
  private:
    Parallel_rotational_sweep_visibility_2 * algo;
  public:
    Parallel_init_orientation ( Parallel_rotational_sweep_visibility_2 * a  )
      : algo( a )
      {}
    void operator () ( const tbb::blocked_range<int>& range ) const
    {
      for ( int i = range.begin(); i != range.end(); i++ )
        algo->init_orientation_each_edge( i );
    }
  };

  class Parallel_do_intersect_edge
  {
  private:
    Parallel_rotational_sweep_visibility_2 * algo;
    const Point_2& dp;
    std::vector<int>& results;
  public:
    Parallel_do_intersect_edge ( Parallel_rotational_sweep_visibility_2 * a,
                                 const Point_2& p,
                                 std::vector<int>& r  )
      : algo( a ), dp( p ), results( r )
      {}
    void operator () ( const tbb::blocked_range<int>& range ) const
    {
      for ( int i = range.begin(); i != range.end(); i++ )
        results[i] = algo->do_intersect_edge( dp, i );
    }
  };
#endif

//private methods
private:
  void add_box ()
  {
    std::vector<Point_2> pts;
    pts.reserve( vs.size()+1 );
    for ( int i = 0; i < vs.size(); i++ )
      pts.push_back( vs[i].point() );
    pts.push_back( query_pt );

    typename Geometry_traits_2::Iso_rectangle_2 bb = CGAL::bounding_box( pts.begin(), pts.end() );

    Number_type xmin = bb.xmin() - 1,
                ymin = bb.ymin() - 1,
                xmax = bb.xmax() + 1,
                ymax = bb.ymax() + 1;
    Point_2 box[4] = { Point_2( xmin, ymin ), Point_2( xmax, ymin ),
                       Point_2( xmax, ymax ), Point_2( xmin, ymax ) };
    Halfedge_handle e1 = arr_box.insert_in_face_interior
                ( Segment_2( box[0], box[1] ), arr_box.unbounded_face() );
    Halfedge_handle e2 = arr_box.insert_from_left_vertex
                ( Segment_2( box[1], box[2] ), e1->target() );
    Halfedge_handle e3 = arr_box.insert_from_right_vertex
                ( Segment_2( box[2], box[3] ), e2->target() );
    arr_box.insert_at_vertices( Segment_2( box[0], box[3] ),
                                e1->source(), e3->target() );

    Circulator circ,curr;
    circ = curr = e1->face()->outer_ccb();
    int edge_base = vs.size();
    do {
      assert( curr->face() == e1->face() );
      vs.push_back( Vertex_type( curr->source() ) );
      es.push_back( Edge_type( vs.size()-1, vs.size() ) );
    } while ( ++curr != circ );
    es.back().set_index( es.back().source_index(), edge_base );
  }

  void compute_shifted_source ()
  {
    if ( query_type != VERTEX_QUERY ) {
      shifted_source = query_pt;
      shifted_quadrant = 0;
    } else {
      Vertex_type vt( query_edge->next()->target() );
      vt.init_quadrant( query_pt );
      if ( is_small_cone ) {
        shifted_source = Point_2( query_pt.x() + query_pt.x() - vt.point().x(),
                                  query_pt.y() + query_pt.y() - vt.point().y() );
        shifted_quadrant = ( vt.quadrant() + 2 ) % 4;
      } else {
        shifted_source = vt.point();
        shifted_quadrant = vt.quadrant();
      }
    }
  }

  void init_quadrant_each_vertex( int i )
    { vs[i].init_quadrant( query_pt ); }
  void init_quadrant_parallel( CGAL::Sequential_tag )
  {
    for ( int i = 0; i < vs.size(); i++ )
      init_quadrant_each_vertex( i );
  }
  void init_quadrant_parallel( CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_EBB
    Parallel_init_quadrant init( this );
    tbb::parallel_for( tbb::blocked_range<int>(0, vs.size() ), init );
#else
    init_quadrant_parallel( CGAL::Sequential_tag() );
#endif
  }

  void init_orientation_each_edge( int i )
  { 
    es[i].init_orientation( query_pt, vs[es[i].source_index()].point(),
                            vs[es[i].target_index()].point() );
  }
  void init_orientation_parallel( CGAL::Sequential_tag )
  {
    for ( int i = 0; i < es.size(); i++ )
      init_orientation_each_edge( i );
  }
  void init_orientation_parallel( CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_EBB
    Parallel_init_orientation init( this );
    tbb::parallel_for( tbb::blocked_range<int>(0, es.size() ), init );
#else
    init_orientation_parallel( CGAL::Sequential_tag() );
#endif
  }

  void sort_vertices( CGAL::Sequential_tag )
  {
    Is_swept_earlier comp ( query_pt, &vs, shifted_source, shifted_quadrant );
    std::sort( good_vdx.begin(), good_vdx.end(), comp );
  }

  void sort_vertices( CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_TBB
    Is_swept_earlier comp ( query_pt, &vs, shifted_source, shifted_quadrant );
    tbb::parallel_sort( good_vdx.begin(), good_vdx.end(), comp );
#else
    sort_vertices( CGAL::Sequential_tag() );
#endif
  }

  void remove_duplicated_vertices()
  {
    // find duplicated vertices
    int last = 0;
    for ( int i = 0; i < good_vdx.size(); i++ ) {
      if ( vs[good_vdx[i]].handle() != vs[good_vdx[last]].handle() ) {
        last++;
        vs[good_vdx[i]].set_sorted_index( last );
        good_vdx[last] = good_vdx[i];
      } else {
        vs[good_vdx[i]].set_sorted_index( last );
        vs[good_vdx[i]].set_alias_index( good_vdx[last] );
      }
    }
    good_vdx.erase( good_vdx.begin()+last+1, good_vdx.end() );
  }

  void sort_incident ( CGAL::Sequential_tag )
    { std::sort( incident.begin(), incident.end() ); }
  void sort_incident ( CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_TBB
    tbb::parallel_sort( incident.begin(), incident.end() );
#else
    sort_incident( CGAL::Sequential_tag() );
#endif
  }

  void construct_incident_map()
  {
    incident.clear();
    incident.reserve( es.size()*2 );
    for ( int i = 0; i < es.size(); i++ ) {
      incident.push_back( std::make_pair( vs[es[i].source_index()].alias_index(), i ) );
      incident.push_back( std::make_pair( vs[es[i].target_index()].alias_index(), i ) );
    }
    sort_incident( ConcurrencyCategory() );
    
    int i = 0;
    while ( i < incident.size() ) {
      int j = i;
      while ( j < incident.size() && incident[j].first == incident[i].first )
        j++;
      vs[incident[i].first].set_incident_index( i, j );
      i = j;
    }
  }

  bool funnel_block_right( int v_idx, int e_idx )
  {
    int s_idx = vs[es[e_idx].source_index()].alias_index();
    int t_idx = vs[es[e_idx].target_index()].alias_index();
    assert( s_idx == v_idx || t_idx == v_idx );
    if ( v_idx == s_idx ) {
      return ( es[e_idx].right_turn() ||
               es[e_idx].outward() ||
               ( es[e_idx].at_source() && !is_small_cone ) );
    } else {
      return ( es[e_idx].left_turn() ||
               es[e_idx].outward() ||
               es[e_idx].in_edge() ||
               es[e_idx].at_source() ||
               ( es[e_idx].at_target() && !is_small_cone ) );
    }
  }

  bool funnel_has_precedessor( int v_idx, int e_idx )
  {
    int s_idx = vs[es[e_idx].source_index()].alias_index();
    int t_idx = vs[es[e_idx].target_index()].alias_index();
    assert( s_idx == v_idx || t_idx == v_idx );
    if ( v_idx == s_idx )
      return es[e_idx].inward();
    else
      return es[e_idx].outward();
  }

  void funnel ( int first, int last )
  {
    std::vector<int> left, right;
    left.reserve( last - first );
    right.reserve( last - first );
    bool block_right = false;
    for ( int i = first; i < last; i++ ) {
      int v_idx = good_vdx[i];
      bool right_v = false, has_precedessor = false;
      for ( int j = vs[v_idx].first_incident();
                j < vs[v_idx].last_incident(); j++ ) {
        int e_idx = incident[j].second;
        if ( funnel_block_right( v_idx, e_idx ) )
          right_v = true;
        if ( funnel_has_precedessor( v_idx, e_idx ) )
          has_precedessor = true;
      }
      if ( has_precedessor )
        block_right = block_right || right_v;
      else
        block_right = right_v;
      if ( block_right )
        right.push_back( good_vdx[i] );
      else
        left.push_back( good_vdx[i] );
    }
    for ( int i = 0; i < right.size(); i++ )
      good_vdx[first+i] = right[i];
    for ( int i = 0; i < left.size(); i++ )
      good_vdx[first+right.size()+i] = left[left.size()-1-i];
    for ( int i = first; i < last; i++ )
      vs[good_vdx[i]].set_sorted_index( i );
  }

  void process_funnel ()
  {
    // TBD: future inmprovement: parallelly compute orientation
    std::vector<CGAL::Orientation> orients( good_vdx.size()-1, CGAL::LEFT_TURN );
    for ( int i = 0; i < good_vdx.size()-1; i++ ) {
      if ( query_type == VERTEX_QUERY ) {
        if ( vs[good_vdx[i+1]].quadrant() == 0 ) {
          // (i+1)-th point is the query point
          orients[i] = CGAL::orientation( query_pt, vs[good_vdx[i]].point(),
                                          shifted_source );
          if ( orients[i] == CGAL::COLLINEAR &&
               shifted_quadrant != vs[good_vdx[i]].quadrant() )
            orients[i] = CGAL::RIGHT_TURN;
        } else if ( vs[good_vdx[i]].quadrant() == 0 ) {
          // i-th point is the query point
          orients[i] = CGAL::orientation( query_pt, shifted_source,
                                          vs[good_vdx[i+1]].point() );
          if ( orients[i] == CGAL::COLLINEAR &&
               shifted_quadrant != vs[good_vdx[i+1]].quadrant() )
            orients[i] = CGAL::RIGHT_TURN;
        } else {
          orients[i] = CGAL::orientation( query_pt, vs[good_vdx[i]].point(),
                                          vs[good_vdx[i+1]].point() );
        }
      } else {
        orients[i] = CGAL::orientation( query_pt, vs[good_vdx[i]].point(),
                                          vs[good_vdx[i+1]].point() );
        if ( orients[i] == CGAL::COLLINEAR &&
             vs[good_vdx[i]].quadrant() != vs[good_vdx[i+1]].quadrant() )
          orients[i] = CGAL::RIGHT_TURN;
      }
    }

    for ( int i = 0; i < good_vdx.size(); i++ ) {
      int j = i + 1;
      while ( j < good_vdx.size() ) {
        if ( orients[j-1] != CGAL::COLLINEAR )
          break;
        j++;
      }
      if ( j - i > 1 )
        funnel ( i, j );
      i = j - 1;
    }
  }

#ifndef NDEBUG
  void check_consistency_after_init()
  {
    for ( int i = 0; i < vs.size(); i++ ) {
      int alias = vs[i].alias_index();
      if ( vs[alias].handle() != vs[i].handle() ) {
        cout << "vs[" << i << "].handle() != vs[" << alias << "].handle()" << endl;
      }
      int vdx = vs[alias].sorted_index();
      if ( vdx < 0 || vdx >= good_vdx.size() ) {
        cout << "vs[" << alias << "].sorted_index() = " << vdx << " : out of range" << endl;
      } else if ( good_vdx[vdx] != alias ) {
        cout << "Inconsistency: vs[" << alias << "].sorted_index = " << vdx << ", but good_vdx[" << vdx  <<"] = " << good_vdx[vdx] << endl;
      }
      for ( int j = vs[alias].first_incident(); j < vs[alias].last_incident(); j++ ) {
        if ( incident[j].first != alias ) {
          cout << "Inconsistency: in vs[" << alias << "].incident() : incident[" << j << "].first = " << incident[j].first << endl;
        }
        if ( j == vs[alias].first_incident() &&
             j != 0 && incident[j-1].first >= alias ) {
          cout << "Inconsistency: in vs[" << alias << "].incident() : previous : incident[" << (j-1) << "].first = " << incident[j-1].first << endl;
        }
        if ( j == vs[alias].last_incident()-1 &&
             j < incident.size()-1 && incident[j+1].first <= alias ) {
          cout << "Inconsistency: in vs[" << alias << "].incident() : next : incident[" << (j+1) << "].first = " << incident[j+1].first << endl;
        }
        int s_idx = es[incident[j].second].source_index();
        int t_idx = es[incident[j].second].target_index();
        if ( vs[s_idx].alias_index() != alias &&
             vs[t_idx].alias_index() != alias  ) {
          cout << "Inconsistency: in vs[" << i << "] : alias = [" << alias << "]; incident[" << j << "].second=[" << incident[j].second << "]" << endl;
        }
      }
    }
  }
#endif

  void keep_consistency_after_init()
  {
    for ( int i = 0; i < vs.size(); i++ ) {
      int alias = vs[i].alias_index();
      if ( alias == i )
        continue;
      vs[i].set_incident_index( vs[alias].first_incident(),
                             vs[alias].last_incident() );
      vs[i].set_sorted_index( vs[alias].sorted_index() );
    }
    for ( int i = 0; i < es.size(); i++ ) {
      int s_idx = es[i].source_index();
      int t_idx = es[i].target_index();
      es[i].set_index( vs[s_idx].alias_index(), vs[t_idx].alias_index() );
    }
  }

  void init_vertices ( const Face_const_handle& fh )
  {
    Circulator circ, curr;
    Hole_const_iterator hi;

    vs.clear();
    es.clear();
    good_vdx.clear();
    arr_box.clear();

    curr = circ = fh->outer_ccb();
    int edge_base = vs.size();
    do {
      assert( curr->face() == fh );
      vs.push_back( Vertex_type( curr->source() ) );
      es.push_back( Edge_type( vs.size()-1, vs.size() ) );
    } while ( ++curr != circ );
    es.back().set_index( es.back().source_index(), edge_base );
    for ( hi = fh->holes_begin(); hi != fh->holes_end(); hi++ ) {
      curr = circ = *hi;
      edge_base = vs.size();
      do {
        assert( curr->face() == fh );
        vs.push_back( Vertex_type( curr->source() ) );
        es.push_back( Edge_type( vs.size()-1, vs.size() ) );
      } while ( ++curr != circ );
      es.back().set_index( es.back().source_index(), edge_base );
    }

    if ( query_type != FACE_QUERY )
      add_box();

    init_quadrant_parallel( ConcurrencyCategory() );

    init_orientation_parallel( ConcurrencyCategory() );

    good_vdx.reserve( vs.size() );
    for ( int i = 0; i < vs.size(); i++ )
      good_vdx.push_back( i );

    // sort vertices by their polar angle
    compute_shifted_source();
    sort_vertices( ConcurrencyCategory() );

    // Build the reverse indexes
    for ( int i = 0; i < good_vdx.size(); i++ ) {
      vs[good_vdx[i]].set_alias_index( good_vdx[i] );
      vs[good_vdx[i]].set_sorted_index( i );
    }
     
    // Remove duplicated vertices
    remove_duplicated_vertices();

    // construct the edge index list
    construct_incident_map();

    // deal with funnel
    process_funnel();

    // Keep consistency
    keep_consistency_after_init();
  }

  // Precondtion: dp != any end point of the edge.
  int do_intersect_edge ( const Point_2& dp, int i )
  {
    CGAL::Orientation orient1, orient2;
    if ( es[i].pass_query_pt() )    // ignore bad edges
      return 0;
    if ( es[i].inward() || es[i].outward() )
      return 0;
    orient1 = CGAL::orientation( query_pt, dp,
                                 vs[es[i].source_index()].point() );
    orient2 = CGAL::orientation( query_pt, dp,
                                 vs[es[i].target_index()].point() );
    if ( orient1 != orient2 && orient1 != es[i].orientation() )
      return 1;
    return 0;
  }
  void do_intersect_parallel ( const Point_2& dp,
                               std::vector<int>& results,
                               CGAL::Sequential_tag )
  {
    for ( int i = 0; i < es.size(); i++ )
      results[i] = do_intersect_edge( dp, i );
  }
  void do_intersect_parallel ( const Point_2& dp,
                               std::vector<int>& results,
                               CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_TBB
    Parallel_do_intersect_edge obj( this, dp, results );
    tbb::parallel_for( tbb::blocked_range<int>(0, es.size() ), obj );
#else
    do_intersect_parallel( dp, results, CGAL::Sequential_tag() );
#endif
  }

  int default_cone_size ( CGAL::Sequential_tag )
    { return vs.size(); }
  int default_cone_size ( CGAL::Parallel_tag )
    { return 256; }

  void partition_cones ()
  {
    Intersection_edges active_edges( es.size() );
    int curr;

    cones.clear();

    // Determine the first ray
    Point_2 dp;
    if ( vs[good_vdx.back()].quadrant() == 0 ) {
      if ( shifted_quadrant < 4 )
        dp = query_pt + Vector_2( 1, -1 );
      else
        dp = shifted_source + Vector_2( 1, 0 );
    } else if ( vs[good_vdx.back()].quadrant() < 4 ) {
      dp = query_pt + Vector_2( 1, -1 );
    } else {
      dp = vs[good_vdx.back()].point() + Vector_2( 1, 0 );
    }

    // Initialize
    std::vector<int> is( es.size(), 0 );
    /*
     * There is a bug here. When I verified the intersection parallelly,
     * the program crashed. Thus, I use the Sequential_tag here.
     * I will solve this one later.
    do_intersect_parallel( dp, is, ConcurrencyCategory() );
    */
    do_intersect_parallel( dp, is, CGAL::Sequential_tag() );
    for ( int i = 0; i < es.size(); i++ ) {
      if ( is[i] != 0 )
        active_edges.insert( i );
    }

    // initialize the first cone
    int step = default_cone_size( ConcurrencyCategory() );
    int cone_base = 0;
    int cone_next = cone_base + step;
    if ( cone_next + 16 >= vs.size() )
      cone_next = vs.size();
    cones.push_back( Cone( cone_base, cone_next ) );
    for ( curr = active_edges.begin(); curr != active_edges.end();
          curr = active_edges.next( curr ) ) {
      cones.back().insert( curr );
    }
    if ( cone_next == vs.size() )
      return;

    // Rotational sweep points to find the intersection edges for other cones
    for ( int i = 0; i < good_vdx.size(); i++ ) {
      if ( i == cone_next ) {
        cone_base = cone_next;
        cone_next = cone_base + step;
        if ( cone_next + 16 >= good_vdx.size() )
          cone_next = good_vdx.size();
        cones.push_back( Cone( cone_base, cone_next ) );
        for ( curr = active_edges.begin(); curr != active_edges.end();
              curr = active_edges.next( curr ) ) {
          cones.back().insert( curr );
        }
        if ( cone_next == good_vdx.size() )
          break;
      }

      int v_idx = good_vdx[i];
      if ( vs[v_idx].quadrant() == 0 )
        continue;   // ignore bad vertex
      for ( int j = vs[v_idx].first_incident();
                j < vs[v_idx].last_incident(); j++ ) {
        int edx = incident[j].second;
        if ( es[edx].pass_query_pt() )
          continue;   // ignore bad edges
        if ( !active_edges.has( edx ) ) {
          active_edges.insert( edx );
        } else {
          active_edges.erase( edx );
        }
      }
    }
  }

  Point_2 ray_edge_intersection( int v_idx, int e_idx )
  {
    const Point_2& dp = vs[v_idx].point();
    const Point_2& s = vs[es[e_idx].source_index()].point();
    const Point_2& t = vs[es[e_idx].target_index()].point();

    if ( CGAL::collinear( query_pt, dp, s ) ) {
      if ( CGAL::collinear( query_pt, dp, t ) ) {
        if ( CGAL::collinear_are_ordered_along_line( query_pt, s, t ) )
          return s;
        else
          return t;
      } else {
        return s;
      }
    }

    Ray_2 ray( query_pt, dp );
    Segment_2 seg( s, t );
    CGAL::Object obj = CGAL::intersection( ray, seg );
    const Point_2 * ipoint = CGAL::object_cast<Point_2> (&obj );
    assert( ipoint != NULL );
    return Point_2( ipoint->x(), ipoint->y() );
  }

  void compute_visibility_partition( int cone_idx )
  {
    const Cone& cone = cones[cone_idx];
    Sub_region& result = sub_regions[cone_idx];
    typedef std::set<int,Closer_edge> edge_container;
    typedef typename edge_container::const_iterator edge_iterator;

    Closer_edge comp( query_pt, es, vs );
    edge_container active_edges( comp );
    edge_iterator eit;

    // Initialize
    for ( int i = 0; i < cone.size(); i++ )
      active_edges.insert( cone[i] );

    // Rotational sweep
    std::vector<int> insert_edx;
    std::vector<edge_iterator> remove_its;
    for ( int i = cone.begin(); i < cone.end(); i++ ) {
      int old_top, new_top;
      if ( active_edges.empty() )
        old_top = -1;
      else
        old_top = *active_edges.begin();
      int v_idx = good_vdx[i];
      if ( vs[v_idx].quadrant() == 0 )
        continue;   // ignore bad vertex

      insert_edx.clear();
      remove_its.clear();
      for ( int j = vs[v_idx].first_incident();
                j < vs[v_idx].last_incident(); j++ ) {
        int edx = incident[j].second;
        if ( es[edx].pass_query_pt() )
          continue;   // ignore bad edges
        eit = active_edges.find( edx );
        if ( eit == active_edges.end() )
          insert_edx.push_back( edx );
        else
          remove_its.push_back( eit );
      }

      if ( insert_edx.empty() && remove_its.empty() )
        continue;

      if ( insert_edx.size() == 1 && remove_its.size() == 1 ) {
        // Simple replacement
        // Use a dirty trick for local replacement
        const int& ctemp = *remove_its.front();
        int& temp = const_cast<int&>( ctemp );
        temp = insert_edx.front();
      } else {
        for ( int i = 0; i < insert_edx.size(); i++ ) {
          active_edges.insert( insert_edx[i] );
        }
        for ( int i = 0; i < remove_its.size(); i++ )
          active_edges.erase( remove_its[i] );
      }

      if ( active_edges.empty() )
        new_top = -1;
      else
        new_top = *active_edges.begin();
      if ( old_top != new_top ) {
        // The closest edge changed
        if ( !insert_edx.empty() && !remove_its.empty() ) {
          // Some edges are added, and some are removed
          // the current vertex is part of visibility region
          result.insert( v_idx );
        } else if ( remove_its.empty() ) {
          // Only add edges, the view ray is blocked by new edges
          if ( old_top != -1 )
            result.insert( v_idx, old_top );
          result.insert( v_idx );
        } else {
          // Only remove edges, a block for the view ray is removed
          result.insert( v_idx );
          if ( new_top != -1 )
            result.insert( v_idx, new_top );
        }
      }
    }
  }

  void compute_visibility_parallel( CGAL::Sequential_tag )
  {
    for ( int i = 0; i < cones.size(); i++ )
      compute_visibility_partition( i );
  }

  void compute_visibility_parallel( CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_TBB
    if ( cones.size() == 1 )
      return compute_visibility_parallel( CGAL::Sequential_tag() );
    Parallel_sweep sweep( this );
    tbb::parallel_for( tbb::blocked_range<int>(0, cones.size() ), sweep );
#else
    return compute_visibility_parallel( CGAL::Sequential_tag() );
#endif
  }

  void merge_result()
  {
    polygon.clear();
    cone_end_idx = cone_start_idx = -1;
    for ( int i = 0; i < sub_regions.size(); i++ ) {
      for ( int j = 0; j < sub_regions[i].size(); j++ ) {
        Point_2 p = sub_regions[i].get_point( query_pt, vs, es, j );
        if ( !polygon.empty() && p == polygon.back() )
          continue;
        polygon.push_back( p );
        if ( query_type != FACE_QUERY ) {
          if ( p == cone_end->point() )
            cone_end_idx = polygon.size() - 1;
          if ( p == cone_start->point() )
            cone_start_idx = polygon.size() - 1;
        }
      }
    }
    assert( polygon.size() > 2 );
  }

  void compute_visibility_impl ( const Face_const_handle& fh )
  {
    assert( !fh->is_unbounded() );

    init_vertices( fh );

    partition_cones();

    // Initial visibility results
    sub_regions.clear();
    for ( int i = 0; i < cones.size(); i++ ) {
      int pts_num = cones[i].end() - cones[i].begin();
      sub_regions.push_back( Sub_region( pts_num ) );
    }

    compute_visibility_parallel( ConcurrencyCategory() );

    merge_result();

    //trace_all( cout );
  }

  template <typename VARR>
  void conditional_regularize( VARR& arr_out, CGAL::Tag_true )
    { regularize_output( arr_out ); }

  template <typename VARR>
  void conditional_regularize( VARR& arr_out, CGAL::Tag_false )
    {}  // do nothing

  template <typename VARR>
  void regularize_output( VARR& arr_out )
  {
    typename VARR::Edge_iterator eit;
    for ( eit = arr_out.edges_begin(); eit != arr_out.edges_end(); eit++ ) {
      if ( eit->face() == eit->twin()->face() ) {
        arr_out.remove_edge( eit );
      }
    }
  }


// private trace mthods
private:
#ifndef NDEBUG
  void trace_all ( ostream& os )
  {
    os << "***********************************" << endl;
    os << "        Trace All" << endl;
    os << "***********************************" << endl;
    os << "query_pt = [" << query_pt << "]" << endl;
    os << "query_type = [";
    switch( query_type ) {
    case VERTEX_QUERY:
      os << "Vertex query";
      break;
    case EDGE_QUERY:
      os << "Edge query";
      break;
    case FACE_QUERY:
      os << "Face query";
      break;
    }
    os << "]" << endl;
    os << "vs = " << endl;
    for ( int i = 0; i < vs.size(); i++ ) {
      os << "  idx = [" << i << "] : ";
      vs[i].trace( os, 0 );
    }
    os << "es = " << endl;
    for ( int i = 0; i < es.size(); i++ )  {
      os << "  idx = [" << i << "] : ";
      es[i].trace( os, 0 );
    }
    os << "good_vdx = " << endl;
    for ( int i = 0; i < good_vdx.size(); i++ )
      os << "  vdx[" << i << "] = [" << good_vdx[i] << "]" << endl;
    os << "incident = " << endl;
    for ( int i = 0; i < incident.size(); i++ ) {
      os << "  incident[" << i << "] : vidx = " << incident[i].first << " , eidx = " << incident[i].second << endl;
    }
    if ( query_type != FACE_QUERY ) {
      os << "query_edge = [" << query_edge->source()->point()
         << " -> " << query_edge->target()->point()
         << ", is_small_cone=[" << is_small_cone
         << "]" << endl;
    }
    os << "cones = " << endl;
    for ( int i = 0; i < cones.size(); i++ ) {
      os << "  " << "cone[" << i << "] =" << endl;
      cones[i].trace( os, 2 );
    }
    os << "sub_regions = " << endl;
    for ( int i = 0; i < sub_regions.size(); i++ ) {
      os << "  " << "sub_regions[" << i << "] =" << endl;
      sub_regions[i].trace( os, 2 );
    }
    os << "polygon = " << endl;
    for ( int i = 0; i < polygon.size(); i++ ) {
      os << "  " << polygon[i] << endl;
    }
    os << "***********************************" << endl;
    os << endl;
  }
#endif

// public methods
public:
  // Constructors
  Parallel_rotational_sweep_visibility_2 ()
    : p_arr(NULL), geom_traits(NULL)
    {}
  Parallel_rotational_sweep_visibility_2 ( const Arrangement_2& arr )
    : p_arr(&arr)
    { geom_traits = p_arr->geometry_traits(); }
  const std::string name ()
    { return std::string("R_visibility_2"); }
  
  // function to compute visibility, query point lies in the interior of a face
  template <typename VARR> 
  typename VARR::Face_handle 
  compute_visibility( const Point_2& q, const Halfedge_const_handle& e,
                      VARR& arr_out )
  {
    if ( q == e->source()->point() )
      return compute_visibility( q, e->prev(), arr_out );
    arr_out.clear();
    query_pt = q;
    query_edge = e;

    if ( query_pt == e->target()->point() ) {
      query_type = VERTEX_QUERY;
      cone_end = e->source();
      cone_start = e->next()->target();
      if ( CGAL::orientation( e->source()->point(),
                              e->target()->point(),
                              e->next()->target()->point() )
           == CGAL::LEFT_TURN )
        is_small_cone = true;
      else
        is_small_cone = false;
    } else {
      query_type = EDGE_QUERY;
      cone_end = e->source();
      cone_start = e->target();
      is_small_cone = false;
    }

    compute_visibility_impl( e->face() );

    // decide which inside of the visibility butterfly
    int small_idx, big_idx;
    if ( cone_end_idx < cone_start_idx ) {
      small_idx = cone_end_idx;
      big_idx = cone_start_idx;
    } else {
      small_idx = cone_start_idx;
      big_idx = cone_end_idx;
    }
    int next_idx = small_idx + 1;
    bool is_between;
    if ( CGAL::right_turn( cone_end->point(), query_pt, cone_start->point() ) ) {
      is_between = false;
      while ( next_idx != big_idx ) {
        if ( CGAL::left_turn( cone_end->point(), query_pt, polygon[next_idx] ) ||
             CGAL::left_turn( query_pt, cone_start->point(), polygon[next_idx] ) ) {
          is_between = true;
          break;
        }
        next_idx++;
      }
    } else {
      is_between = true;
      while ( next_idx != big_idx ) {
        if ( CGAL::right_turn( cone_end->point(), query_pt, polygon[next_idx] ) ||
             CGAL::right_turn( query_pt, cone_start->point(), polygon[next_idx] ) ) {
          is_between = false;
          break;
        }
        next_idx++;
      }
    }

    std::vector<Point_2> polygon_out;
    typename std::vector<Point_2>::iterator first = polygon.begin() + small_idx;
    typename std::vector<Point_2>::iterator last = polygon.begin() + big_idx;
    if ( is_between ) {
      polygon_out.assign( first, last+1 );
      if ( query_type == VERTEX_QUERY )
        polygon_out.push_back( query_pt );
    } else {
      polygon_out.assign( polygon.begin(), first+1 );
      if ( query_type == VERTEX_QUERY )
        polygon_out.push_back( query_pt );
      for ( int i = big_idx; i != polygon.size(); i++ )
        polygon_out.push_back( polygon[i] );
    }

    CGAL::Visibility_2::report_while_handling_needles
          < Parallel_rotational_sweep_visibility_2 >
          ( geom_traits, query_pt, polygon_out, arr_out );

    conditional_regularize( arr_out, Regularization_category() );

    if ( arr_out.faces_begin()->is_unbounded() )
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }

  // function to compute visibility, query point lies on an edge or at a vertex
  template <typename VARR> 
  typename VARR::Face_handle 
  compute_visibility( const Point_2& q, const Face_const_handle f,
                      VARR& arr_out )
  {
    arr_out.clear();
    query_pt = q;
    query_type = FACE_QUERY;
    is_small_cone = false;

    compute_visibility_impl( f );

    CGAL::Visibility_2::report_while_handling_needles
          < Parallel_rotational_sweep_visibility_2 >
          ( geom_traits, query_pt, polygon, arr_out );

    conditional_regularize( arr_out, Regularization_category() );

    if ( arr_out.faces_begin()->is_unbounded() )
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }

  bool is_attached () const
    { return ( p_arr != NULL ); }
  void attach( const Arrangement_2& arr )
    { p_arr = &arr; geom_traits = p_arr->geometry_traits(); }
  void detach ()
    { p_arr = NULL; geom_traits = NULL; }
  const Arrangement_2& arr () const
    { return *p_arr; }

// Private data members
private:
  const Geometry_traits_2 * geom_traits;
  const Arrangement_2 * p_arr;

  Point_2 query_pt;
  Halfedge_const_handle query_edge;
  enum { VERTEX_QUERY, EDGE_QUERY, FACE_QUERY } query_type;
  Vertex_vector vs;
  Edge_vector es;
  std::vector<int> good_vdx;
  std::vector< std::pair<int, int> > incident;
  std::vector<Cone> cones;
  std::vector<Sub_region> sub_regions;
  Point_vector polygon;

  Vertex_const_handle cone_end;
  Vertex_const_handle cone_start;
  int cone_end_idx;
  int cone_start_idx;
  Arrangement_2 arr_box;

  bool is_small_cone;
  Point_2 shifted_source;
  int shifted_quadrant;
};

} // end namespace CGAL



#endif


