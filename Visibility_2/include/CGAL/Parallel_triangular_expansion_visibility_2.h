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
// Author(s):  Michael Hemmer <michael.hemmer@cgal.org>
//             Ning Xu <longyin0904@gmail.com>
//             

#ifndef CGAL_PARALLEL_TRIANGULAR_EXPANSION_VISIBILITY_2__H
#define CGAL_PARALLEL_TRIANGULAR_EXPANSION_VISIBILITY_2__H

#include <CGAL/Arrangement_2.h>
#include <boost/shared_ptr.hpp>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <queue>

#ifdef CGAL_LINKED_WITH_TBB
  #include "tbb/parallel_do.h"
  #include "tbb/concurrent_vector.h"
#endif

namespace CGAL {

template < class Arrangement_2_,
           class RegularizationCategory = CGAL::Tag_true,
           class ConcurrencyCategory = CGAL::Parallel_tag >
class Parallel_triangular_expansion_visibility_2
{
  typedef typename Arrangement_2_::Geometry_traits_2    Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            K;
public:
  // Currently only consider with same type for both
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                        Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;

  typedef typename K::Point_2                           Point_2;
  typedef typename Geometry_traits_2::Ray_2             Ray_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Line_2            Line_2;
  typedef typename Geometry_traits_2::Vector_2          Vector_2;
  typedef typename Geometry_traits_2::Direction_2       Direction_2;
  typedef typename Geometry_traits_2::FT                Number_type;
  typedef typename Geometry_traits_2::Object_2          Object_2;

  // TODO 
  typedef RegularizationCategory                             Regularization_category;
  
  typedef CGAL::Tag_true                                Supports_general_polygon_category;
  typedef CGAL::Tag_true                                Supports_simple_polygon_category;    

private:
  typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::No_intersection_tag                                Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;

  typedef typename CDT::Face_handle                     CDT_Face_handle;
  typedef std::pair<CDT_Face_handle,int>                CDT_Edge;
  typedef typename CDT::Vertex_handle                   CDT_Vertex_handle;

private:
  class Edge
  {
  private:
    CDT_Face_handle fh;
    int idx;
    CDT_Vertex_handle lvh;
    CDT_Vertex_handle rvh;
  public:
    Edge ()
      : fh( NULL ), idx( -1 ), lvh( NULL ), rvh( NULL )
      {}
    Edge ( CDT_Face_handle f, int i, CDT_Vertex_handle l, CDT_Vertex_handle r )
      : fh( f ), idx( i ), lvh( l ), rvh( r )
      {}
    CDT_Face_handle face () const
      { return fh; }
    int index () const
      { return idx; }
    CDT_Vertex_handle left_source () const
      { return lvh; }
    CDT_Vertex_handle right_source () const
      { return rvh; }
    CDT_Vertex_handle left_vertex () const
      { return fh->vertex( fh->cw( idx ) ); }
    CDT_Vertex_handle right_vertex () const
      { return fh->vertex( fh->ccw( idx ) ); }
    bool left_ray_degenerate () const
      { return ( lvh->point() == left_vertex()->point() ); }
    bool right_ray_degenerate () const
      { return ( rvh->point() == right_vertex()->point() ); }

    CDT_Edge edge () const
      { return CDT_Edge( fh, idx ); }

    void trace ( ostream& os, int level ) const
    {
      for ( int i = 0; i < level; i++ ) os << "  ";
      os << "lvh = " << lvh->point() << ", rvh = " << rvh->point() << " , left = " << left_vertex()->point() << " , right = " << right_vertex()->point() << endl;
    } 
  };

  class Edge_bundle
  {
  private:
    std::queue<Edge> data;
    std::vector<Edge> constrained;
  public:
    Edge_bundle ()
      {}

    int capacity () const
      { return 128; }
    int size () const
      { return data.size(); }
    bool empty () const
      { return data.empty(); }
    bool full () const
      { return ( data.size() == capacity() ); }
    int constrained_size () const
      { return constrained.size(); }
    const Edge& get_constrained( int pos ) const
      { return constrained[pos]; }

    const Edge& front () const
      { return data.front(); }
    void push( const Edge& e )
      { data.push( e ); }
    void pop ()
      { data.pop(); }
    void clear ()
    {
      while ( !data.empty() ) data.pop();
      constrained.clear();
    }
    void push_constrained ( const Edge& e )
    {
      constrained.push_back( e );
    }

    void trace ( ostream& os, int level ) const
    {
      std::queue<Edge> tmp = data;
      int i = 0;
      while ( !tmp.empty() ) {
        for ( int j = 0; j < level; j++ ) os << "  ";
        os << "edge[" << i << "] = ";
        tmp.front().trace( os, level+1 );
        i++;
        tmp.pop();
      }
    }
  };

private:
  typedef std::vector<Edge_bundle>                      EB_vector;
#ifdef CGAL_LINKED_WITH_TBB
  typedef tbb::concurrent_vector<Edge>                  Result_container;
#else
  typedef std::vector<Edge>                             Result_container;
#endif

private:

#ifdef CGAL_LINKED_WITH_TBB
  class Parallel_expand 
  {
  private:
    Parallel_triangular_expansion_visibility_2 * algo;
    Result_container * p_result;
  public:
    Parallel_expand( Parallel_triangular_expansion_visibility_2 * a,
                     Result_container * pr )
      : algo( a ), p_result( pr )
      {}
    void operator () ( Edge_bundle& bundle, tbb::parallel_do_feeder<Edge_bundle>& feeder ) const
    {
      algo->compute_visibility_bundle( bundle );
      for ( int i = 0; i < bundle.constrained_size(); i++ )
        p_result->push_back( bundle.get_constrained( i ) );
      if ( bundle.empty() )
        return;
      Edge_bundle res;
      while ( !bundle.empty() ) {
        res.push( bundle.front() );
        bundle.pop();
        if ( res.full() ) {
          feeder.add( res );
          res.clear();
        }
      }
      if ( !res.empty() )
        feeder.add( res );
    }
  };
#endif

private:
  Edge left_edge ( CDT_Face_handle f, int idx )
    { return Edge( f, f->cw(idx), f->vertex(f->ccw(idx)), f->vertex(idx) ); }
  Edge right_edge ( CDT_Face_handle f, int idx )
    { return Edge( f, f->ccw(idx), f->vertex(idx), f->vertex(f->cw(idx)) ); }

  CDT_Face_handle find_face_having_halfedge( CDT_Face_handle fh, int index )
  {
    CDT_Face_handle curr, circ, nfh;
    int curr_idx, nindex;

    curr = circ = fh;
    curr_idx = index;
    do {
      if ( curr->vertex( curr->cw( curr_idx ) )->point() == query_edge->source()->point() )
        return curr;
      nfh = curr->neighbor( curr->cw(curr_idx) );
      nindex = nfh->index( curr );
      curr = nfh;
      curr_idx = nfh->cw( nindex );
    } while ( curr != circ && !p_cdt->is_infinite( curr ) );

    if ( p_cdt->is_infinite( curr ) ) {
      curr = circ = fh;
      curr_idx = index;
      do {
        if ( curr->vertex( curr->cw( curr_idx ) )->point() == query_edge->source()->point() )
          return curr;
        nfh = curr->neighbor( curr->ccw(curr_idx) );
        nindex = nfh->index( curr );
        curr = nfh;
        curr_idx = nfh->ccw( nindex );
      } while ( curr != circ && !p_cdt->is_infinite( curr ) );
    }
    // Something must be wrong here!
    return NULL;
  }
  int find_index_of_halfedge( CDT_Face_handle fh )
  {
    for ( int i = 0; i < 3; i++ ) {
      if ( fh->vertex( i )->point() == query_pt )
        return i;
    }
    // Something must be wrong here!
    return -1;
  }
  void insert_edges_of_vertex_query( Edge_bundle& edges, CDT_Face_handle fh, int index )
  {
    CDT_Face_handle curr = fh;
    int curr_idx = index;
    do {
      edges.push( left_edge( curr, curr->ccw( curr_idx ) ) );
      CDT_Face_handle nfh = curr->neighbor( curr->cw(curr_idx) );
      int nindex = nfh->index( curr );
      curr = nfh;
      curr_idx = nfh->cw( nindex );
    } while ( curr->vertex(curr->cw(curr_idx))->point() != query_edge->next()->target()->point() );
  }

  int expand_threshold ()
    { return 1024; }

  void compute_visibility_bundle( Edge_bundle& bundle )
  {
    int count = 0;
    int threshold = expand_threshold();
    while ( !bundle.empty() && count < threshold ) {
      const Edge& e = bundle.front();

      if ( p_cdt->is_constrained( e.edge() ) ) {
        // A constrained edge
        bundle.push_constrained( e );
      } else {
        // An unconstrained edge, should be expanded
        // deal with the left edge

        CDT_Face_handle nfh = e.face()->neighbor( e.index() );
        int nindex = nfh->index( e.face() );
        int rindex = p_cdt->ccw( nindex );
        int lindex = p_cdt->cw( nindex );

        CDT_Vertex_handle nvh = nfh->vertex( nindex );
        CDT_Vertex_handle rvh = nfh->vertex( lindex );
        CDT_Vertex_handle lvh = nfh->vertex( rindex );

        CGAL::Orientation ro = CGAL::orientation( query_pt,
                                                  e.right_source()->point(),
                                                  nvh->point() );
        CGAL::Orientation lo = CGAL::orientation( query_pt,
                                                  e.left_source()->point(),
                                                  nvh->point() );

        if ( ro == CGAL::LEFT_TURN ) {
          // The right edge is seen
          if ( lo == CGAL::LEFT_TURN ) {
            // No split is needed
            bundle.push( Edge( nfh, rindex, e.left_source(), e.right_source() ) );
            count++;
          } else if ( lo == CGAL::RIGHT_TURN ) {
            // split at the new vertex
            bundle.push( Edge( nfh, rindex, nvh, e.right_source() ) );
            count++;
          } else {
            // the right edge is not seen, but the opposite vertex is seen
            bundle.push( Edge( nfh, rindex, nvh, e.right_source() ) );
            count++;
            bundle.push_constrained( Edge( nfh, rindex, e.left_source(), e.left_source() ) );
          }
        }

        if ( lo == CGAL::RIGHT_TURN ) {
          // The left edge is seen
          if ( ro == CGAL::RIGHT_TURN ) {
            // No split is needed
            bundle.push( Edge( nfh, lindex, e.left_source(), e.right_source() ) );
            count++;
          } else if ( ro == CGAL::LEFT_TURN ) {
            // split at the new vertex
            bundle.push( Edge( nfh, lindex, e.left_source(), nvh ) );
            count++;
          } else {
            // the left edge is not seen, but the opposite vertex is seen
            bundle.push_constrained( Edge( nfh, lindex, e.right_source(), e.right_source() ) );
          }
        }
      }

      bundle.pop();
    }
  }

  void compute_visibility_parallel ( Edge_bundle& bundle,
                                     Result_container& result,
                                     CGAL::Sequential_tag )
  {
  }

  void compute_visibility_parallel ( Edge_bundle& bundle,
                                     Result_container& result,
                                     CGAL::Parallel_tag )
  {
#ifdef CGAL_LINKED_WITH_TBB
    std::vector<Edge_bundle> bundles;
    bundles.push_back( bundle );
    Parallel_expand exp( this, &result );
    tbb::parallel_do( bundles.begin(), bundles.end(), exp );
#else
    compute_visibility_parallel( bundle, CGAL::Sequential_tag() );
#endif
  }


  Point_2 ray_seg_intersection ( const Point_2& dp,
                                 const Point_2& s,
                                 const Point_2& t )
  {
    Ray_2 ray( query_pt, dp );
    Segment_2 seg( s, t );
    assert( CGAL::do_intersect( ray, seg ) );
    CGAL::Object obj = CGAL::intersection( ray, seg );
    const Point_2 * ipoint = CGAL::object_cast<Point_2>( &obj );
    assert( ipoint != NULL );
    return *ipoint;
  }

  template < class VARR >
  void compute_visibility_impl ( VARR& arr_out)
  {
    arr_out.clear();
    // Locate the query point
    typename CDT::Locate_type location;
    int index;
    CDT_Face_handle fh = p_cdt->locate( query_pt, location, index );

    Edge_bundle edges;

    switch( location ) {
    case CDT::FACE:
      // The query point locates in a face
      // It only occurs in FACE_QUERY
      assert( query_type == FACE_QUERY );
      edges.push( Edge( fh, 0, fh->vertex(2), fh->vertex(1) ) );
      edges.push( Edge( fh, 1, fh->vertex(0), fh->vertex(2) ) );
      edges.push( Edge( fh, 2, fh->vertex(1), fh->vertex(0) ) );
      break;
    case CDT::VERTEX:
      // The query point locates at a vertex
      // It only occurs in VERTEX_QUERY
      // Find the face containing the halfedge in query
      assert( query_type == VERTEX_QUERY );
      fh = find_face_having_halfedge( fh, index );
      assert( fh != NULL );
      index = find_index_of_halfedge( fh );
      assert( index >= 0 );

      insert_edges_of_vertex_query( edges, fh, index );
      break;
    case CDT::EDGE:
      // The query point locates on an edge
      // It may occur in FACE_QUERY and EDGE_QUERY
      assert( query_type != VERTEX_QUERY );
      if ( !p_cdt->is_constrained( std::make_pair(fh, index) ) ) {
        // the edge is not contrained
        CDT_Face_handle nfh = fh->neighbor( index );
        int nindex = nfh->index( fh );
        assert( !p_cdt->is_infinite( fh ) );
        assert( !p_cdt->is_infinite( nfh ) );
        edges.push( left_edge( fh, index ) );
        edges.push( right_edge( fh, index ) );
        edges.push( left_edge( nfh, nindex  ) );
        edges.push( right_edge( nfh, nindex  ) );
      } else {
        // the edge is contrained
        if ( fh->vertex(fh->cw(index))->point() == query_edge->target()->point() ) {
          // the locating edge is the halfedge in query
          assert( !p_cdt->is_infinite( fh ) );
          edges.push( left_edge( fh, index ) );
          edges.push( right_edge( fh, index ) );
        } else {
          // the locating edge is not the halfedge in query
          CDT_Face_handle nfh = fh->neighbor( index );
          int nindex = nfh->index( fh );
          assert( !p_cdt->is_infinite( nfh ) );
          edges.push( left_edge( nfh, nindex ) );
          edges.push( right_edge( nfh, nindex ) );
        }
      }
      break;
    }

    Result_container results;

    compute_visibility_parallel( edges, results, ConcurrencyCategory() );

    // Generate results
    std::vector<Segment_2> segs;
    for ( int i = 0; i < results.size(); i++ ) {
      // push segments on the edge
      if ( results[i].left_source()->point() == results[i].right_source()->point() ) {
        // Special case
        // the constrained edge is introduced by a single ray
        if ( CGAL::orientation( query_pt, results[i].left_source()->point(),
                                results[i].left_vertex()->point() )
                                == CGAL::COLLINEAR ) { 
          // the left vertex is shoot by the ray
          if ( !results[i].left_ray_degenerate() ) {
            segs.push_back( Segment_2( results[i].left_vertex()->point(),
                                       results[i].left_source()->point() ) );
          }
        }  else {
          // the right vertex is shoot by the ray
          if ( !results[i].right_ray_degenerate() ) {
            segs.push_back( Segment_2( results[i].right_source()->point(),
                                       results[i].right_vertex()->point() ) );
          }
        }
        continue;
      }

      Point_2 lp = ray_seg_intersection( results[i].left_source()->point(),
                                         results[i].left_vertex()->point(),
                                         results[i].right_vertex()->point() );
      Point_2 rp = ray_seg_intersection( results[i].right_source()->point(),
                                         results[i].left_vertex()->point(),
                                         results[i].right_vertex()->point() );
      if ( lp != rp )
        segs.push_back( Segment_2( rp, lp ) );

      //push segments on the left side
      if ( results[i].left_source() != results[i].left_vertex() ) {
        segs.push_back( Segment_2( lp, results[i].left_source()->point() ) );
      }

      //push segments on the right side
      if ( results[i].right_source() != results[i].right_vertex() ) {
        segs.push_back( Segment_2( results[i].right_source()->point(), rp ) );
      }
    }

    CGAL::insert( arr_out, segs.begin(), segs.end() );
  }

  void init_cdt ()
  {
    std::vector< std::pair<Point_2, Point_2> > constraints;
    typename Arrangement_2::Edge_const_iterator eit;
    for ( eit = p_arr->edges_begin(); eit != p_arr->edges_end(); eit++ ) {
      constraints.push_back( std::make_pair( eit->source()->point(),
                                             eit->target()->point() ) );
    }
    CDT * p = new CDT( constraints.begin(), constraints.end() );
    assert( p != NULL );
    p_cdt = boost::shared_ptr<CDT>( p );
  }
  
public: 
  Parallel_triangular_expansion_visibility_2 ()
    : p_arr( NULL )
    {}
  Parallel_triangular_expansion_visibility_2 ( const Arrangement_2& arr )
    : p_arr( &arr )
    { init_cdt(); }
  void attach ( const Arrangement_2& arr )
  {
    if ( p_arr != &arr ) {
      p_arr = arr;
      init_cdt();
    }
  }
  void detach ()
  {
    p_arr = NULL;
    p_cdt = boost::shared_ptr<CDT>();
  }
  bool is_attached () const
    { return ( p_arr != NULL ); }
  const Arrangement_2& arr () const
    { return &p_arr; }
  const std::string name () const
    { return std::string( "PT_visibility_2" ); }


  template < class VARR >
  typename VARR::Face_handle
  compute_visibility ( const Point_2& q,
                       const Face_const_handle& fh,
                       VARR& arr_out )
  {
    query_pt = q;
    query_type = FACE_QUERY;


    compute_visibility_impl( arr_out );

    assert( arr_out.number_of_faces() == 2 );

    if ( arr_out.faces_begin()->is_unbounded() )
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }

  template < class VARR >
  typename VARR::Face_handle
  compute_visibility ( const Point_2& q,
                       const Halfedge_const_handle& eh,
                       VARR& arr_out )
  {
    if ( q == eh->source()->point() )
      return compute_visibility( q, eh->prev(), arr_out );

    query_pt = q;
    query_edge = eh;
    if ( query_pt == eh->target()->point() ) {
      query_type = VERTEX_QUERY;
    } else {
      query_type = EDGE_QUERY;
    }
    compute_visibility_impl( arr_out );

    assert( arr_out.number_of_faces() == 2 );

    if ( arr_out.faces_begin()->is_unbounded() )
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }


private:
  const Arrangement_2* p_arr;
  boost::shared_ptr<CDT> p_cdt; 

  Point_2 query_pt;
  enum { VERTEX_QUERY, EDGE_QUERY, FACE_QUERY } query_type;
  Halfedge_const_handle query_edge;
};

} // namespace CGAL

#endif //CGAL_TRIANGULAR_EXPANSION_VISIBILITY_2__H
