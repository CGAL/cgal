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
// $URL$`
// $Id$
//
//
// Author(s):  Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>
//             Ning Xu <longyin0904@gmail.com>

#ifndef CGAL_SIMPLE_POLYGON_VISIBILITY_2_H
#define CGAL_SIMPLE_POLYGON_VISIBILITY_2_H

#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/assertions.h>
#include <stack>
#include <map>

namespace CGAL {

template<class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true> 
class Simple_polygon_visibility_2 {

public:
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            K;

  typedef typename Arrangement_2::Halfedge_const_handle       
                                                        Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
                                        Halfedge_around_vertex_const_circulator;

  typedef typename Geometry_traits_2::Point_2           Point_2;
  typedef typename Geometry_traits_2::Ray_2             Ray_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Line_2            Line_2;
  typedef typename Geometry_traits_2::Vector_2          Vector_2;
  typedef typename Geometry_traits_2::Direction_2       Direction_2;
  typedef typename Geometry_traits_2::FT                Number_type;
  typedef typename Geometry_traits_2::Object_2          Object_2;

  typedef RegularizationCategory              Regularization_category;
  typedef CGAL::Tag_false                     Supports_general_polygon_category;
  typedef CGAL::Tag_true                      Supports_simple_polygon_category;

  Simple_polygon_visibility_2() : p_arr(NULL), traits(NULL) {}

  /*! Constructor given an arrangement and the Regularization tag. */
  Simple_polygon_visibility_2(const Arrangement_2& arr): 
    p_arr(&arr) {
    traits = p_arr->geometry_traits();
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
  }

  
  std::string name() const { return std::string("S_visibility_2"); }

  /*! Method to check if the visibility object is attached or not to
      an arrangement*/
  bool is_attached() const {
    return (p_arr != NULL);
  }

  /*! Attaches the visibility object to the 'arr' arrangement */
  void attach(const Arrangement_2& arr) {
    if(p_arr != &arr){
        detach();
        p_arr = &arr;
        traits = p_arr->geometry_traits();
    }
  }

  /*! Detaches the visibility object from the arrangement it is
      attached to*/
  void detach() {
    p_arr = NULL;
    traits = NULL;
    vertices.clear();
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
    p_cdt.reset();
  }

  /*! Getter method for the input arrangement*/
  const Arrangement_2& arrangement_2() const {
    return *p_arr;
  }

  /*! Computes the visibility object from the query point 'q' in the face
      'face' and constructs the output in 'out_arr'*/
  template <typename VARR> 
  typename VARR::Face_handle 
  compute_visibility(const Point_2& q,
                     const Face_const_handle face,
                     VARR& out_arr) const
  {

    CGAL_precondition(!face->is_unbounded());

    out_arr.clear();
    
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;

    // Now retrieve the circulator to first visible vertex from triangulation
    Ccb_halfedge_const_circulator circ = find_visible_start(face, q);
    Ccb_halfedge_const_circulator curr = circ;

    do {
      vertices.push_back(curr->source()->point());
    } while(++curr != circ);

    vertices.push_back(vertices[0]);

    visibility_region_impl(q);

    return output(q, out_arr);
  }

  /*! Computes the visibility region of the query point 'q' located on the
      halfedge 'he' and constructs the output in 'out_arr'*/
  template <typename VARR> 
  typename VARR::Face_handle 
  compute_visibility(
      const Point_2& q, 
      const Halfedge_const_handle he,
      VARR& out_arr ) const
  {

    out_arr.clear();

    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
    bool query_on_target = false;

    if (q != he->source()->point()) {
      if (q != he->target()->point()) {
        vertices.push_back(he->target()->point());
        query_pt_is_on_halfedge = true;
      }
      else {
        query_pt_is_vertex = true;
        query_on_target = true;
      }
    } else {
      vertices.push_back( he->target()->point() );
      query_pt_is_vertex = true;
    }

    Ccb_halfedge_const_circulator circ = he;
    ++circ;
    Ccb_halfedge_const_circulator curr = circ;

    do {
      const Point_2& curr_vertex = curr->target()->point();
      vertices.push_back(curr_vertex);
    } while (++curr != circ);

    if ( query_on_target ) {
      vertices.push_back( vertices[0] );
    }

    visibility_region_impl(q);

    return output(q, out_arr);

  }

private:
  typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::No_intersection_tag                                Itag;
  typedef CGAL::Constrained_triangulation_2<K, TDS, Itag>          CDT;

private:
  const Arrangement_2 *p_arr;
  const Geometry_traits_2 *traits;

  /*! Boost pointer to the constrained Delaunay triangulation object*/
  mutable boost::shared_ptr<CDT> p_cdt;
  /*! Mapping of the vertices of the input to the corresponding circulator
      needed for finding the first visible vertex in case of face queries*/
  mutable std::map<Point_2, Ccb_halfedge_const_circulator> point_itr_map;
  /*! Stack of visibile points; manipulated when going through the sequence
      of input vertices; contains the vertices of the visibility region after 
      the run of the algorithm*/
  mutable std::stack<Point_2> s;
  /*! Sequence of input vertices*/
  mutable std::vector<Point_2> vertices;
  /*! State of visibility region algorithm*/
  mutable enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} upcase;
  mutable bool query_pt_is_vertex;
  mutable bool query_pt_is_on_halfedge;

  /*! Regularize output if flag is set to true*/
  template <typename VARR> 
  void conditional_regularize(VARR& out_arr, CGAL::Tag_true) const {
    regularize_output(out_arr);
  }

  /*! No need to regularize output if flag is set to false*/
  template <typename VARR> 
  void conditional_regularize(VARR&, CGAL::Tag_false) const {
    //do nothing
  }


  /*! Regularizes the output - removes edges that have the same face on both
      sides */
  template <typename VARR> 
  void regularize_output(VARR& out_arr) const {
    typename VARR::Edge_iterator e_itr;
    for (e_itr = out_arr.edges_begin(); e_itr != out_arr.edges_end(); ++e_itr) {

      if (e_itr->face() == e_itr->twin()->face()) {
        out_arr.remove_edge(e_itr);
      }
    }
  }


  /*! Initialized the constrained Delaunay triangulation using the edges of
      the outer boundary of 'face' */
  void init_cdt(const Face_const_handle &face) const {

    point_itr_map.clear();

    std::vector<std::pair<Point_2,Point_2> > constraints; 
    Ccb_halfedge_const_circulator circ = face->outer_ccb();
    Ccb_halfedge_const_circulator curr = circ;

    do {
      const Point_2& source = curr->source()->point();
      const Point_2& target = curr->target()->point();
      point_itr_map.insert(std::make_pair(source, curr));
      constraints.push_back(std::make_pair(source, target));
    } while(++curr != circ);

    p_cdt = boost::shared_ptr<CDT>(new CDT(constraints.begin(),
                                           constraints.end()));
  }


  template <typename VARR>
  typename VARR::Face_handle
  output(const Point_2& q, VARR& out_arr) const {

      std::vector<Point_2> points;
      while (!s.empty()) {
        points.push_back(s.top());
        s.pop();
      }

      Visibility_2::report_while_handling_needles<Simple_polygon_visibility_2>(
                  traits, q, points, out_arr);

      CGAL_postcondition(out_arr.number_of_isolated_vertices() == 0);
      CGAL_postcondition(s.empty());

      conditional_regularize(out_arr, Regularization_category());
      vertices.clear();

      if (out_arr.faces_begin()->is_unbounded()) {
        return ++out_arr.faces_begin();
      }
      else {
        return out_arr.faces_begin();
      }
  }


  /*! Finds a visible vertex from the query point 'q' in 'face' 
      to start the algorithm from*/
  Ccb_halfedge_const_circulator find_visible_start(Face_const_handle face,
                                                   const Point_2 &q) const {
    init_cdt(face);
    typename CDT::Face_handle fh = p_cdt->locate(q);
    const Point_2& start_point = fh->vertex(0)->point();

    // Now retrieve the circulator to first visible vertex from triangulation
    Ccb_halfedge_const_circulator circ = point_itr_map[start_point];

    Halfedge_around_vertex_const_circulator incident_circ =
            circ->source()->incident_halfedges();
    Halfedge_around_vertex_const_circulator incident_curr = incident_circ;

    Ccb_halfedge_const_circulator incident_next;

    do {

      if (incident_curr->face() == face) {
        incident_next = incident_curr;
        ++incident_next;

        if (Visibility_2::orientation_2(traits,
                                        incident_curr->source()->point(),
                                        incident_curr->target()->point(),
                                        q) == LEFT_TURN
         || Visibility_2::orientation_2(traits,
                                        incident_next->source()->point(),
                                        incident_next->target()->point(),
                                        q) == LEFT_TURN)
        {
          break;
        }
      }
    } while (++incident_curr != incident_circ);

    return incident_next;
  }


  /*! Main method of the algorithm - initializes the stack and variables
      and calles the corresponding methods acc. to the algorithm's state;
      'q' - query point;
      'i' - current vertex' index
      'w' - endpoint of ray shot from query point */
  void visibility_region_impl(const Point_2& q) const {
    int i = 0;
    Point_2 w;
    Orientation orient =
            Visibility_2::orientation_2(traits, q, vertices[0], vertices[1]);

    if ( orient != RIGHT_TURN ) {
      upcase = LEFT;
      i = 1;
      w = vertices[1];
      s.push(vertices[0]);
      s.push(vertices[1]);
    }
    else {
      upcase = SCANA;
      i = 1;
      w = vertices[1];
      s.push(vertices[0]);
    }

    Ray_2 ray_origin( q, vertices[0] );
    do {
      switch(upcase) {
        case LEFT: 
          left(i, w, q);
          break;
        case RIGHT:
          right(i, w, q);
          break;
        case SCANA:
          scana(i, w, q);
          break;
        case SCANB:
          scanb(i, w);
          break;
        case SCANC:
          scanc(i, w);
          break;
        case SCAND:
          scand(i, w);
          break;
        case FINISH:
          break;
      }
      if ( upcase == LEFT ) {
        Point_2 s_t = s.top();
        s.pop();
        if ( Visibility_2::orientation_2(traits, q, vertices[0], s.top() )
             == RIGHT_TURN
             &&
             Visibility_2::orientation_2(traits, q, vertices[0], s_t)
             == LEFT_TURN )
        {
          Segment_2 seg( s.top(), s_t );
          if ( Visibility_2::do_intersect_2(traits, seg, ray_origin) )
          {
            Object_2 result = Visibility_2::intersect_2(traits, seg,ray_origin);
            const Point_2 * ipoint = object_cast<Point_2>(&result);
            assert( ipoint != NULL );
            s_t = *ipoint;
            upcase = SCANB;
          }
        }
        s.push( s_t );
      }
    } while(upcase != FINISH);
  }

  /*! Method that handles the left turns in the vertex algorithm */
  void left(int& i, Point_2& w, const Point_2& q) const {
    if (i >= vertices.size() - 1) {
      upcase = FINISH;
    }
    else {
       Point_2 s_t = s.top();
       s.pop();
       Point_2 s_t_prev = s.top();
       s.push( s_t );
       Orientation orient1 = Visibility_2::orientation_2
                                   ( traits,
                                     q,
                                     vertices[i],
                                     vertices[i+1] );

       if ( orient1 != RIGHT_TURN ) {
         // Case L2
         upcase = LEFT;
         s.push( vertices[i+1] );
         w = vertices[i+1];
         i++;
       } else {
         Orientation orient2 = Visibility_2::orientation_2
                                     ( traits,
                                       s_t_prev,
                                       vertices[i],
                                       vertices[i+1] );

         if ( orient2 == RIGHT_TURN ) {
           // Case L3
           upcase = SCANA;
           w = vertices[i+1];
           i++;
         } else {
           // Case L4
           upcase = RIGHT;
           w = vertices[i];
           i++;
         }
       } 
    }
  }

  /*! Scans the stack such that all vertices that were pushed before to the 
      stack and are now not visible anymore. */
  void right(int& i, Point_2& w, const Point_2& q) const {
     Point_2 s_j;
     Point_2 s_j_prev;
     Point_2 u;
     int mode = 0;
     Orientation orient1, orient2;

     s_j_prev = s.top();
     orient2 = Visibility_2::orientation_2( traits, q, s_j_prev, vertices[i] );

     while ( s.size() > 1 ) {
       s_j = s_j_prev;
       orient1 = orient2;
       s.pop();
       s_j_prev = s.top();

       orient2 = Visibility_2::orientation_2( traits, q, s_j_prev, vertices[i]);
       if ( orient1 != LEFT_TURN && orient2 != RIGHT_TURN ) {
         mode = 1;
         break;
       }

       Segment_2 seg2( vertices[i-1], vertices[i] );
       Segment_2 seg( s_j_prev, s_j );
       if ( vertices[i-1] != s_j &&
            Visibility_2::do_intersect_2(traits, seg, seg2) )
       {
         Object_2 result = Visibility_2::intersect_2(traits, seg, seg2);
         const Point_2 * ipoint = object_cast<Point_2>(&result);
         assert( ipoint != NULL );
         u = *ipoint;
         mode = 2;
         break;
       }
     }

     assert( mode != 0 );
     if ( mode == 1 ) {
       orient1 = Visibility_2::orientation_2
                 ( traits, q, vertices[i], vertices[i+1] );
       orient2 = Visibility_2::orientation_2
                 ( traits, vertices[i-1], vertices[i], vertices[i+1] );

       if ( orient1 == RIGHT_TURN ) {
         // Case R1
         // Since the next action is RIGHT, we do not compute the intersection
         // of (s_j,s_j_prev) and the ray (query_pt, vertices[i]),
         // thus, (s_j,s_j_prev) is not shortcutted, but it is harmless
         upcase = RIGHT;
         s.push( s_j );
         w = vertices[i];
         i++;
       } else if ( orient2 == RIGHT_TURN ) {
         // Case R2
         Ray_2 ray( q, vertices[i] );
         Segment_2 seg( s_j_prev, s_j );

         Object_2 result = Visibility_2::intersect_2( traits, seg, ray );
         const Point_2 * ipoint = object_cast<Point_2>(&result);

         assert( ipoint != NULL );

         u = *ipoint;
         if ( s.top() != u ) {
           s.push( u );
         }
         upcase = LEFT;
         s.push( vertices[i] );
         s.push( vertices[i+1] );
         w = vertices[i+1];
         i++;
       } else {
         // Case R3
         Ray_2 ray( q, vertices[i] );
         Segment_2 seg( s_j_prev, s_j );

         Object_2 result = Visibility_2::intersect_2( traits, seg, ray );
         const Point_2 * ipoint = object_cast<Point_2>(&result);

         assert( ipoint != NULL );

         u = *ipoint;
         if ( s.top() != u ) {
           s.push( u );
         }
         upcase = SCANC;
         w = vertices[i];
         i++;
       }
     } else if ( mode == 2 ) {
       // Case R4
       upcase = SCAND;
       w = u;
     }
  }

  /*! Scans the vertices starting from index 'i' for the first visible vertex
      out of the back hidden window */
  void scana(int& i, Point_2& w, const Point_2& q) const {
    // Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
    Point_2 u;
    int k = scan_edges( i, q, s.top(), u, true );

    Orientation orient1 =
            Visibility_2::orientation_2(traits, q, vertices[k], vertices[k+1] );

    if ( orient1 == RIGHT_TURN ) {
      bool fwd = Visibility_2::collinear_are_ordered_along_line_2
                 ( traits, q, s.top(), u );

      if ( !fwd ) {
        // Case A1
        upcase = RIGHT;
        i = k+1;
        w = u;
      } else {
        // Case A2
        upcase = SCAND;
        i = k+1;
        w = u;
      }
    } else {
      // Case A3
      upcase = LEFT;
      i = k+1;
      s.push( u );
      if ( u != vertices[k+1] ) {
        s.push( vertices[k+1] );
      }
      w = vertices[k+1];
    }
  }

  /*! Find the first edge interecting the segment (v_0, s_t) */
  void scanb(int& i, Point_2& w) const {
    if ( i == vertices.size() - 1 ) {
      upcase = FINISH;
      return;
    }
    Point_2 u;
    int k = scan_edges( i, s.top(), vertices[0], u, false );
    if ( (k+1 == vertices.size()-1) && (vertices[0] == u) ) {
      // Case B1
      upcase = FINISH;
      s.push( vertices[0] );
    } else {
      // Case B2
      upcase = RIGHT;
      i = k+1;
      w = u;
    }
  }

  /*! Finds the exit from a general front hidden window by finding the first
      vertex to the right of the ray defined by the query_point and w*/
  void scanc(int& i, Point_2& w) const {
    Point_2 u;
    int k = scan_edges( i, s.top(), w, u, false );
    upcase = RIGHT;
    i = k+1;
    w = u;
  }

  /*! find the first edge intersecting the given window (s_t, w) */
  void scand(int& i, Point_2& w) const {
    Point_2 u;
    int k = scan_edges( i, s.top(), w, u, false );
    upcase = LEFT;
    i = k+1;
    s.push( u );
    if ( u != vertices[k+1] ) {
      s.push( vertices[k+1] );
    }
    w = vertices[k+1];
  }
  
  /*! Scan edges v_i,v_{i+1},...,v_n, until find an edge intersecting given ray
      or given segment. is_ray = true -> ray, false -> segment.
      The intersection point is returned by u */
  int scan_edges( int i,
                  const Point_2& ray_begin,
                  const Point_2& ray_end,
                  Point_2& u,
                  bool is_ray ) const
  {
    Orientation old_orient = RIGHT_TURN;
    Ray_2 ray( ray_begin, ray_end );
    Segment_2 s2( ray_begin, ray_end );
    int k;
    Object_2 result;
    for ( k = i; k+1 < vertices.size(); k++ ) {
      Orientation curr_orient = Visibility_2::orientation_2
                                      ( traits,
                                        ray_begin,
                                        ray_end,
                                        vertices[k+1] );
      if ( curr_orient != old_orient ) {
        // Orientation switch, an intersection may occur
        Segment_2 seg( vertices[k], vertices[k+1] );
        if ( is_ray ) {
          if (Visibility_2::do_intersect_2(traits, seg, ray) )
          {
            result = Visibility_2::intersect_2( traits, seg, ray );
            break;
          }
        } else {
          if (Visibility_2::do_intersect_2(traits, seg, s2) )
          {
            result = Visibility_2::intersect_2( traits, seg, s2 );
            break;
          }
        }
      }
      old_orient = curr_orient;
	}
    assert( k+1<vertices.size() );
    const Point_2 * ipoint = object_cast<Point_2>( &result );
    if ( ipoint ) {
       u = *ipoint;
    } else {
       u = vertices[k+1];
    }
	return k;
  }
};

} // namespace CGAL
#endif
