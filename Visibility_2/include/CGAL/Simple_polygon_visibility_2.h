// Copyright (c) 2013 Technical University Braunschweig (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):  Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>
//             Ning Xu <longyin0904@gmail.com>

#ifndef CGAL_SIMPLE_POLYGON_VISIBILITY_2_H
#define CGAL_SIMPLE_POLYGON_VISIBILITY_2_H

#include <CGAL/license/Visibility_2.h>


#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/assertions.h>
#include <stack>

// TODO:
// * fix handle needles = O(nlogn)

namespace CGAL {

  template<class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true>
    class Simple_polygon_visibility_2 {

  public:
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            K;

  typedef typename K::Intersect_2                       Intersect_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
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
  typedef typename Geometry_traits_2::Object_2          Object_2;

  typedef RegularizationCategory              Regularization_category;
  typedef CGAL::Tag_false                     Supports_general_polygon_category;
  typedef CGAL::Tag_true                      Supports_simple_polygon_category;

  Simple_polygon_visibility_2() : p_arr(nullptr), traits(nullptr) {}

  /*! Constructor given an arrangement and the Regularization tag. */
  Simple_polygon_visibility_2(const Arrangement_2& arr):
  p_arr(&arr) {
    traits = p_arr->geometry_traits();
    point_location.attach(arr);
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
    inserted_artificial_starting_vertex = false;
  }


  std::string name() const { return std::string("S_visibility_2"); }

  /*! Method to check if the visibility object is attached or not to
    an arrangement*/
  bool is_attached() const {
    return (p_arr != nullptr);
  }

  /*! Attaches the visibility object to the 'arr' arrangement */
  void attach(const Arrangement_2& arr) {
    if(p_arr != &arr){
      detach();
      p_arr = &arr;
      traits = p_arr->geometry_traits();
      point_location.attach(arr);
    }
  }

  /*! Detaches the visibility object from the arrangement it is
    attached to*/
  void detach() {
    point_location.detach();
    p_arr = nullptr;
    traits = nullptr;
    vertices.clear();
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
    inserted_artificial_starting_vertex = false;
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
    inserted_artificial_starting_vertex = false;

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
  typedef Arr_walk_along_line_point_location<Arrangement_2>  Arr_point_location;
  typedef typename Arr_point_location::result_type           Location_result;

  typedef std::vector<Point_2>                                 Vertex_container;
  typedef typename Vertex_container::size_type                 Size_type;

  const Arrangement_2 *p_arr;
  const Geometry_traits_2 *traits;

  mutable Arr_point_location point_location;

  /*! Stack of visibile points; manipulated when going through the sequence
    of input vertices; contains the vertices of the visibility region after
    the run of the algorithm*/
  mutable std::stack<Point_2> stack;
  /*! Sequence of input vertices*/
  mutable Vertex_container vertices;
  /*! State of visibility region algorithm*/
  mutable enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} upcase;
  mutable bool query_pt_is_vertex;
  mutable bool query_pt_is_on_halfedge;
  mutable bool inserted_artificial_starting_vertex;


  template <typename VARR>
  typename VARR::Face_handle
  output(const Point_2& q, VARR& out_arr) const {

    if(inserted_artificial_starting_vertex)
      stack.pop();

    std::vector<Point_2> points;
    while(!stack.empty()) {
      const Point_2& top = stack.top();
      if (top != q || query_pt_is_vertex) {
        points.push_back(top);
      }
      stack.pop();
    }

    if(inserted_artificial_starting_vertex) {
      points.back() = points[0];
      inserted_artificial_starting_vertex = false;
    }


    // Quick fix for now. Can be done faster
    bool is_degenerate = false;

    for(typename std::vector<Point_2>::size_type i = 0; i < points.size()-2;i++){
      if(CGAL::orientation(points[i],points[i+1],points[i+2]) == CGAL::COLLINEAR){
        is_degenerate = true;
        break;
      }
    }
    if(is_degenerate){
      //std::cout << is_degenerate << std::endl;
      std::vector<Segment_2> segments;

      for(typename std::vector<Point_2>::size_type i = 0;i < points.size() - 1; ++i)
        {
          segments.push_back(Segment_2(points[i], points[i+1]));
        }
      CGAL::insert(out_arr, segments.begin(), segments.end());
    }else{
      points.pop_back();
      //std::cout << " ordanary " << std::endl;
      typename VARR::Vertex_handle v_last, v_first;
      v_last = v_first =
        out_arr.insert_in_face_interior(points[0],out_arr.unbounded_face());

      for(unsigned int i = 0; i < points.size()-1; i++){
        if(points[i] < points[(i+1)]){
          v_last = out_arr.insert_from_left_vertex (
                                                    Segment_2(points[i], points[i+1]), v_last
                                                    )->target();
        } else {
          v_last = out_arr.insert_from_right_vertex(
                                                    Segment_2(points[i], points[i+1]), v_last
                                                    )->target();
        }
      }
      out_arr.insert_at_vertices(
                                 Segment_2(points.front(), points.back()),
                                 v_last, v_first
                                 );

    }

    CGAL_postcondition(out_arr.number_of_isolated_vertices() == 0);
    CGAL_postcondition(stack.empty());

    Visibility_2::conditional_regularize(out_arr, Regularization_category());
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
                                                   const Point_2 &q) const
  {
    Location_result result = point_location.ray_shoot_up(q);

    if(const Halfedge_const_handle* e =
       boost::get<Halfedge_const_handle>(&(result)))
      {
        CGAL_assertion((*e)->face() == face);
        Point_2 p(q.x(),
                  traits->compute_y_at_x_2_object()(
                                                    Line_2((*e)->source()->point(),
                                                           (*e)->target()->point()) ,
                                                    q.x()));

        vertices.push_back(p);
        inserted_artificial_starting_vertex = true;

        return (*e)->next()->ccb();
      }
    else if (const Vertex_const_handle* v =
             boost::get<Vertex_const_handle>(&(result)))
      {
        Halfedge_around_vertex_const_circulator cir =
        (*v)->incident_halfedges();

        while(face != cir->face()) {
          ++cir;
        }
        return cir->next()->ccb();
      }
    else
      {
        CGAL_assertion_msg(false, "Should not be reachable.");
        return Ccb_halfedge_const_circulator();
      }
  }


  /*! Main method of the algorithm - initializes the stack and variables
    and calles the corresponding methods acc. to the algorithm's state;
    'q' - query point;
    'i' - current vertex' index
    'w' - endpoint of ray shot from query point */
  void visibility_region_impl(const Point_2& q) const {
    Size_type i = 0;
    Point_2 w;
    Orientation o = traits->orientation_2_object()(q, vertices[0], vertices[1]);

    if ( o != RIGHT_TURN ) {
      upcase = LEFT;
      i = 1;
      w = vertices[1];
      stack.push(vertices[0]);
      stack.push(vertices[1]);
    }
    else {
      upcase = SCANA;
      i = 1;
      w = vertices[1];
      stack.push(vertices[0]);
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
        Point_2 s_t = stack.top();
        stack.pop();
        if (traits->orientation_2_object()(q, vertices[0], stack.top() )
            == RIGHT_TURN
            &&
            traits->orientation_2_object()(q, vertices[0], s_t)
            == LEFT_TURN )
          {
            Segment_2 seg( stack.top(), s_t );
            if (Object_2 result = Intersect_2()(seg, ray_origin) )
              {
                const Point_2 * ipoint = object_cast<Point_2>(&result);
                CGAL_assertion( ipoint != nullptr );
                s_t = *ipoint;
                upcase = SCANB;
              }
          }
        stack.push( s_t );
      }
    } while(upcase != FINISH);
  }

  /*! Method that handles the left turns in the vertex algorithm */
  void left(Size_type& i, Point_2& w, const Point_2& q) const {
    if (i >= vertices.size() - 1) {
      upcase = FINISH;
    }
    else {
      Point_2 s_t = stack.top();
      stack.pop();
      Point_2 s_t_prev = stack.top();
      stack.push( s_t );
      Orientation orient1 = traits->orientation_2_object()(
                                                           q,
                                                           vertices[i],
                                                           vertices[i+1] );

      if ( orient1 != RIGHT_TURN ) {
        // Case L2
        upcase = LEFT;
        stack.push( vertices[i+1] );
        w = vertices[i+1];
        i++;
      } else {
        Orientation orient2 = traits->orientation_2_object()(
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
  void right(Size_type& i, Point_2& w, const Point_2& q) const {
    Point_2 s_j;
    Point_2 s_j_prev;
    Point_2 u;
    int mode = 0;
    Orientation orient1, orient2;

    s_j_prev = stack.top();
    orient2 = traits->orientation_2_object()( q, s_j_prev, vertices[i] );

    while ( stack.size() > 1 ) {
      s_j = s_j_prev;
      orient1 = orient2;
      stack.pop();
      s_j_prev = stack.top();

      orient2 = traits->orientation_2_object()( q, s_j_prev, vertices[i]);
      if ( orient1 != LEFT_TURN && orient2 != RIGHT_TURN ) {
        mode = 1;
        break;
      }

      Segment_2 seg2( vertices[i-1], vertices[i] );
      Segment_2 seg( s_j_prev, s_j );
      if ( vertices[i-1] != s_j )
        {
          Object_2 result = Intersect_2()( seg, seg2 );
          if(result) {
            const Point_2 * ipoint = object_cast<Point_2>(&result);
            CGAL_assertion( ipoint != nullptr );
            u = *ipoint;
            mode = 2;
            break;
          }
        }
    }

    CGAL_assertion( mode != 0 );
    if ( mode == 1 ) {
      orient1 = traits->orientation_2_object()(q, vertices[i], vertices[i+1] );

      orient2 = traits->orientation_2_object()(vertices[i-1],
                                               vertices[i],
                                               vertices[i+1] );

      if ( orient1 == RIGHT_TURN ) {
        // Case R1
        // Since the next action is RIGHT, we do not compute the intersection
        // of (s_j,s_j_prev) and the ray (query_pt, vertices[i]),
        // thus, (s_j,s_j_prev) is not shortcutted, but it is harmless
        upcase = RIGHT;
        stack.push( s_j );
        w = vertices[i];
        i++;
      } else if ( orient2 == RIGHT_TURN ) {
        // Case R2
        Ray_2 ray( q, vertices[i] );
        Segment_2 seg( s_j_prev, s_j );

        Object_2 result = Intersect_2()( seg, ray );
        const Point_2 * ipoint = object_cast<Point_2>(&result);

        CGAL_assertion( ipoint != nullptr );

        u = *ipoint;
        if ( stack.top() != u ) {
          stack.push( u );
        }
        upcase = LEFT;
        stack.push( vertices[i] );
        stack.push( vertices[i+1] );
        w = vertices[i+1];
        i++;
      } else {
        // Case R3
        Ray_2 ray( q, vertices[i] );
        Segment_2 seg( s_j_prev, s_j );

        Object_2 result = Intersect_2()( seg, ray );
        const Point_2 * ipoint = object_cast<Point_2>(&result);

        CGAL_assertion( ipoint != nullptr );

        u = *ipoint;
        if ( stack.top() != u ) {
          stack.push( u );
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
  void scana(Size_type& i, Point_2& w, const Point_2& q) const {
    // Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
    Point_2 u;
    Size_type k = scan_edges( i, q, stack.top(), u, true );

    Orientation orient1 =
    traits->orientation_2_object()(q, vertices[k], vertices[k+1] );

    if ( orient1 == RIGHT_TURN ) {
      bool fwd = traits->
        collinear_are_ordered_along_line_2_object()(q, stack.top(), u );

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
      stack.push( u );
      if ( u != vertices[k+1] ) {
        stack.push( vertices[k+1] );
      }
      w = vertices[k+1];
    }
  }

  /*! Find the first edge interecting the segment (v_0, s_t) */
  void scanb(Size_type& i, Point_2& w) const {
    if ( i == vertices.size() - 1 ) {
      upcase = FINISH;
      return;
    }
    Point_2 u;
    Size_type k = scan_edges( i, stack.top(), vertices[0], u, false );
    if ( (k+1 == vertices.size()-1) && (vertices[0] == u) ) {
      // Case B1
      upcase = FINISH;
      stack.push( vertices[0] );
    } else {
      // Case B2
      upcase = RIGHT;
      i = k+1;
      w = u;
    }
  }

  /*! Finds the exit from a general front hidden window by finding the first
    vertex to the right of the ray defined by the query_point and w*/
  void scanc(Size_type& i, Point_2& w) const {
    Point_2 u;
    Size_type k = scan_edges( i, stack.top(), w, u, false );
    upcase = RIGHT;
    i = k+1;
    w = u;
  }

  /*! find the first edge intersecting the given window (s_t, w) */
  void scand(Size_type& i, Point_2& w) const {
    Point_2 u;
    Size_type k = scan_edges( i, stack.top(), w, u, false );
    upcase = LEFT;
    i = k+1;
    stack.push( u );
    if ( u != vertices[k+1] ) {
      stack.push( vertices[k+1] );
    }
    w = vertices[k+1];
  }

  /*! Scan edges v_i,v_{i+1},...,v_n, until find an edge intersecting given ray
    or given segment. is_ray = true -> ray, false -> segment.
    The intersection point is returned by u */
  Size_type scan_edges( Size_type i,
                        const Point_2& ray_begin,
                        const Point_2& ray_end,
                        Point_2& u,
                        bool is_ray ) const
  {
    Orientation old_orient = RIGHT_TURN;
    Ray_2 ray( ray_begin, ray_end );
    Segment_2 s2( ray_begin, ray_end );
    Size_type k;
    Object_2 result;
    for ( k = i; k+1 < vertices.size(); k++ ) {
      Orientation curr_orient = traits->orientation_2_object()(
                                                               ray_begin,
                                                               ray_end,
                                                               vertices[k+1] );
      if ( curr_orient != old_orient ) {
        // Orientation switch, an intersection may occur
        Segment_2 seg( vertices[k], vertices[k+1] );
        if ( is_ray ) {
          result = Intersect_2()( seg, ray );
          if(result)
            break;
        } else {
          result = Intersect_2()( seg, s2 );
          if(result)
            break;
        }
      }
      old_orient = curr_orient;
    }
    CGAL_assertion( k+1<vertices.size() );
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
