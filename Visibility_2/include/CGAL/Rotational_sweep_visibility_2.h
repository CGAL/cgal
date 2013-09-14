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
//

#ifndef CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H
#define CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/tags.h>
#include <CGAL/bounding_box.h>
#include <CGAL/enum.h>

#include <CGAL/Timer.h>
namespace CGAL {



template <typename Arrangement_2, typename RegularizationTag>
class Rotational_sweep_visibility_2 {



public:
  typedef Arrangement_2                                 Input_arrangement_2;
  typedef Arrangement_2                                 Output_arrangement_2;
  typedef typename Input_arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Vertex_const_handle         Vertex_const_handle;
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

  typedef RegularizationTag                             Regularization_tag;
  typedef CGAL::Tag_true                                Supports_general_polygon_tag;
  typedef CGAL::Tag_true                                Supports_simple_polygon_tag;

  //profile
  Timer timer;
  static double input_t;
  static double sweep_t;
  static double cut_from_butterfly_t;
  static double heap_insert_t;
  static double heap_remove_t;
  static double heap_swap_t;
  static double input_v_t;
  static double quicksort_t;

private:
  typedef std::vector<Point_2>          Points;
  typedef Vertex_const_handle   Vertex;
  typedef std::vector<Vertex>   Vertices;
  typedef Halfedge_const_handle   Edge;
  typedef std::vector<Edge> Edges;

  const Geometry_traits_2 *geom_traits;

  class Less_edge: public std::binary_function<Edge, Edge, bool> {
    const Geometry_traits_2* traits;
  public:
//    Less_edge(const Geometry_traits_2* traits):traits(traits) {}
    bool operator() (const Edge e1, const Edge e2) const {
      if (e1 != e2) {
        if (e1->source() != e2->source())
          return e1->source()->point() < e2->source()->point();
    //        return Visibility_2::compare_xy_2(traits, e1->source()->point(), e2->source()->point()) == SMALLER;
        else
          return e1->target()->point() < e2->target()->point();
    //        return Visibility_2::compare_xy_2(traits, e1->target()->point(), e2->target()->point()) == SMALLER;
        }
      else
        return false;
    }
  };

  class Less_vertex: public std::binary_function<Vertex, Vertex, bool> {
//    const Geometry_traits_2* geom_traits;
  public:
//    Compare_vertex(const Geometry_traits_2* traits):geom_traits(traits) {}
    bool operator() (const Vertex v1, const Vertex v2) const {
      if (v1 != v2)
//        return Visibility_2::compare_xy_2(geom_traits, v1->point(), v2->point()) == SMALLER;
        return v1->point() < v2->point();
      else
        return false;
    }
  };

  const Input_arrangement_2 *p_arr;
  Point_2         q;
  Points polygon;                       //visibility polygon
  std::map<Vertex, Vertices, Less_vertex> neighbors;  //vertex and its neighbours that are relevant to visibility polygon
  std::map<Vertex, Edges, Less_vertex>  incident_edges;
  std::map<Edge, int, Less_edge> edx;            //index of edge in the heap
  Edges  active_edges;    //a heap of edges that interset the current vision ray.

  Vertices vs;                            //angular sorted vertices
  bool is_vertex_query;
  bool is_edge_query;
  bool is_big_cone;                   //whether the angle of visibility_cone is greater than pi.

  Edges bad_edge;
  Vertex query_vertex;
  Point_2  source;                    //one end of visibility cone
  Point_2  target;                    //another end of visibility cone





public:
  Rotational_sweep_visibility_2(): p_arr(NULL), geom_traits(NULL) {}
  Rotational_sweep_visibility_2(const Input_arrangement_2& arr): p_arr(&arr) {
    geom_traits = p_arr->geometry_traits();
  }

  Face_handle compute_visibility(const Point_2& q, const Halfedge_const_handle e, Arrangement_2& arr_out) {
    arr_out.clear();
    bad_edge.clear();
    this->q = q;

    if (Visibility_2::compare_xy_2(geom_traits, q, e->target()->point())==EQUAL) {
      query_vertex = e->target();
      is_vertex_query = true;
      is_edge_query = false;
      source = e->source()->point();
      target = e->next()->target()->point();
      is_big_cone = CGAL::right_turn(source, q, target);

      typename Input_arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
      first = curr = e->target()->incident_halfedges();
      do {
        if (curr->face() == e->face())
          bad_edge.push_back(curr);
        else if (curr->twin()->face() == e->face())
          bad_edge.push_back(curr->twin());
      } while (++curr != first);
    }
    else {
      is_vertex_query = false;
      is_edge_query = true;
      source = e->source()->point();
      target = e->target()->point();
      bad_edge.push_back(e);
      is_big_cone = false;
    }
    visibility_region_impl(e->face(), q);

    timer.reset();
    timer.start();
    //Decide which inside of the visibility butterfly is needed.
    int source_idx(-1), target_idx(-1) ;
    for (int i = 0; i != polygon.size(); i++) {
      if ( Visibility_2::compare_xy_2(geom_traits, polygon[i], source)==EQUAL ) {
        source_idx = i;
      }
      else if ( Visibility_2::compare_xy_2(geom_traits, polygon[i], target)==EQUAL ) {
        target_idx = i;
      }
      if (source_idx != -1 && target_idx != -1)
        break;
    }
    int small_idx, big_idx;
    if ( source_idx < target_idx ) {
      small_idx = source_idx;
      big_idx = target_idx;
    }
    else {
      small_idx = target_idx;
      big_idx = source_idx;
    }
    int next_idx = small_idx + 1;
    bool is_between;
    if (CGAL::right_turn(source, q, target)) {
      is_between = false;
      while (next_idx != big_idx) {
        if (CGAL::left_turn(source, q, polygon[next_idx]) || CGAL::left_turn(q, target, polygon[next_idx])) {
          is_between = true;
          break;
        }
        next_idx++;
      }
    }
    else {
      is_between = true;
      while (next_idx != big_idx) {
        if (CGAL::right_turn(source, q, polygon[next_idx]) || CGAL::right_turn(q, target, polygon[next_idx])) {
          is_between = false;
          break;
        }
        next_idx++;
      }
    }
    timer.stop();
    cut_from_butterfly_t+=timer.time();

    typename Points::iterator first = polygon.begin() + small_idx;
    typename Points::iterator last = polygon.begin() + big_idx;
    if (is_between) {
      Points polygon_out(first, last+1);
      if (is_vertex_query)
        polygon_out.push_back(q);
      Visibility_2::report_while_handling_needles_<Rotational_sweep_visibility_2>(geom_traits, q, polygon_out, arr_out);
    }
    else {
      Points polygon_out(polygon.begin(), first+1);
      if (is_vertex_query) polygon_out.push_back(q);
      for (int i = big_idx; i != polygon.size(); i++) {
        polygon_out.push_back(polygon[i]);
      }
      Visibility_2::report_while_handling_needles_<Rotational_sweep_visibility_2>(geom_traits, q, polygon_out, arr_out);
    }

    conditional_regularize(arr_out, Regularization_tag());

    if (arr_out.faces_begin()->is_unbounded())
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();

  }

  Face_handle compute_visibility(const Point_2& q, const Face_const_handle f, Output_arrangement_2& arr_out) {
    arr_out.clear();
    this->q = q;
    is_vertex_query = false;
    is_edge_query = false;

    visibility_region_impl(f, q);
    Visibility_2::report_while_handling_needles_<Rotational_sweep_visibility_2>(geom_traits, q, polygon, arr_out);
    conditional_regularize(arr_out, Regularization_tag());
    if (arr_out.faces_begin()->is_unbounded())
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }

bool is_attached() {
  return (p_arr != NULL);
}

void attach(const Input_arrangement_2& arr) {
  p_arr = &arr;
  geom_traits = p_arr->geometry_traits();
}

void detach() {
  p_arr = NULL;
  geom_traits = NULL;
}

const Input_arrangement_2& arr() {
  return *p_arr;
}


private:
  bool do_intersect_ray(const Point_2& q,
                    const Point_2& dp,
                    const Point_2& p1,
                    const Point_2& p2) {
    return (CGAL::orientation(q, dp, p1) != CGAL::orientation(q, dp, p2) && CGAL::orientation(q, p1, dp) == CGAL::orientation(q, p1, p2));
  }

  void funnel(int i, int j) {
    Vertices right, left;
    bool block_left(false), block_right(false);
    Vertex former = vs[i], neib;
    for (int l=i; l<j; l++) {
      bool left_v(false), right_v(false), has_predecessor(false);
      for (int k=0; k<neighbors[vs[l]].size(); k++) {
        neib= neighbors[vs[l]][k];
        if ( neib == former )  {
          has_predecessor = true;
          continue;
        }
        if (CGAL::left_turn(q, vs[l]->point(), neib->point()))
          left_v = true;
        else
          right_v = CGAL::right_turn(q, vs[l]->point(), neib->point());
      }
      if (has_predecessor) {
          block_left = block_left || left_v;
          block_right = block_right || right_v;
      }
      else {
          block_left = left_v;
          block_right = right_v;
      }
      if (block_left && block_right) {
        right.push_back(vs[l]);
        break;
      }
      else {
        if (block_left)
          left.push_back(vs[l]);
        else
          right.push_back(vs[l]);
      }
      former = vs[l];
    }
    for (int l=0; l!=right.size(); l++)
      vs[i+l] = right[l];
    for (int l=0; l!=left.size(); l++)
      vs[i+l+right.size()] = left[left.size()-1-l];
  }

  void visibility_region_impl(const Face_const_handle f, const Point_2& q) {
    vs.clear();
    polygon.clear();
    active_edges.clear();
    incident_edges.clear();
    neighbors.clear();
    edx.clear();

    Edges good_edges;
    Input_arrangement_2 bbox;
    if (is_vertex_query || is_edge_query)
      input_face(f, good_edges, bbox);
    else
      input_face(f);
    //initiation of vision ray
    Vector_2 dir;
    if (Direction_2(-1, 0) < Direction_2(Vector_2(q, vs.back()->point())))
    {
      dir = Vector_2(1, 0) + Vector_2(q, vs.back()->point());
    }
    else {
        dir = Vector_2(0, -1);
    }

    Point_2 dp = q + dir;

//    std::vector<Edge> heapc;
//    heapc.clear();

    //initiation of active_edges
    if (is_vertex_query || is_edge_query) {
      for (int i=0; i!=good_edges.size(); i++) {
        if (do_intersect_ray(q, dp, good_edges[i]->source()->point(), good_edges[i]->target()->point())) {
          heap_insert(good_edges[i]);
//          heapc.push_back(good_edges[i]);
        }
      }
//      std::make_heap(heapc.begin(), heapc.end(), Is_closer(q, geom_traits));
//      for (int i=0; i!=heap.size(); i++) {
//        edx[heap[i]] = i;
//      }
//      compare_heap(heap, heapc);
    }
    else {
      Ccb_halfedge_const_circulator curr = f->outer_ccb();
      Ccb_halfedge_const_circulator circ = curr;
      do {
//        Point_2 p1 = curr->target()->point();
//        Point_2 p2 = curr->source()->point();
        if (do_intersect_ray(q, dp, curr->target()->point(), curr->source()->point()))
          heap_insert(curr);
      } while (++curr != circ);

      typename Arrangement_2::Hole_const_iterator hi;
      for (hi = f->holes_begin(); hi != f->holes_end(); ++hi) {
        Ccb_halfedge_const_circulator curr = *hi, circ = *hi;
        do {
//          Point_2 p1 = curr->target()->point();
//          Point_2 p2 = curr->source()->point();
          if (do_intersect_ray(q, dp, curr->target()->point(), curr->source()->point()))
            heap_insert(curr);
        } while (++curr != circ);
      }
    }

    //angular sweep begins
    for (int i=0; i!=vs.size(); i++) {
      Vertex vh = vs[i];
      Edge closest_e = active_edges.front();   //save the closest edge;
      int insert_cnt(0), remove_cnt(0);
      Edges& edges = incident_edges[vh];
      Edges insert_es, remove_es;

      for (int j=0; j!=edges.size(); j++) {
        Edge e = edges[j];
//        Orientation o=Visibility_2::orientation_2(geom_traits, q, dp, nei);
/*        if (o==RIGHT_TURN ||
            (o==COLLINEAR && i>0 && Visibility_2::compare_xy_2(geom_traits, nei, vs[i-1])==EQUAL))*/

        if (edx.count(e)){
          remove_es.push_back(e);
        }
        else {
          insert_es.push_back(e);
        }

      }
      insert_cnt = insert_es.size();
      remove_cnt = remove_es.size();
      if (remove_es.size()==1 && insert_es.size()==1) {
        int remove_idx = edx[remove_es.front()];
        active_edges[remove_idx] = insert_es.front();
        edx[insert_es.front()] = remove_idx;
        edx.erase(remove_es.front());
      }
      else {
        for (int j=0; j!=remove_es.size(); j++) {
          heap_remove(edx[remove_es[j]]);
        }
        for (int j=0; j!=insert_es.size(); j++) {
          heap_insert(insert_es[j]);
        }
      }

      if (closest_e != active_edges.front()) {
        //when the closest edge changed
        if (remove_cnt > 0 && insert_cnt > 0) {
          //some edges are added and some are deleted, which means the vertice sweeped is a vertice of visibility polygon.
          update_visibility(vh->point());
        }
        if (remove_cnt == 0 && insert_cnt > 0) {
          //only add some edges, means the view ray is blocked by new edges.
          //therefore first add the intersection of view ray and former closet edge, then add the vertice sweeped.
          update_visibility(ray_seg_intersection(q,
                                                 vh->point(),
                                                 closest_e->target()->point(),
                                                 closest_e->source()->point()));
          update_visibility(vh->point());
        }
        if (remove_cnt > 0 && insert_cnt == 0) {
          //only delete some edges, means some block is moved and the view ray can reach the segments after the block.
          update_visibility(vh->point());
          update_visibility(ray_seg_intersection(q,
                                                 vh->point(),
                                                 active_edges.front()->target()->point(),
                                                 active_edges.front()->source()->point()));
        }
      }
    }
  }

//  Edge create_pair(const Point_2& p1, const Point_2& p2) const{
//    assert(p1 != p2);
//    if (Visibility_2::compare_xy_2(geom_traits, p1, p2)==SMALLER)
//      return Edge(p1, p2);
//    else
//      return Edge(p2, p1);
//  }

  void heap_insert(const Edge e) {
    active_edges.push_back(e);
    int i = active_edges.size()-1;
    edx[e] = i;
    int parent = (i-1)/2;
    while (i!=0 && is_closer(q, active_edges[i], active_edges[parent])){
      heap_swap(i, parent);
      i = parent;
      parent = (i-1)/2;
    }
  }

  void heap_remove(int i) {
    edx.erase(active_edges[i]);
    if (i == active_edges.size()-1)
    {
      active_edges.pop_back();
    }
    else {
      active_edges[i] = active_edges.back();
      edx[active_edges[i]] = i;
      active_edges.pop_back();
      int i_before_swap = i;

      int parent = (i-1)/2;
      while (i!=0 && is_closer(q, active_edges[i], active_edges[parent])){
        heap_swap(i, parent);
        i = parent;
        parent = (i-1)/2;
      }
      if (i==i_before_swap) {
        bool swapped;
        do {
          int left_son = i*2+1;
          int right_son = i*2+2;
          int closest_idx = i;
          if (left_son < active_edges.size() && is_closer(q, active_edges[left_son], active_edges[i])) {
            closest_idx = left_son;
          }
          if (right_son < active_edges.size() && is_closer(q, active_edges[right_son], active_edges[closest_idx])) {
            closest_idx = right_son;
          }
          swapped = false;
          if (closest_idx != i) {
            heap_swap(i, closest_idx);
            i = closest_idx;
            swapped = true;
          }
        } while(swapped);
      }
    }

  }

  void heap_swap(int i, int j) {
    edx[active_edges[i]] = j;
    edx[active_edges[j]] = i;
    Edge temp = active_edges[i];
    active_edges[i] = active_edges[j];
    active_edges[j] = temp;
  }
  void print_point(const Point_2& p) {
    std::cout<<p.x()<<','<<p.y()<<std::endl;
  }

  bool is_closer(const Point_2& q,
                 const Edge& e1,
                 const Edge& e2) {
    const Point_2& s1=e1->target()->point(),
        t1=e1->source()->point(),
        s2=e2->target()->point(),
        t2=e2->source()->point();
    Orientation e1q = Visibility_2::orientation_2(geom_traits, s1, t1, q);
    switch (e1q)
    {
    case COLLINEAR:
      if (Visibility_2::collinear(geom_traits, q, s2, t2)) {
        //q is collinear with e1 and e2.
        return (Visibility_2::less_distance_to_point_2(geom_traits, q, s1, s2)
                || Visibility_2::less_distance_to_point_2(geom_traits, q, t1, t2));
      }
      else {
        //q is not collinear with e2. q is collinear with e1.
        if (Visibility_2::collinear(geom_traits, s2, t2, s1))
          return (Visibility_2::orientation_2(geom_traits, s2, t2, q)
                  == Visibility_2::orientation_2(geom_traits, s2, t2, t1));
        else
          return (Visibility_2::orientation_2(geom_traits, s2, t2, q)
                  == Visibility_2::orientation_2(geom_traits, s2, t2, s1));
      }
    case RIGHT_TURN:
      switch (Visibility_2::orientation_2(geom_traits, s1, t1, s2)) {
      case COLLINEAR:
        return Visibility_2::orientation_2(geom_traits, s1, t1, t2)!=e1q;
      case RIGHT_TURN:
        if (Visibility_2::orientation_2(geom_traits, s1, t1, t2) == LEFT_TURN)
          return Visibility_2::orientation_2(geom_traits, s2, t2, q)
              == Visibility_2::orientation_2(geom_traits, s2, t2, s1);
        else
          return false;
      case LEFT_TURN:
        if (Visibility_2::orientation_2(geom_traits, s1, t1, t2) == RIGHT_TURN)
          return Visibility_2::orientation_2(geom_traits, s2, t2, q)
              == Visibility_2::orientation_2(geom_traits, s2, t2, s1);
        else
          return true;
      }
    case LEFT_TURN:
      switch (Visibility_2::orientation_2(geom_traits, s1, t1, s2)) {
      case COLLINEAR:
        return Visibility_2::orientation_2(geom_traits, s1, t1, t2)!=e1q;
      case LEFT_TURN:
        if (Visibility_2::orientation_2(geom_traits, s1, t1, t2) == RIGHT_TURN)
          return Visibility_2::orientation_2(geom_traits, s2, t2, q)
              == Visibility_2::orientation_2(geom_traits, s2, t2, s1);
        else
          return false;
      case RIGHT_TURN:
        if (Visibility_2::orientation_2(geom_traits, s1, t1, t2) == LEFT_TURN)
          return Visibility_2::orientation_2(geom_traits, s2, t2, q)
              == Visibility_2::orientation_2(geom_traits, s2, t2, s1);
        else
          return true;
      }
    }
  }

  Point_2 ray_seg_intersection(
      const Point_2& q, const Point_2& dp,  // the ray
      const Point_2& s, const Point_2& t)   // the segment
  {
    if (CGAL::collinear(q, dp, s)) {
      if (CGAL::collinear(q, dp, t)) {
        if (CGAL::compare_distance_to_point(q, s, t)==CGAL::SMALLER)
          return s;
        else
          return t;
      }
      else
        return s;
    }
    Ray_2 ray(q,dp);
    Segment_2 seg(s,t);
    CGAL::Object result = CGAL::intersection(ray, seg);
    return *(CGAL::object_cast<Point_2>(&result));
//    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
//      return *ipoint;
//    }
//    else {
//      if (const Segment_2 *iseg = CGAL::object_cast<Segment_2 >(&result)) {
//        switch (CGAL::compare_distance_to_point(ray.source(), iseg->source(), iseg->target())) {
//        case (CGAL::SMALLER):
//          return iseg->source();
//          break;
//        case (CGAL::LARGER) :
//          return iseg->target();
//          break;
//        }
//      } else {
//        assert(false);
//      }
//    }
  }

  void update_visibility(const Point_2& p){
    if (polygon.empty())
      polygon.push_back(p);
    else
    {
      if (Visibility_2::compare_xy_2(geom_traits, polygon.back(), p) != EQUAL) {
        polygon.push_back(p);
      }
    }
  }
  class Is_sweeped_first:public std::binary_function<Vertex, Vertex, bool> {
    const Point_2& q;
    const Geometry_traits_2* geom_traits;
  public:
    Is_sweeped_first(const Point_2& q, const Geometry_traits_2* traits):q(q){
      geom_traits = traits;
    }
    bool operator() (const Vertex v1, const Vertex v2) const {
      const Point_2& p1 = v1->point();
      const Point_2& p2 = v2->point();
      int qua1 = quadrant(q, p1);
      int qua2 = quadrant(q, p2);
      if (qua1 < qua2)
        return true;
      if (qua1 > qua2)
        return false;
      if (collinear(q, p1, p2))
        return (CGAL::compare_distance_to_point(q, p1, p2) == CGAL::SMALLER);
      else
        return CGAL::right_turn(p1, q, p2);
    }
    int quadrant(const Point_2& o, const Point_2& p) const {
      typename Geometry_traits_2::Compare_x_2 compare_x = geom_traits->compare_x_2_object();
      typename Geometry_traits_2::Compare_y_2 compare_y = geom_traits->compare_y_2_object();

      Comparison_result dx = compare_x(p, o);
      Comparison_result dy = compare_y(p, o);

      if (dx==LARGER && dy!=SMALLER)
        return 1;
      if (dx!=LARGER && dy==LARGER)
        return 2;
      if (dx==SMALLER && dy!=LARGER)
        return 3;
      if (dx!=SMALLER && dy==SMALLER)
        return 4;
      return 0;
    }

  };

  //when query is in face, every edge is good.
  void input_neighbor_f( const Halfedge_const_handle e) {
    Vertex v = e->target();
    if (!neighbors.count(v))
      vs.push_back(v);
      neighbors[v].push_back(e->source());
      neighbors[v].push_back(e->next()->target());
      incident_edges[v].push_back(e);
      incident_edges[v].push_back(e->next());
  }

  bool is_in_cone(const Point_2& p) const{
    if (is_big_cone)
      return (!CGAL::right_turn(source, q, p)) || (!CGAL::left_turn(target, q, p));
    else
      return (!CGAL::right_turn(source, q, p)) && (!CGAL::left_turn(target, q, p));
  }

  //for vertex and edge query: the visibility is limited in a cone.
  void input_edge(const Halfedge_const_handle e,
                  Edges& good_edges) {
    for (int i=0; i<bad_edge.size(); i++)
      if (e == bad_edge[i])
        return;

    Vertex v1 = e->target();
    Vertex v2 = e->source();
    if (is_in_cone(v1->point()) || is_in_cone(v2->point()) || do_intersect_ray(q, source, v1->point(), v2->point())) {
      good_edges.push_back(e);
      if (!neighbors.count(v1))
        vs.push_back(v1);
      neighbors[v1].push_back(v2);
      incident_edges[v1].push_back(e);
      if (!neighbors.count(v2))
        vs.push_back(v2);
      neighbors[v2].push_back(v1);
      incident_edges[v2].push_back(e);
    }
  }

  //for face query: traverse the face to get all edges and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh)
  {
    Ccb_halfedge_const_circulator curr = fh->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      assert(curr->face() == fh);
      input_neighbor_f(curr);
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator curr = *hi, circ = *hi;
      do {
        assert(curr->face() == fh);
        input_neighbor_f(curr);
      } while (++curr != circ);
    }

    std::sort(vs.begin(), vs.end(), Is_sweeped_first(q, geom_traits));

    for (int i=0; i!=vs.size(); i++) {
      int j = i+1;
      while (j != vs.size()) {
        if (!CGAL::collinear(q, vs[i]->point(), vs[j]->point()))
          break;
        j++;
      }
      if (j-i>1)
        funnel(i, j);
      i = j-1;
    }
  }
  //for vertex or edge query: traverse the face to get all edges and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh,
                   Edges& good_edges,
                   Input_arrangement_2& bbox)
  {
//    timer.reset();
//    timer.start();

    Ccb_halfedge_const_circulator curr = fh->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      assert(curr->face() == fh);
      input_edge(curr, good_edges);
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator curr = *hi, circ = *hi;
      do {
        assert(curr->face() == fh);
        input_edge(curr, good_edges);
      } while (++curr != circ);
    }
    //todo
    Points points;
    for (int i=0; i<vs.size(); i++) {
      points.push_back(vs[i]->point());
    }
    points.push_back(q);
    typename Geometry_traits_2::Iso_rectangle_2 bb = bounding_box(points.begin(), points.end());
//    points.pop_back();

    Number_type xmin, xmax, ymin, ymax;
    typename Geometry_traits_2::Compute_x_2 compute_x = geom_traits->compute_x_2_object();
    typename Geometry_traits_2::Compute_y_2 compute_y = geom_traits->compute_y_2_object();
    xmin = compute_x(bb.min())-1;
    ymin = compute_y(bb.min())-1;
    xmax = compute_x(bb.max())+1;
    ymax = compute_y(bb.max())+1;
    Point_2 box[4] = {Point_2(xmin, ymin), Point_2(xmax, ymin),
                      Point_2(xmax, ymax), Point_2(xmin, ymax)};
    Halfedge_handle e1 = bbox.insert_in_face_interior(Segment_2(box[0], box[1]), bbox.unbounded_face());
    Halfedge_handle e2 = bbox.insert_from_left_vertex(Segment_2(box[1], box[2]), e1->target());
    Halfedge_handle e3 = bbox.insert_from_right_vertex(Segment_2(box[2], box[3]), e2->target());
    bbox.insert_at_vertices(Segment_2(box[0], box[3]), e1->source(), e3->target());

//    Face_const_handle f = e1->face();
    circ = curr = e1->face()->outer_ccb();
    do {
      Vertex v = curr->target();
      vs.push_back(v);
      neighbors[v].push_back(curr->source());
      neighbors[v].push_back(curr->next()->target());
      incident_edges[v].push_back(curr);
      incident_edges[v].push_back(curr->next());
      good_edges.push_back(curr);
    } while(++curr != circ);


//    for (int i=0; i<4; i++) {
//      vs.push_back(box[i]);
//      neighbors[box[i]].push_back(box[(i+3)%4]);
//      neighbors[box[i]].push_back(box[(i+1)%4]);
//      good_edges.push_back(create_pair(box[i], box[(i+1)%4]));
//    }

    std::sort(vs.begin(), vs.end(), Is_sweeped_first(q, geom_traits));

    for (int i=0; i!=vs.size(); i++) {
      int j = i+1;
      while (j != vs.size()) {
        if (!CGAL::collinear(q, vs[i]->point(), vs[j]->point()))
          break;
        j++;
      }
      if (j-i>1)
        funnel(i, j);
      i = j-1;
    }
  }


//  void build_arr(const Pvec& polygon, Output_arrangement_2& arr ) {
//      for (int i = 0; i != polygon.size()-1; i++ ) {
//          CGAL::insert(arr, Segment_2(polygon[i], polygon[i+1]));
//      }
//      //print_vectex(polygon);
//      CGAL::insert(arr, Segment_2(polygon.front(), polygon.back()));
//  }


  void conditional_regularize(Output_arrangement_2& arr_out, CGAL::Tag_true) {
    regularize_output(arr_out);
  }

  void conditional_regularize(Output_arrangement_2& arr_out, CGAL::Tag_false) {
    //do nothing
  }

  void regularize_output(Output_arrangement_2& arr_out) {
    typename Output_arrangement_2::Edge_iterator e_itr;
    for (e_itr = arr_out.edges_begin();
         e_itr != arr_out.edges_end();
         e_itr++) {

      Halfedge_handle he = e_itr;
      Halfedge_handle he_twin = he->twin();
      if (he->face() == he_twin->face()) {
        arr_out.remove_edge(he);
      }
    }
  }

};
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::sweep_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::cut_from_butterfly_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::input_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::heap_insert_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::heap_remove_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::heap_swap_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::input_v_t = 0;
template <typename Arrangement_2, typename RegularizationTag>
double CGAL::Rotational_sweep_visibility_2<Arrangement_2, RegularizationTag>::quicksort_t = 0;

} // end namespace CGAL



#endif


