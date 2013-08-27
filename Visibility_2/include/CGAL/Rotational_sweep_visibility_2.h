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
#include <list>
#include <map>
#include <utility>
#include <CGAL/Arrangement_2.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>

namespace CGAL {

template <typename Arrangement_2, typename RegularizationTag>
class Rotational_sweep_visibility_2 {
public:
  typedef Arrangement_2                                 Input_arrangement_2;
  typedef Arrangement_2                                 Output_arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Vertex_const_handle         Vertex_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                        Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

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

  enum Intersection_type { UNBOUNDED, CORNER, INNER };

  Rotational_sweep_visibility_2(): p_arr(NULL) {}
  Rotational_sweep_visibility_2(Input_arrangement_2& arr): p_arr(&arr) {}


  Face_handle compute_visibility(const Point_2& q, const Halfedge_const_handle e, Arrangement_2& out_arr) {
    out_arr.clear();
    this->q = q;
    Point_2 source, target;
    if (q == e->target()->point()) {
      is_vertex_query = true;
      is_edge_query = false;
      source = e->source()->point();
      target = e->next()->target()->point();
    }
    else {
      is_vertex_query = false;
      is_edge_query = true;
      source = e->source()->point();
      target = e->target()->point();
    }
    visibility_region_impl(e->face(), q, source, target);

    //Decide which inside of the visibility butterfly is needed.
    int source_i, target_i ;
    for (int i = 0; i != polygon.size(); i++) {
      if ( polygon[i]== source ) {
          source_i = i;
      }
      else if ( polygon[i] == target ) {
          target_i = i;
      }
    }
    int small, big;
    if ( source_i < target_i ) {
      small = source_i;
      big = target_i;
    }
    else {
      small = target_i;
      big = source_i;
    }
    int next_i = small + 1;
    bool is_between;
    if (CGAL::right_turn(source, q, target)) {
      is_between = false;
      while (next_i != big) {
        if (CGAL::left_turn(source, q, polygon[next_i]) || CGAL::left_turn(q, target, polygon[next_i])) {
          is_between = true;
          break;
        }
        next_i++;
      }
    }
    else {
      is_between = true;
      while (next_i != big) {
        if (CGAL::right_turn(source, q, polygon[next_i]) || CGAL::right_turn(q, target, polygon[next_i])) {
          is_between = false;
          break;
        }
        next_i++;
      }
    }
    typename std::vector<Point_2>::iterator first = polygon.begin() + small;
    typename std::vector<Point_2>::iterator last = polygon.begin() + big;
    if (is_between) {
      std::vector<Point_2> polygon1(first, last+1);
      if (is_vertex_query) polygon1.push_back(q);
      build_arr(polygon1, out_arr);
    }
    else {
      std::vector<Point_2> polygon1(polygon.begin(), first+1);
      if (is_vertex_query) polygon1.push_back(q);
      for (int i = big; i != polygon.size(); i++) {
        polygon1.push_back(polygon[i]);
      }
      build_arr(polygon1, out_arr);
    }

    conditional_regularize(out_arr, Regularization_tag());

    if (out_arr.faces_begin()->is_unbounded())
      return ++out_arr.faces_begin();
    else
      return out_arr.faces_begin();

  }

  Face_handle compute_visibility(const Point_2& q, const Face_const_handle f, Output_arrangement_2& out_arr) {
    out_arr->clear();
    this->q = q;
    is_vertex_query = false;
    is_edge_query = false;

    visibility_region_impl(f, q, q, q);
    build_arr(polygon, out_arr);

    conditional_regularize(out_arr, Regularization_tag());
    if (out_arr.faces_begin()->is_unbounded())
      return ++out_arr.faces_begin();
    else
      return out_arr.faces_begin();
  }

bool is_attached() {
  return (p_arr != NULL);
}

void attach(Input_arrangement_2& arr) {
  p_arr = &arr;
//  geom_traits = p_arr->geometry_traits();
}

void detach() {
  p_arr = NULL;
//  geom_traits = NULL;
  vs.clear();
}

const Input_arrangement_2& arr() {
  return *p_arr;
}


private:
  //members
  typedef std::vector<Point_2> Pvec;
  typedef std::pair<Point_2, Point_2> Pair;
  const Input_arrangement_2 *p_arr;
  bool            attach_tag;
  Point_2         q;
  std::vector<Point_2> polygon;   //visibility polygon
  std::map<Point_2, Pvec> vmap;   //vertex and two edges incident to it that might block vision
  std::list<Point_2> vs;          //angular sorted vertices
  bool is_vertex_query;
  bool is_edge_query;

  int quadrant(Number_type x, Number_type y) {
    if (x>0 && y>=0)
      return 1;
    if (x<=0 && y>0)
      return 2;
    if (x<0 && y<=0)
      return 3;
    if (x>=0 && y<0)
      return 4;
    return 0;
  }

  bool if_intersect(const Point_2& q,
                    const Point_2& d_p,
                    const Point_2& p1,
                    const Point_2& p2) {
    if (CGAL::collinear(q, d_p, p1))
      return (quadrant(p1.x()-q.x(), p1.y()-q.y()) == quadrant((d_p.x()-q.x(), d_p.y()-q.y())));

    if (CGAL::collinear(q, d_p, p2))
      return (quadrant(p2.x()-q.x(), p2.y()-q.y()) == quadrant((d_p.x()-q.x(), d_p.y()-q.y())));

    return (CGAL::orientation(q, d_p, p1) != CGAL::orientation(q, d_p, p2) && CGAL::orientation(q, p1, d_p) == CGAL::orientation(q, p1, p2));

  }

  //compute visibility region between qa and qb
  void visibility_region_impl(Face_const_handle f, const Point_2& q, const Point_2& a, const Point_2& b) {
    input_face(f, q, a, b);
    vs.sort(compare_angle);

//        debug
//        for (int i = 0; i<vertices.size(); i++) {
//          print(vertices[i]->point());
//        }

    //initiation of vision ray
    Vector_2 dir;
    if (Direction_2(-1, 0) < Direction_2(Vector_2(q, vs.back())))
    {
        dir = Vector_2(1, 0) + Vector_2(q, vs.back());
    }
    else {
        dir = Vector_2(0, -1);
    }

    Point_2 direct_p(q.x()+dir.x(), q.y()+dir.y());
    //initiation of active_edges
    std::vector<Pair> active_edges;
    Ccb_halfedge_const_circulator curr = f->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      Point_2 p1 = curr->target()->point();
      Point_2 p2 = curr->source()->point();
      if (q != p1 && q != p2 && if_intersect(q, direct_p, p1, p2))
        active_edges.push_back(create_pair(p1, p2));
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = f->holes_begin(); hi != f->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator c1 = *hi, c2 = *hi;
      do {
        Point_2 p1 = c1->target()->point();
        Point_2 p2 = c1->source()->point();
        if (q != p1 && q != p2 && if_intersect(q, direct_p, p1, p2))
          active_edges.push_back(create_pair(p1, p2));
      } while (++c1 != c2);
    }

    //angular sweep begins

//    typename Pvec::iterator vit = vs.begin(), begin_it, end_it;

//    while (vit != vs.end())
//    {
//      Point_2 right_p, left_p, mid_p;
//      begin_it = vit;
//      end_it = vit + 1;
//      right_p = intersection_point(direct_p, active_edges[0]);

//      //find end_it such that all vertices between begin_it and end_it(not included) are collinear with query point.
//      while (end_it != vs.end()) {
//        if (!CGAL::collinear(q, *begin_it, *end_it))
//          break;
//        insert_edge(*end_it, active_edges, direct_p);
//        end_it++;
//      }

////                std::cout<<"after adding\n";
////                print_edges(active_edges);
//      mid_p = intersection_point(direct_p, active_edges[0]);
//      Pvec collinear_vertices;
//      Intersection_type i_type = needle(active_edges, direct_p, collinear_vertices);
//      switch (i_type) {
//      case UNBOUNDED :
//              //todo:this part is not finished.
//              //remove right and collinear;
//              remove_edges(active_edges, direct_p);
//              update_visibility(right_p, polygon);
//              update_visibility(mid_p, polygon);
//              //todo CGAL::insert_curve();
//              if (!active_edges.empty()) {
//                  left_p = intersection_point(direct_p, active_edges[0]);
//                  update_visibility(left_p, polygon);
//              }
//              break;
//            case CORNER :
//                //remove right and collinear;
//                remove_edges(active_edges, direct_p);

////                    std::cout<<"after removing\n";
////                    print_edges(active_edges);

//                left_p = intersection_point(direct_p, active_edges[0]);
//                update_visibility(right_p, polygon);
//                if (right_p == collinear_vertices[0]) {
//                  insert_needle(collinear_vertices, polygon, true);
//                }
//                else if (left_p == collinear_vertices[0]) {
//                  insert_needle(collinear_vertices, polygon, false);
//                }
//                update_visibility(left_p, polygon);
//                break;
//            case INNER :
//                //remove right and collinear;
//                remove_edges(active_edges, direct_p);
//                if (collinear_vertices.size() < 2) {
//                    //this means mid_p = left_p = right_p = furthest_p. no new vertex is found.
//                }
//                else {
//                  //debug
////                      std::cout<<"print a needle:\n";
////                      print(collinear_vertices);
////                      std::cout<<"the left_p is "<<left_p.x()<<' '<<left_p.y()<<std::endl;
////                      std::cout<<"the right_p is "<<right_p.x()<<' '<<right_p.y()<<std::endl;
////                      std::cout<<"the front_p is "<<collinear_vertices[0].x()<<' '<<collinear_vertices[0].y()<<std::endl;
//                  left_p = intersection_point(direct_p, active_edges[0]);
//                  update_visibility(right_p, polygon);
//                  if (right_p == collinear_vertices[0]) {
//                    insert_needle(collinear_vertices, polygon, true);
//                  }
//                  else if (left_p == collinear_vertices[0]) {
//                    insert_needle(collinear_vertices, polygon, false);
//                  }
//                  update_visibility(left_p, polygon);
//                }
//                break;
//            }
//        }
//        vit = end_it;
  }

//todo add edge location record
  void heap_insert(std::vector<Pair>& heap, const Pair& e, const Point_2& dp) {
    heap.push_back(e);
    int i = heap.size()-1;
    int parent = (i-1)/2;
    while (i!=0 && is_closer(heap[i], heap[parent], dp)){
      heap_swap(heap, i, parent);
      i = parent;
      parent = (i-1)/2;
    }
  }

  void heap_remove(std::vector<Pair>& heap, int i, const Point_2& dp) {
    heap[i] = heap.back();
    heap.pop_back();
    bool swapped;
    do {
      int left_son = i*2+1;
      int right_son = i*2+2;
      int closest = i;
      if (left_son < heap.size() && is_closer(heap[left_son], heap[i], dp)) {
        closest = left_son;
      }
      if (right_son < heap.size() && is_closer(heap[right_son], heap[closest], dp)) {
        closest = right_son;
      }
      swapped = false;
      if (closest != i) {
        heap_swap(heap, i, closest);
        i = closest;
        swapped = true;
      }
    } while(swapped);
  }

  void heap_swap(std::vector<Pair>& heap, int i, int j) {

  }

  bool is_closer(Pair e1, Pair e2, Point_2 dp) {

  }

  Point_2 intersection_point(Ray_2 ray, Segment_2 seg )
  {

      CGAL::Object result = CGAL::intersection(ray, seg);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          return *ipoint;
      } else
          //if result is a segment, return the end closer to the source of ray.
          if (const Segment_2 *iseg = CGAL::object_cast<Segment_2 >(&result)) {
              switch (CGAL::compare_distance_to_point(ray.source(), iseg->source(), iseg->target())) {
              case (CGAL::SMALLER):
                  return iseg->source();
                  break;
              case (CGAL::LARGER) :
                  return iseg->target();
                  break;
              }

          } else {
              // if no intersection, return the source of ray.
              return ray.source();
          }
  }


  //insert newly-discovered edges into active_edges according to its intersection with the view ray.
  void insert_halfedge(std::vector<Halfedge_const_handle>& active_edges, const Ray_2& ray, Halfedge_const_handle edge)
  {
      Point_2 cross_of_e = intersection_point(ray, edge);
      if (cross_of_e != ray.source())
      {
          typename std::vector<Halfedge_const_handle>::iterator curr = active_edges.begin();
          while (curr != active_edges.end())
          {

              Point_2 cross_of_curr = intersection_point(ray, *curr);
              if (CGAL::compare_distance_to_point(ray.source(), cross_of_e, cross_of_curr) == CGAL::SMALLER)
                  break;
              if (cross_of_curr == cross_of_e && is_closer(ray, edge, *curr))
                  break;
              ++curr;
          }
          active_edges.insert(curr, edge);
      }
  }

  bool compare_angle(const Point_2& p1, const Point_2& p2)
  {
    Direction_2 d1(Ray_2(q, p1));
    Direction_2 d2(Ray_2(q, p2));
    if (d1 < d2)
      return true;
    if (d1 > d2)
      return false;
    return (CGAL::compare_distance_to_point(q, p1, p2) == CGAL::SMALLER);
  }


  //traverse the face to get all edges and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh,
                   const Point_2& q,
                   const Point_2& a,
                   const Point_2& b)
  {
    Ccb_halfedge_const_circulator curr = fh->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      Point_2 v = curr->target()->point();
      if (v == q)
        continue;
      vs.push_back(v);
      Pvec neighbor;
      if (curr->source()->point() != q)
        neighbor.push_back(curr->source()->point());
      if (curr->next()->target()->point() != q)
        neighbor.push_back(curr->next()->target()->point());
      vmap[v] = neighbor;
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      typename Arrangement_2::Ccb_halfedge_const_circulator c1 = *hi, c2 = *hi;
      do {
        Point_2 v = curr->target()->point();
        if (v == q)
          continue;
        vs.push_back(v);
        Pvec neighbor;
        if (c1->source()->point() != q)
          neighbor.push_back(c1->source()->point());
        if (c1->next()->target()->point() != q)
          neighbor.push_back(c1->next()->target()->point());
        vmap[v] = neighbor;
      } while (++c1 != c2);
    }
    if (q != a) {
      Number_type xmin, xmax, ymin, ymax;
      Point_2 q1 = vs[0];
      xmax = xmin = q1.x();
      ymin = ymax = q1.y();
      for (int i=0; i<vs.size(); i++) {
        Point_2 q1 = vs[i];
        if (q1.x() < xmin)
          xmin = q1.x();
        if (q1.x() > xmax)
          xmax = q1.x();
        if (q1.y() < ymin)
          ymin = q1.y();
        if (q1.y() > ymax)
          ymax = q1.y();
      }
      xmin -= 10;
      xmax += 10;
      ymin -= 10;
      ymax += 10;
      Point_2 box[4] = {Point_2(xmin, ymin), Point_2(xmax, ymin), Point_2(xmax, ymax), Point_2(xmin, ymax)};
      for (int i=0; i<4; i++) {
        vs.push_back(box[i]);
        Pvec pvec;
        pvec.push_back(box[(i-1)%4]);
        pvec.push_back(box[(i+1)%4]);
        vmap[box[i]] = pvec;
      }
    }
  }




  //insert new vertice to polygon. before insertion, check if this vertice has been added before.
  void update_visibility(const Point_2 p, std::vector<Point_2>& polygon){
      if (polygon.empty())
          polygon.push_back(p);
      else
      {
          if (polygon.back() != p) {
              polygon.push_back(p);
          }
      }
  }

  void insert_needle(const std::vector<Point_2>& points, std::vector<Point_2>& polygon, bool is_right_close){
    if (is_right_close) {
      for (int i = 0; i != points.size(); i++) {
        update_visibility(points[i], polygon);
      }
    }
    else {
      for (int i = points.size()-1; i != -1; i--) {
        update_visibility(points[i], polygon);
      }
    }
  }


  //add a new edge when vision ray passes a vertex
  void insert_edge(const Point_2& p,
                 std::vector<Pair>& edges,
                 const Point_2& dp) {

  }


  //remove edges that are not active any longer
  void remove_edges(std::vector<Halfedge_const_handle>& edges, const Ray_2& r) {
      typename std::vector<Halfedge_const_handle>::iterator eit = edges.begin();
      while (eit != edges.end()) {
          Point_2 p1 = (*eit)->target()->point();
          Point_2 p2 = (*eit)->source()->point();
          bool is_incident(false);
          if (is_on_ray(r, p1)) {
              is_incident = true;
          }
          else if (is_on_ray(r, p2)) {
              Point_2 tmp = p1;
              p1 = p2;
              p2 = tmp;
              is_incident = true;
          }
          if ( (is_incident && !CGAL::left_turn(r.source(), p1, p2)) || intersection_point(r, *eit) == r.source() )
          {
              eit = edges.erase(eit);
              continue;
          }

          else {
              eit++;
          }
      }
  }

  bool is_on_ray(const Ray_2& r, const Point_2& p) {
      return Direction_2(Vector_2(r.source(), p)) == Direction_2(r);
  }
  //return the type of the needle.
  //the vertices on the needle will be saved in collinear_vertices.
  Intersection_type needle(std::vector<Halfedge_const_handle>& edges, Ray_2& r, std::vector<Point_2>& collinear_vertices) {
      typename std::vector<Halfedge_const_handle>::iterator curr = edges.begin();
//        Point_2 p = r.source(), end1, end2;
      Vertex_const_handle vertex1;
      //flag shows whether the left side or right side of needle is blocked.
      bool block_left, block_right;
      do {
          Point_2 cross = intersection_point(r, *curr);
          if (cross != (*curr)->source()->point() && cross != (*curr)->target()->point()) {
              collinear_vertices.push_back(cross);
              return INNER;
          }
          if (CGAL::orientation(r.source(), (*curr)->source()->point(), (*curr)->target()->point()) == CGAL::COLLINEAR) {
            if (CGAL::compare_distance_to_point(r.source(), (*curr)->source()->point(), (*curr)->target()->point()) == CGAL::SMALLER)
              vertex1 = (*curr)->source();
            else
              vertex1 = (*curr)->target();
          }
          else {
              if (cross == (*curr)->source()->point()) {
                  vertex1 = (*curr)->source();
              }
              else {
                  vertex1 = (*curr)->target();
              }
          }
          if (collinear_vertices.empty() || vertex1->point() != collinear_vertices.back()) {
              collinear_vertices.push_back(vertex1->point());
              //flag shows whether the left side or right side of current vertex is blocked.
              //has_predecessor indicates whether this vertex is incident to an edge whose another end is between the source of ray and it,
              //because that will effect the value of block_left, block_right.
              bool left_v(false), right_v(false), has_predecessor(false);

              typename Arrangement_2::Halfedge_around_vertex_const_circulator first_edge, curr_edge;
              first_edge = curr_edge = vertex1->incident_halfedges();
              do {
                  switch (CGAL::orientation(r.source(), curr_edge->target()->point(), curr_edge->source()->point())) {
                  case CGAL::RIGHT_TURN :
                      right_v = true;
                      break;
                  case CGAL::LEFT_TURN :
                      left_v = true;
                      break;
                  case CGAL::COLLINEAR :
                      if (CGAL::compare_distance_to_point(r.source(), curr_edge->target()->point(), curr_edge->source()->point()) == CGAL::LARGER) {
                          has_predecessor = true;
                      }
                  }

              } while (++curr_edge != first_edge);
              if (has_predecessor) {
                  block_left = block_left || left_v;
                  block_right = block_right || right_v;
              }
              else {
                  block_left = left_v;
                  block_right = right_v;
              }
              if (block_left && block_right) {
                  return CORNER;
              }
          }
      } while (++curr != edges.end());
      return UNBOUNDED;
  }
  //debug
  void print_edges(std::vector<Halfedge_const_handle>& edges){
      for (int i = 0; i != edges.size(); i++) {
          Point_2 p1, p2;
          p1 = edges[i]->source()->point();
          p2 = edges[i]->target()->point();
          std::cout<<p1<<"->"<<p2<<std::endl;
      }
  }

  void print_vertex(const std::vector<Point_2>& polygon) {
    for (int i = 0; i != polygon.size(); i++) {
      std::cout<<polygon[i]<<std::endl;
    }
  }

  void build_arr(const std::vector<Point_2>& polygon, Arrangement_2& arr ) {
      for (int i = 0; i != polygon.size()-1; i++ ) {
          CGAL::insert(arr, Segment_2(polygon[i], polygon[i+1]));
      }
      //print_vectex(polygon);
      CGAL::insert(arr, Segment_2(polygon.front(), polygon.back()));
  }


  void conditional_regularize(Output_arrangement_2& out_arr, CGAL::Tag_true) {
    regularize_output(out_arr);
  }

  void conditional_regularize(Output_arrangement_2& out_arr, CGAL::Tag_false) {
    //do nothing
  }

  void regularize_output(Arrangement_2& out_arr) {
    typename Output_arrangement_2::Edge_iterator e_itr;
    for (e_itr = out_arr.edges_begin() ;
         e_itr != out_arr.edges_end() ; e_itr++) {
      Halfedge_handle he = e_itr;
      Halfedge_handle he_twin = he->twin();
      if (he->face() == he_twin->face()) {
        out_arr.remove_edge(he);
      }
    }
  }

};

} // end namespace CGAL



#endif


