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

#ifndef CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H
#define CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H

#include <iostream>
#include <vector>
#include <map>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Ray_2.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>

namespace CGAL {

//debug
template<typename Point_2>
void print(Point_2 p){
  std::cout<<p.x()<<","<< p.y()<<std::endl;
}

template<typename Point_2>
void print(std::vector<Point_2> ps){
    for (int i=0; i<ps.size(); i++)
    {
      print<Point_2>(ps[i]);
    }
}


template<typename Halfedge_handle>
void print_edges(std::vector<Halfedge_handle> es){
    for (int i=0; i<es.size(); i++)
    {
      std::cout << es->curve() << std::endl;
    }
}
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
  Rotational_sweep_visibility_2(const Input_arrangement_2 &arr): p_arr(&arr) {}


  Face_handle compute_visibility(const Point_2 &q, const Halfedge_const_handle e, Arrangement_2 &out_arr) {
      Arrangement_2 arrc = arr_in ; //copy of arr;
      Halfedge_handle ec; //copy of edge;
      for (Halfedge_handle eh = arrc.edges_begin(); eh != arrc.edges_end(); eh++) {
          if (eh->source()->point() == e-> source()->point() && eh->target()->point() == e->target()->point()) {
              ec = eh;
              break;
          }
      }
      //print_arrangement(arrc);
      //Insert a bounding box;
      Number_type xmin, xmax, ymin, ymax;
      Point_2 q1 = arrc.vertices_begin()->point();
      xmax = xmin = q1.x();
      ymin = ymax = q1.y();
      for (Vertex_const_handle vh = arrc.vertices_begin(); vh != arrc.vertices_end(); vh++) {
        if (vh->point().x() < xmin)
          xmin = vh->point().x();
        if (vh->point().x() > xmax)
          xmax = vh->point().x();
        if (vh->point().y() < ymin)
          ymin = vh->point().y();
        if (vh->point().y() > ymax)
          ymax = vh->point().y();
      }
      xmin -= 10;
      xmax += 10;
      ymin -= 10;
      ymax += 10;
      Point_2 p1(xmin, ymin), p2(xmax, ymin), p3(xmax, ymax), p4(xmin, ymax);
      CGAL::insert(arrc, Segment_2(p1, p2));
      CGAL::insert(arrc, Segment_2(p2, p3));
      CGAL::insert(arrc, Segment_2(p3, p4));
      CGAL::insert(arrc, Segment_2(p4, p1));

      if (ec->target()->point() == q) {
        Point_2 source = ec->source()->point();
        Point_2 target = ec->next()->target()->point();
        Halfedge_handle prev = ec->prev();
//          arrc.remove_edge(ec->next());
//          arrc.remove_edge(ec);

        bool check_remove;
        do {
          check_remove = false;
          for (Halfedge_handle eh = arrc.edges_begin(); eh != arrc.edges_end(); eh++) {
            if (eh->target()->point() == q || eh->source()->point() == q) {
                  arrc.remove_edge(eh);
                  check_remove = true;
                  break;
              }
          }
        } while (check_remove);

        std::vector<Point_2> polygon;
        visibility_region_impl(q, prev->face(), polygon);
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
                next_i ++;
            }
        }
        else {
            is_between = true;
            while (next_i != big) {
                if (CGAL::right_turn(source, q, polygon[next_i]) || CGAL::right_turn(q, target, polygon[next_i])) {
                    is_between = false;
                    break;
                }
                next_i ++;
            }
        }
        typename std::vector<Point_2>::iterator first = polygon.begin() + small;
        typename std::vector<Point_2>::iterator last = polygon.begin() + big;
        if (is_between) {
            std::vector<Point_2> polygon1(first, last+1);
            polygon1.push_back(q);
            build_arr(polygon1, out_arr);
        }
        else {
            std::vector<Point_2> polygon1(polygon.begin(), first+1);
            polygon1.push_back(q);
            for (int i = big; i != polygon.size(); i++) {
                polygon1.push_back(polygon[i]);
            }
            build_arr(polygon1, out_arr);
        }

      }
      else {
          Point_2 source = ec->source()->point();
          Point_2 target = ec->target()->point();
          Halfedge_handle eh1 = ec->next();
          arrc.remove_edge(ec);
          Face_handle fh = eh1->face();
          std::vector<Point_2> polygon;
          visibility_region_impl(q, fh, polygon);
          int source_i, target_i;
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
          while (CGAL::collinear(source, polygon[next_i], target))
              next_i++;
          typename std::vector<Point_2>::iterator first = polygon.begin() + small;
          typename std::vector<Point_2>::iterator last = polygon.begin() + big;
          if (CGAL::left_turn(source, target, polygon[next_i])) {
              std::vector<Point_2> polygon1(first, last+1);
              build_arr(polygon1, out_arr);
          }
          else {
              std::vector<Point_2> polygon1(polygon.begin(), first+1);
              for (int i = big; i != polygon.size(); i++) {
                  polygon1.push_back(polygon[i]);
              }
              build_arr(polygon1, out_arr);
          }

      }
      conditional_regularize(out_arr, Regularization_tag());

      if (out_arr.faces_begin()->is_unbounded())
          return ++out_arr.faces_begin();
      else
          return out_arr.faces_begin();

  }

  Face_handle compute_visibility(const Point_2 &q, const Face_const_handle f, Output_arrangement_2 &out_arr) {
    out_arr->clear();

    std::vector<Point_2> polygon;
    visibility_region_impl(q, f, polygon);
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

void attach(const Input_arrangement_2 &arr) {
  p_arr = &arr;
//  geom_traits = p_arr->geometry_traits();
}

void detach() {
  p_arr = NULL;
//  geom_traits = NULL;
  vertices.clear();
}

const Input_arrangement_2& arr() {
  return *p_arr;
}


private:
  //members
  typedef std::vector<Point_2> Vend;
  const Input_arrangement_2 *p_arr;
  bool            attach_tag;
  std::vector<Point_2> polygon;   //visibility polygon
  std::map<Point_2, Vend> vmap;   //vertex and two edges incident to it that might block vision


  //methods
  //compute visibility region between qa and qb
  void visibility_region_impl(const Point_2& q, const Point_2& a, const Point_2& b) {
      std::vector<Vertex_const_handle> vertices;                    //all vertices of the face.
      std::vector<Halfedge_const_handle> edges, active_edges;       //edges stores all halfedges of the face; and active_edges stores all halfedges that is currently intersected by the view ray.
      //preprocess the face
      input_face(fh, vertices, edges, q);

      //debug
//        for (int i = 0; i<vertices.size(); i++) {
//          print(vertices[i]->point());
//        }

      //initiation of vision ray
      Vector_2 dir;
      if (Direction_2(-1, 0) < Direction_2(Vector_2(q, (*vertices.rbegin())->point())))
      {
          dir = Vector_2(1, 0) + Vector_2(q, (*vertices.rbegin())->point());
      }
      else {
          dir = Vector_2(0, -1);
      }
      Ray_2 init_vision_ray(q, dir);
      //initiation of active_edges
      typename std::vector<Halfedge_const_handle>::iterator iter1;
      for (iter1 = edges.begin(); iter1 != edges.end(); iter1++)
      {
          insert_halfedge(active_edges, init_vision_ray, *iter1);
      }

      //angular sweep begins
      Ray_2 curr_vision_ray = init_vision_ray;
      typename std::vector<Vertex_const_handle>::iterator vit = vertices.begin(), begin_it, end_it;
      Halfedge_const_handle closest_edge;
      while (vit != vertices.end())
      {
          if (active_edges.empty())
          {
              begin_it = vit;
              end_it = vit;
              curr_vision_ray= Ray_2(q, (*begin_it)->point());
              Direction_2 d1(curr_vision_ray), d2;
              do {
                  d2 = Direction_2(Ray_2(q, (*end_it)->point()));
                  if (d1 != d2) break;
              } while (++end_it != vertices.end());
              add_edges(begin_it, end_it, active_edges, curr_vision_ray);

              //since active_edges is empty, that means adding new edges will bring in an unbounded edge to arrangement.
              Point_2 p1 = intersection_point(curr_vision_ray, active_edges[0]);
              polygon.push_back(p1);
//todo                CGAL::insert_curve(out_arr, Ray_2(p1, d1));
              closest_edge = active_edges[0];
              remove_edges(active_edges, curr_vision_ray);

              if (active_edges.empty()) {
                  //this means there is no edge that intersects curr_vision_ray
                  //except those collinear ones.
                  //because one unbounded ray has been inserted before,
                  //we don't need to do anything here.
              }
              else {
                  Point_2 p2 = intersection_point(curr_vision_ray, active_edges[0]);
                  update_visibility(p2, polygon);
              }
          }
          else {
              Point_2 right_p, left_p, mid_p;
              begin_it = vit;
              end_it = vit;
              curr_vision_ray= Ray_2(q, (*begin_it)->point());
              right_p = intersection_point(curr_vision_ray, active_edges[0]);
              Direction_2 d1(curr_vision_ray), d2;
              //find end_it such that all vertices between begin_it and end_it(not included) are collinear with query point.
              do {
                  d2 = Direction_2(Ray_2(q, (*end_it)->point()));
                  if (d1 != d2) break;
              } while (++end_it != vertices.end());
              add_edges(begin_it, end_it, active_edges, curr_vision_ray);
//                std::cout<<"after adding\n";
//                print_edges(active_edges);
              mid_p = intersection_point(curr_vision_ray, active_edges[0]);
              std::vector<Point_2> collinear_vertices;
              Intersection_type i_type = needle(active_edges, curr_vision_ray, collinear_vertices);
              switch (i_type) {
              case UNBOUNDED :
                  //todo:this part is not finished.
                  //remove right and collinear;
                  remove_edges(active_edges, curr_vision_ray);
                  update_visibility(right_p, polygon);
                  update_visibility(mid_p, polygon);
                  //todo CGAL::insert_curve();
                  if (!active_edges.empty()) {
                      left_p = intersection_point(curr_vision_ray, active_edges[0]);
                      update_visibility(left_p, polygon);
                  }
                  break;
              case CORNER :
                  //remove right and collinear;
                  remove_edges(active_edges, curr_vision_ray);

//                    std::cout<<"after removing\n";
//                    print_edges(active_edges);

                  left_p = intersection_point(curr_vision_ray, active_edges[0]);
                  update_visibility(right_p, polygon);
                  if (right_p == collinear_vertices[0]) {
                    insert_needle(collinear_vertices, polygon, true);
                  }
                  else if (left_p == collinear_vertices[0]) {
                    insert_needle(collinear_vertices, polygon, false);
                  }
                  update_visibility(left_p, polygon);
                  break;
              case INNER :
                  //remove right and collinear;
                  remove_edges(active_edges, curr_vision_ray);
                  if (collinear_vertices.size() < 2) {
                      //this means mid_p = left_p = right_p = furthest_p. no new vertex is found.
                  }
                  else {
                    //debug
//                      std::cout<<"print a needle:\n";
//                      print(collinear_vertices);
//                      std::cout<<"the left_p is "<<left_p.x()<<' '<<left_p.y()<<std::endl;
//                      std::cout<<"the right_p is "<<right_p.x()<<' '<<right_p.y()<<std::endl;
//                      std::cout<<"the front_p is "<<collinear_vertices[0].x()<<' '<<collinear_vertices[0].y()<<std::endl;
                    left_p = intersection_point(curr_vision_ray, active_edges[0]);
                    update_visibility(right_p, polygon);
                    if (right_p == collinear_vertices[0]) {
                      insert_needle(collinear_vertices, polygon, true);
                    }
                    else if (left_p == collinear_vertices[0]) {
                      insert_needle(collinear_vertices, polygon, false);
                    }
                    update_visibility(left_p, polygon);
                  }
                  break;
              }
          }
          vit = end_it;
      }

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

  Point_2 intersection_point(Ray_2 ray, Halfedge_const_handle seg) {
      return intersection_point(ray, halfedge2seg(seg));
  }

  //convertor for halfedge to segment
  Segment_2 halfedge2seg(Halfedge_const_handle e){
      return Segment_2(e->source()->point(), e->target()->point());
  }


  //given two edges incident to a vision ray at the same point, find which one is first seen in sweeping.
  bool is_closer(const Ray_2 &ray, Halfedge_const_handle seg1, Halfedge_const_handle seg2) {
      Point_2 shared = intersection_point(ray, seg1);
      Point_2 end1, end2;
      if (shared == seg1->source()->point())
          end1 = seg1->target()->point();
      else
          end1 = seg1->source()->point();

      if (shared == seg2->source()->point())
          end2 = seg2->target()->point();
      else
          end2 = seg2->source()->point();
      if (CGAL::right_turn(ray.source(), shared, end1) && !CGAL::right_turn(ray.source(), shared, end2))
          return true;
      if (CGAL::right_turn(ray.source(), shared, end2) && !CGAL::right_turn(ray.source(), shared, end1))
          return false;
      switch (CGAL::orientation(ray.source(), shared, end1)) {
      case CGAL::COLLINEAR:
          return (CGAL::left_turn(ray.source(), shared, end2));
      case CGAL::RIGHT_TURN:
          return (CGAL::right_turn(end1, shared, end2));
      case CGAL::LEFT_TURN:
          return (CGAL::left_turn(end1, shared, end2));
      }
  }

  //insert newly-discovered edges into active_edges according to its intersection with the view ray.
  void insert_halfedge(std::vector<Halfedge_const_handle> &active_edges, const Ray_2 &ray, Halfedge_const_handle edge)
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

  //insert vh into vertices by the angle between ray<p, vh> and positive x-ray.
  //if the angles are the same, compare their distances to p.
  void sort_vertex(std::vector<Vertex_const_handle>& vertices, Vertex_const_handle vh, const Point_2& p)
  {
      typename std::vector<Vertex_const_handle>::iterator first = vertices.begin();
      Vector_2 vector_of_v(p, vh->point());
      Direction_2 dir_of_v(vector_of_v);
      while (first != vertices.end())
      {
          if (vh->point() == (*first)->point()) {
              return;
          }
          Vector_2 vector1(p, (*first)->point());
          Direction_2 d1(vector1);
          if (dir_of_v < d1)
              break;
          //if same angles are the same, put the vertex which is closer to query point before.
          if (dir_of_v == d1 && CGAL::compare_distance_to_point(p, vh->point(), (*first)->point()) == CGAL::SMALLER)
              break;
          ++first;
      }
      vertices.insert(first, vh);
  }


  //traverse the face to get all edges and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh,
                   std::vector<Vertex_const_handle>& vertices,
                   std::vector<Halfedge_const_handle>& edges,
                   const Point_2& p)
  {
      typename Arrangement_2::Ccb_halfedge_const_circulator curr = fh->outer_ccb();
      typename Arrangement_2::Ccb_halfedge_const_circulator circ = curr;
      do {
          sort_vertex(vertices, curr->source(), p);
          edges.push_back(curr);
      } while (++curr != circ);
      typename Arrangement_2::Hole_const_iterator hi;
      for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
          typename Arrangement_2::Ccb_halfedge_const_circulator c1 = *hi, c2 = *hi;
          do {
              sort_vertex(vertices, c1->source(), p);
              edges.push_back(c1);
          } while (++c1 != c2);
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
  void add_edge(Vertex_const_handle vh,
                 std::vector<Halfedge_const_handle>& edges,
                 const Ray_2& r) {
      typename Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
      first = curr = vh->incident_halfedges();
      do {
          if (!CGAL::right_turn(r.source(), vh->point(), curr->source()->point()))
          {
              insert_halfedge(edges, r, curr);
          }
      } while (++curr != first);

  }
  //add new edges
  void add_edges(typename std::vector<Vertex_const_handle>::iterator begin_it,
                 typename std::vector<Vertex_const_handle>::iterator end_it,
                 std::vector<Halfedge_const_handle>& edges,
                 const Ray_2& r)
  {
      do {
          add_edge(*begin_it, edges, r);
      } while (++begin_it != end_it);

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




  //angular sweep a vertice of face.
  void sweep_vertex(std::vector<Halfedge_const_handle> &active_edges, const Point_2 &query, Vertex_const_handle vh, std::vector<Point_2> &polygon )
  {
      //closest_edge_copy is a copy of the closest edge to query point in active_edges before sweeping.
      Halfedge_const_handle closest_edge_copy = active_edges[0];
      Ray_2 ray(query, vh->point());
      int add_count(0);
      int del_count(0);

      //delete all edges in active_edges which is incident to v, because they has been sweeped over
      typename std::vector<Halfedge_const_handle>::iterator edge_iter = active_edges.begin();
      while (edge_iter != active_edges.end()) {
          if (((*edge_iter)->source()->point() == vh->point()) || ((*edge_iter)->target()->point() == vh->point()))
          {
              edge_iter = active_edges.erase(edge_iter);
              ++del_count;
          }
          else {
              ++edge_iter;
          }
      }

      //add all edges which is incident to v but not in active_edges before to active_edges
      typename Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
      first = curr = vh->incident_halfedges();
      do {
          if (CGAL::left_turn(query, vh->point(), curr->source()->point()))
          {
              insert_halfedge(active_edges, ray, curr);
              ++add_count;
          }
          else if (CGAL::collinear(query, vh->point(), curr->source()->point()) &&
                  CGAL::compare_distance_to_point(query, vh->point(), curr->source()->point()) == CGAL::SMALLER)
          {
              insert_halfedge(active_edges, ray, curr);
              ++add_count;
          }
      } while (++curr != first);

      //update the visibility region
      if (closest_edge_copy != active_edges[0]) {
          //when the closest edge changed
          if (del_count > 0 && add_count > 0) {
              //some edges are added and some are deleted, which means the vertice sweeped is a vertice of visibility polygon.
              update_visibility(vh->point(), polygon);
          }
          if (del_count == 0 && add_count > 0) {
              //only add some edges, means the view ray is blocked by new edges.
              //therefore first add the intersection of view ray and former closet edge, then add the vertice sweeped.
              update_visibility(intersection_point(ray, halfedge2seg(closest_edge_copy)), polygon);
              update_visibility(vh->point(), polygon);
          }
          if (del_count > 0 && add_count == 0) {
              //only delete some edges, means some block is moved and the view ray can reach the segments after the block.
              update_visibility(vh->point(), polygon);
              update_visibility(intersection_point(ray, halfedge2seg(active_edges[0])), polygon);
          }
      }

  }

  void conditional_regularize(Output_arrangement_2 &out_arr, CGAL::Tag_true) {
    regularize_output(out_arr);
  }

  void conditional_regularize(Output_arrangement_2 &out_arr, CGAL::Tag_false) {
    //do nothing
  }

  void regularize_output(Arrangement_2 &out_arr) {
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


//For debug. Print all edges of arrangements into console.
template <typename Arrangement_2>
void print_arrangement(const Arrangement_2 &arr) {
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  Edge_const_iterator eit;
  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
    std::cout << "[" << eit->curve() << "]" << std::endl;
}

} // end namespace CGAL



#endif


