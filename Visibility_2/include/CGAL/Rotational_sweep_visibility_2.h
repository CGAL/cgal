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
    //debug
//    print_vertex(polygon);

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
    out_arr.clear();
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
  typedef std::vector<Point_2>          Pvec;
  typedef std::pair<Point_2, Point_2>   Pair;

  Input_arrangement_2 *p_arr;
  Point_2         q;
  Point_2         dp;
  Pvec polygon;   //visibility polygon
  std::map<Point_2, Pvec> vmap;   //vertex and two edges incident to it that might block vision
  std::map<Pair, int> edx;   //index of edge in the heap
  std::vector<Pair>  heap;
  std::vector<Point_2> vs;          //angular sorted vertices
  bool is_vertex_query;
  bool is_edge_query;

  template <class NT>
  int quadrant(NT x, NT y) {
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

  bool do_intersect(const Point_2& q,
                    const Point_2& dp,
                    const Point_2& p1,
                    const Point_2& p2) {
    if (CGAL::collinear(q, dp, p1))
      return (quadrant(p1.x()-q.x(), p1.y()-q.y()) == quadrant(dp.x()-q.x(), dp.y()-q.y()));

    if (CGAL::collinear(q, dp, p2))
      return (quadrant(p2.x()-q.x(), p2.y()-q.y()) == quadrant(dp.x()-q.x(), dp.y()-q.y()));

    return (CGAL::orientation(q, dp, p1) != CGAL::orientation(q, dp, p2) && CGAL::orientation(q, p1, dp) == CGAL::orientation(q, p1, p2));

  }

  void funnel(int i, int j) {
    Pvec right, left;
    bool block_left(false), block_right(false);
    Point_2 former = vs[i];
    for (int l=i; l<j; l++) {
      bool left_v(false), right_v(false), has_predecessor(false);
      for (int k=0; k<vmap[vs[l]].size(); k++) {
        Point_2 temp= vmap[vs[l]][k];
        if (temp == former)
          has_predecessor = true;
        if (CGAL::left_turn(q, vs[l], temp))
          left_v = true;
        else
          right_v = CGAL::right_turn(q, vs[l], temp);
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


  void visibility_region_impl(const Face_const_handle f, const Point_2& q, const Point_2& a, const Point_2& b) {
    vs.clear();
    polygon.clear();
    heap.clear();
    vmap.clear();
    edx.clear();

    std::vector<Pair> bbox;
    input_face(f, q, a, b, bbox);

    //debug
    std::cout<<"new query point: " <<q<<std::endl;
    print_vertex(vs);

    //initiation of vision ray
    Vector_2 dir;
    if (Direction_2(-1, 0) < Direction_2(Vector_2(q, vs.back())))
    {
        dir = Vector_2(1, 0) + Vector_2(q, vs.back());
    }
    else {
        dir = Vector_2(0, -1);
    }

    dp = Point_2(q.x()+dir.x(), q.y()+dir.y());

    //initiation of active_edges

    Ccb_halfedge_const_circulator curr = f->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      Point_2 p1 = curr->target()->point();
      Point_2 p2 = curr->source()->point();
      if (q != p1 && q != p2 && do_intersect(q, dp, p1, p2))
        heap_insert(create_pair(p1, p2));
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = f->holes_begin(); hi != f->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator c1 = *hi, c2 = *hi;
      do {
        Point_2 p1 = c1->target()->point();
        Point_2 p2 = c1->source()->point();
        if (q != p1 && q != p2 && do_intersect(q, dp, p1, p2))
          heap_insert(create_pair(p1, p2));
      } while (++c1 != c2);
    }
    for (int i=0; i!=bbox.size(); i++) {
      if (do_intersect(q, dp, bbox[i].first, bbox[i].second))
        heap_insert(bbox[i]);
    }
    std::cout<<"blew is initial active edges.\n";
    print_edges(heap);
    //angular sweep begins
    for (int i=0; i!=vs.size(); i++) {
      dp = vs[i];
      Point_2 v = dp;
      Pair ce = heap.front(); //save closest edge;
      int insert_cnt(0), remove_cnt(0);
      for (int j=0; j!=vmap[v].size(); j++) {
        Pair e = create_pair(v, vmap[v][j]);
        if (edx.count(e)) {
          heap_remove(edx[e]);
          remove_cnt++;
        }
        else {
          heap_insert(e);
          insert_cnt++;
        }
      }
      if (ce != heap.front()) {
        //when the closest edge changed
        if (remove_cnt > 0 && insert_cnt > 0) {
            //some edges are added and some are deleted, which means the vertice sweeped is a vertice of visibility polygon.
            update_visibility(v);
        }
        if (remove_cnt == 0 && insert_cnt > 0) {
            //only add some edges, means the view ray is blocked by new edges.
            //therefore first add the intersection of view ray and former closet edge, then add the vertice sweeped.
          update_visibility(ray_seg_intersection(q, dp, ce.first, ce.second));
          update_visibility(v);
        }
        if (remove_cnt > 0 && insert_cnt == 0) {
          //only delete some edges, means some block is moved and the view ray can reach the segments after the block.
          update_visibility(v);
          update_visibility(ray_seg_intersection(q, dp, heap.front().first, heap.front().second));
        }
      }
      //debug
      //print_edges(heap);
    }
    //print_vertex(polygon);
  }

  Pair create_pair(const Point_2& p1, const Point_2& p2){
    assert(p1 != p2);
    if (p1 < p2)
      return Pair(p1, p2);
    else
      return Pair(p2, p1);
  }

//todo add edge location record
  void heap_insert(const Pair& e) {
    heap.push_back(e);
    int i = heap.size()-1;
    edx[e] = i;
    int parent = (i-1)/2;
    while (i!=0 && is_closer(q, dp, heap[i], heap[parent])){
      heap_swap(i, parent);
      i = parent;
      parent = (i-1)/2;
    }
  }

  void heap_remove(int i) {
    edx.erase(heap[i]);
    if (i== heap.size()-1)
    {
      heap.pop_back();
    }
    else {
      heap[i] = heap.back();
      edx[heap[i]] = i;
      heap.pop_back();
      bool swapped;
      do {
        int left_son = i*2+1;
        int right_son = i*2+2;
        int closest = i;
        if (left_son < heap.size() && is_closer(q, dp, heap[left_son], heap[i])) {
          closest = left_son;
        }
        if (right_son < heap.size() && is_closer(q, dp, heap[right_son], heap[closest])) {
          closest = right_son;
        }
        swapped = false;
        if (closest != i) {
          heap_swap(i, closest);
          i = closest;
          swapped = true;
        }
      } while(swapped);
    }
  }

  void heap_swap(int i, int j) {
    edx[heap[i]] = j;
    edx[heap[j]] = i;
    Pair temp = heap[i];
    heap[i] = heap[j];
    heap[j] = temp;
  }

  //todo logically it's wrong.
  bool is_closer(const Point_2& q, const Point_2& dp, const Pair& e1, const Pair& e2) {
    Point_2 p1 = ray_seg_intersection(q, dp, e1.first, e1.second);
    Point_2 p2 = ray_seg_intersection(q, dp, e2.first, e2.second);
    if (p1 == p2) {
      Point_2 end1, end2;
      if (p1 == e1.first)
        end1 = e1.second;
      else
        end1 = e1.first;
      if (p2 == e2.first)
        end2 = e2.second;
      else
        end2 = e2.first;
      if (CGAL::left_turn(q, p1, end1))
        return CGAL::left_turn(end1, p1, end2);
      else
        return false;
    }
    else {
      return CGAL::compare_distance_to_point(q, p1, p2)==CGAL::SMALLER;
    }
  }

  Point_2 ray_seg_intersection(
      const Point_2& q, const Point_2& dp, // the ray
      const Point_2& s, const Point_2& t) // the segment
  {
    Ray_2 ray(q,dp);
    Segment_2 seg(s,t);
    CGAL::Object result = CGAL::intersection(ray, seg);
    if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        return *ipoint;
    }
    else {
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
        std::cout<<"doesn't intersect\n";
        std::cout<<q<<','<<dp<<"   "<<s<<','<<t<<std::endl;
        assert(false);
      }
    }
  }


//  {
//    Ray_2 ray(q,dp);
//    Segment_2 seg(s,t);
//    if (typename K::Do_intersect_2()(ray,seg)) {
//      std::cout<<"doesn't intersect\n";
//      std::cout<<q<<','<<dp<<"   "<<s<<','<<t<<std::endl;
//      assert(false);
//    }
//    CGAL::Object obj = typename K::Intersect_2()(ray,seg);
//    Point_2 result =  object_cast<Point_2>(obj);
//    return result;
//  }

  void update_visibility(const Point_2& p){
    if (polygon.empty())
      polygon.push_back(p);
    else
    {
      if (polygon.back() != p){
        polygon.push_back(p);
      }
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

  bool is_good_edge(const Point_2& q,
                    const Point_2& a,
                    const Point_2& b,
                    const Point_2& v,
                    const Point_2& n) {
    if (is_vertex_query)
      return !(n==q);
    if (is_edge_query)
      return !((v==a && n==b)||(v==b && n==a));
  }

  //traverse the face to get all edges and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh,
                   const Point_2& q,
                   const Point_2& a,
                   const Point_2& b,
                   std::vector<Pair>& bbox)
  {

    Ccb_halfedge_const_circulator curr = fh->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      Point_2 v = curr->target()->point();
      if (v == q)
        continue;
      vs.push_back(v);
      //debug
//      std::cout<<"after one push back"<<std::endl;
//      print_vertex(vs);
      Point_2 nei[2] = {curr->source()->point(), curr->next()->target()->point()};
      Pvec neighbor;
      for (int i=0; i<2; i++){
        if ((!is_vertex_query && !is_edge_query) || is_good_edge(q,a,b,v,nei[i]))
          neighbor.push_back(nei[i]);
      }
      vmap[v] = neighbor;
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator c1 = *hi, c2 = *hi;
      do {
        Point_2 v = c1->target()->point();
        if (v == q)
          continue;
        vs.push_back(v);
        Point_2 nei[2] = {c1->source()->point(), c1->next()->target()->point()};
        Pvec neighbor;
        for (int i=0; i<2; i++){
          if ((!is_vertex_query && !is_edge_query) || is_good_edge(q,a,b,v,nei[i]))
            neighbor.push_back(nei[i]);
        }
        vmap[v] = neighbor;
      } while (++c1 != c2);
    }

    if (q != a) {
      Number_type xmin, xmax, ymin, ymax;
      Point_2 q1 = vs.front();
      xmax = xmin = q1.x();
      ymin = ymax = q1.y();
      for (typename Pvec::iterator it= vs.begin(); it!=vs.end(); it++) {
        Point_2 q1 = *it;
        if (q1.x() < xmin)    xmin = q1.x();
        if (q1.x() > xmax)    xmax = q1.x();
        if (q1.y() < ymin)    ymin = q1.y();
        if (q1.y() > ymax)    ymax = q1.y();
      }
      xmin -= 10;
      xmax += 10;
      ymin -= 10;
      ymax += 10;
      Point_2 box[4] = {Point_2(xmin, ymin), Point_2(xmax, ymin),
                        Point_2(xmax, ymax), Point_2(xmin, ymax)};
      for (int i=0; i<4; i++) {
        vs.push_back(box[i]);
        Pvec pvec;
        pvec.push_back(box[(i+3)%4]);
        pvec.push_back(box[(i+1)%4]);
        vmap[box[i]] = pvec;

        bbox.push_back(create_pair(box[i], box[(i+1)%4]));
      }
    }
    //debug
    std::cout<<"before quick sort, vs size: "<<vs.size()<<std::endl;

    quick_sort(vs, 0, vs.size()-1);

   //debug print_vertex(vs);

    for (int i=0; i!=vs.size(); i++) {
      int j = i+1;
      while (j != vs.size()) {
        if (!CGAL::collinear(q, vs[i], vs[j]))
          break;
        j++;
      }
      if (j-i>1)
        funnel(i, j);
      i = j-1;
    }
  }

  void qs_swap(Pvec& vs, int i, int j) {
    Point_2 temp = vs[i];
    vs[i] = vs[j];
    vs[j] = temp;
  }

  int partition(Pvec& vs, int left, int right, int pivotIndex) {
    Point_2 pivot_p = vs[pivotIndex];
    qs_swap(vs, pivotIndex, right);
    int storeIndex = left;
    for (int i=left; i<right; i++) {
      if (compare_angle(vs[i], pivot_p)) {
        qs_swap(vs, i, storeIndex);
        storeIndex += 1;
      }
    }
    qs_swap(vs, storeIndex, right);
    return storeIndex;
  }

  void quick_sort(Pvec& vs, int left, int right) {
    if (left < right) {
      int pivotIndex = left;
      int pivotNewIndex = partition(vs, left, right, pivotIndex);
      quick_sort(vs, left, pivotNewIndex-1);
      quick_sort(vs, pivotNewIndex+1, right);
    }
  }

  bool is_on_ray(const Ray_2& r, const Point_2& p) {
      return Direction_2(Vector_2(r.source(), p)) == Direction_2(r);
  }


  //debug
  void print_edges(std::vector<Pair>& edges){
    std::cout<<edges.size()<<" edges in the heap now.\n";
    for (int i = 0; i != edges.size(); i++) {
          std::cout<<edges[i].first<<"->"<<edges[i].second<<std::endl;
      }
  }

  void print_vertex(const Pvec& polygon) {
    std::cout<<"print points in vector\n";
    for (int i = 0; i != polygon.size(); i++) {
      std::cout<<polygon[i]<<std::endl;
    }
  }
  void print_edx() {
    typename std::map<Pair, int>::iterator map_it = edx.begin();
    std::cout<<"print edx\n";
    while (map_it != edx.end()) {
      std::cout<<map_it->first.first<<"->"<<map_it->first.second<<":"<<map_it->second<<std::endl;
      map_it++;
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


