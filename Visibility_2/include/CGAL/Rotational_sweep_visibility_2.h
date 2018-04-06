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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s):  Kan Huang <huangkandiy@gmail.com>
//

#ifndef CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H
#define CGAL_ROTATIONAL_SWEEP_VISIBILITY_2_H

#include <CGAL/license/Visibility_2.h>


#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <boost/unordered_map.hpp> 
#include <iterator>


namespace CGAL {

template<class Arrangement_2_ , class RegularizationCategory = CGAL::Tag_true >
class Rotational_sweep_visibility_2 {
public:
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::
    Ccb_halfedge_const_circulator                 Ccb_halfedge_const_circulator;
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

  typedef RegularizationCategory              Regularization_category;
  typedef CGAL::Tag_true                      Supports_general_polygon_category;
  typedef CGAL::Tag_true                      Supports_simple_polygon_category;

private:
  typedef std::vector<Point_2>                          Points;
  typedef Vertex_const_handle                           VH;
  typedef std::vector<VH>                               VHs;
  typedef Halfedge_const_handle                         EH;
  typedef std::vector<EH>                               EHs;

  class Less_edge: public CGAL::binary_function<EH, EH, bool> {
    const Geometry_traits_2* geom_traits;
  public:
    Less_edge() {}
    Less_edge(const Geometry_traits_2* traits):geom_traits(traits) {}
    bool operator() (const EH e1, const EH e2) const {
      if (e1 == e2)
        return false;
      else {
         return &(*e1)<&(*e2);
     }
// if (e1->source() == e2->source())
//   return Visibility_2::compare_xy_2(geom_traits,
//      e1->target()->point(), e2->target()->point()) == SMALLER;
// else
//   return Visibility_2::compare_xy_2(geom_traits,
//      e1->source()->point(), e2->source()->point()) == SMALLER;

    }
  };

  class Less_vertex: public CGAL::binary_function<VH, VH, bool> {
    const Geometry_traits_2* geom_traits;
  public:
    Less_vertex() {}
    Less_vertex(const Geometry_traits_2* traits):geom_traits(traits) {}
    bool operator() (const VH v1, const VH v2) const {
      if (v1 == v2)
        return false;
      else
        // I know this is dirty but it speeds up by 25%. Michael 
        return &(*v1)<&(*v2);
//        return Visibility_2::
//          compare_xy_2(geom_traits, v1->point(), v2->point()) == SMALLER;
    }
  };

  class Closer_edge: public CGAL::binary_function<EH, EH, bool> {
    const Geometry_traits_2* geom_traits;
    Point_2 q;
  public:
    Closer_edge() {}
    Closer_edge(const Geometry_traits_2* traits, const Point_2& q) :
        geom_traits(traits), q(q) {}

    int vtype(const Point_2& c, const Point_2& p) const {
      switch(Visibility_2::orientation_2(geom_traits, q, c, p)) {
      case COLLINEAR:
        if (Visibility_2::less_distance_to_point_2(geom_traits, q, c, p))
          return 0;
        else
          return 3;
      case RIGHT_TURN:
        return 1;
      case LEFT_TURN:
        return 2;
      default: CGAL_assume(false);
      }
      return -1;
    }

    bool operator() (const EH& e1, const EH& e2) const {
      if (e1 == e2)
        return false;
      const Point_2& s1=e1->source()->point(),
                     t1=e1->target()->point(),
                     s2=e2->source()->point(),
                     t2=e2->target()->point();
      if (e1->source() == e2->source()) {

        int vt1 = vtype(s1, t1),
            vt2 = vtype(s1, t2);
        if (vt1 != vt2)
          return vt1 > vt2;
        else
          return (Visibility_2::orientation_2(geom_traits, s1, t2, t1)==
                  Visibility_2::orientation_2(geom_traits, s1, t2, q));
      }

      if (e1->target() == e2->source()) {
//          const Point_2& p1 = s1,
//                   p2 = t2,
//                   c = s2;
        int vt1 = vtype(t1, s1),
            vt2 = vtype(t1, t2);
        if (vt1 != vt2)
          return vt1 > vt2;
        else
          return (Visibility_2::orientation_2(geom_traits, s2, t2, s1)==
                  Visibility_2::orientation_2(geom_traits, s2, t2, q));
      }

      if (e1->source() == e2->target()) {
//            const Point_2& p1 = t1,
//                     p2 = s2,
//                     c = s1;
        int vt1 = vtype(s1, t1),
            vt2 = vtype(s1, s2);
        if (vt1 != vt2)
          return vt1 > vt2;
        else return (Visibility_2::orientation_2(geom_traits, s1, s2, t1)==
                     Visibility_2::orientation_2(geom_traits, s1, s2, q));
      }

      if (e1->target() == e2->target()) {
//              const Point_2& p1 = s1,
//                       p2 = s2,
//                       c = t1;
        int vt1 = vtype(t1, s1),
            vt2 = vtype(t1, s2);
        if (vt1 != vt2)
          return vt1 > vt2;
        else return (Visibility_2::orientation_2(geom_traits, t1, s2, s1)==
                     Visibility_2::orientation_2(geom_traits, t1, s2, q));
      }

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
          //q is collinear with e1 not with e2.
          if (Visibility_2::collinear(geom_traits, s2, t2, s1))
            return (Visibility_2::orientation_2(geom_traits, s2, t2, q)
                    == Visibility_2::orientation_2(geom_traits, s2, t2, t1));
          else
            return (Visibility_2::orientation_2(geom_traits, s2, t2, q)
                    == Visibility_2::orientation_2(geom_traits, s2, t2, s1));
        }
        break;
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
          if(Visibility_2::orientation_2(geom_traits, s1, t1, t2) == RIGHT_TURN)
            return Visibility_2::orientation_2(geom_traits, s2, t2, q)
                == Visibility_2::orientation_2(geom_traits, s2, t2, s1);
          else
            return true;
        default: CGAL_assume(false);
        }
        break;
      case LEFT_TURN:
        switch (Visibility_2::orientation_2(geom_traits, s1, t1, s2)) {
        case COLLINEAR:
          return Visibility_2::orientation_2(geom_traits, s1, t1, t2)!=e1q;
        case LEFT_TURN:
          if(Visibility_2::orientation_2(geom_traits, s1, t1, t2) == RIGHT_TURN)
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
        default: CGAL_assume(false);
        }
      }

      CGAL_assume(false);
      return false;
    }

  };
  
  const Arrangement_2 *p_arr;
  const Geometry_traits_2 *geom_traits;


  mutable Point_2 q;                               //query point
  mutable Points polygon;                          //visibility polygon

  mutable std::map<VH, EHs, Less_vertex> incident_edges;

  mutable std::map<EH, int, Less_edge> edx;        //index of active edges in
                                                   //the heap

  mutable std::set<EH, Closer_edge> active_edges;  //a set of edges that
                                                   //intersect the current
                                                   //vision ray.

  mutable VHs vs;                           //angular sorted vertices
  mutable EHs bad_edges;                    //edges that pass the query point
  mutable VH  cone_end1;                    //an end of visibility cone
  mutable VH  cone_end2;                    //another end of visibility cone

  mutable typename Points::size_type cone_end1_idx;
                                            //index of cone_end1->point() in
                                            //visibility polygon

  mutable typename Points::size_type cone_end2_idx;
                                            //index of cone_end2->point() in
                                            //visibility polygon

  mutable bool is_vertex_query;
  mutable bool is_edge_query;
  mutable bool is_face_query;
  mutable bool is_big_cone;               //whether the angle of
                                          //visibility_cone is greater than pi.

public:
  Rotational_sweep_visibility_2(): p_arr(NULL), geom_traits(NULL) {}
  Rotational_sweep_visibility_2(const Arrangement_2& arr): p_arr(&arr) {
    geom_traits = p_arr->geometry_traits();
  }

  const std::string name() const { return std::string("R_visibility_2"); }
  
  template <typename VARR> 
  typename VARR::Face_handle 
  compute_visibility(
          const Point_2& q, const Halfedge_const_handle e, VARR& arr_out) const
  {
    arr_out.clear();
    bad_edges.clear();
    this->q = q;

    if (Visibility_2::compare_xy_2(geom_traits, q, e->target()->point())==EQUAL)
    {
      is_vertex_query = true;
      is_edge_query = false;
      is_face_query = false;
      cone_end1 = e->source();
      cone_end2 = e->next()->target();
      is_big_cone = CGAL::right_turn(cone_end1->point(), q, cone_end2->point());

      typename Arrangement_2::
                Halfedge_around_vertex_const_circulator first, curr;
      first = curr = e->target()->incident_halfedges();
      do {
        if (curr->face() == e->face())
          bad_edges.push_back(curr);
        else if (curr->twin()->face() == e->face())
          bad_edges.push_back(curr->twin());
      } while (++curr != first);
    }
    else {
      is_vertex_query = false;
      is_edge_query = true;
      is_face_query = false;
      cone_end1 = e->source();
      cone_end2 = e->target();
      bad_edges.push_back(e);
      is_big_cone = false;
    }
    visibility_region_impl(e->face(), q);

    //decide which inside of the visibility butterfly is needed.
    typename Points::size_type small_idx, big_idx;
    if ( cone_end1_idx < cone_end2_idx ) {
      small_idx = cone_end1_idx;
      big_idx = cone_end2_idx;
    }
    else {
      small_idx = cone_end2_idx;
      big_idx = cone_end1_idx;
    }
    typename Points::size_type next_idx = small_idx + 1;
    bool is_between;
    //indicate whether the shape between small_idx and big_idx is the visibility
    //region required.
    if (CGAL::right_turn(cone_end1->point(), q, cone_end2->point())) {
      is_between = false;
      while (next_idx != big_idx) {
        if (CGAL::left_turn(cone_end1->point(), q, polygon[next_idx]) ||
            CGAL::left_turn(q, cone_end2->point(), polygon[next_idx]))
        {
          is_between = true;
          break;
        }
        next_idx++;
      }
    }
    else {
      is_between = true;
      while (next_idx != big_idx) {
        if (CGAL::right_turn(cone_end1->point(), q, polygon[next_idx]) ||
            CGAL::right_turn(q, cone_end2->point(), polygon[next_idx]))
        {
          is_between = false;
          break;
        }
        next_idx++;
      }
    }

    typename Points::iterator first = polygon.begin();
    std::advance(first, small_idx);
    typename Points::iterator last = polygon.begin();
    std::advance(last, big_idx);

    if (is_between) {
      Points polygon_out(first, last+1);
      if (is_vertex_query)
        polygon_out.push_back(q);
      Visibility_2::report_while_handling_needles<Rotational_sweep_visibility_2>
        (geom_traits, q, polygon_out, arr_out);
    }
    else {
      Points polygon_out(polygon.begin(), first+1);
      if (is_vertex_query) polygon_out.push_back(q);
      for (typename Points::size_type i = big_idx; i != polygon.size(); i++) {
        polygon_out.push_back(polygon[i]);
      }
      Visibility_2::report_while_handling_needles<Rotational_sweep_visibility_2>
        (geom_traits, q, polygon_out, arr_out);
    }

    Visibility_2::conditional_regularize(arr_out, Regularization_category());

    if (arr_out.faces_begin()->is_unbounded())
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }

  template <typename VARR> 
  typename VARR::Face_handle 
  compute_visibility(
          const Point_2& q, const Face_const_handle f, VARR& arr_out) const
  {
    arr_out.clear();
    this->q = q;
    is_vertex_query = false;
    is_edge_query = false;
    is_face_query = true;

    visibility_region_impl(f, q);

    Visibility_2::report_while_handling_needles<Rotational_sweep_visibility_2>
      (geom_traits, q, polygon, arr_out);

    Visibility_2::conditional_regularize(arr_out, Regularization_category());

    if (arr_out.faces_begin()->is_unbounded())
      return ++arr_out.faces_begin();
    else
      return arr_out.faces_begin();
  }

bool is_attached() const {
  return (p_arr != NULL);
}

void attach(const Arrangement_2& arr) {
  p_arr = &arr;
  geom_traits = p_arr->geometry_traits();
}

void detach() {
  p_arr = NULL;
  geom_traits = NULL;
}

const Arrangement_2& arrangement_2() const {
  return *p_arr;
}

private:
  //get the neighbor of v along edge e
  VH get_neighbor(const EH e, const VH v) const {
    if (e->source() == v)
      return e->target();
    else
      return e->source();
  }

  //check whether ray(q->dp) intersects segment(p1, p2)
  bool do_intersect_ray(const Point_2& q,
                        const Point_2& dp,
                        const Point_2& p1,
                        const Point_2& p2) const
  {
    return (CGAL::orientation(q, dp, p1) != CGAL::orientation(q, dp, p2) &&
            CGAL::orientation(q, p1, dp) == CGAL::orientation(q, p1, p2));
  }

  //arrange vertices that on a same vision ray in a 'funnel' order
  void funnel(typename VHs::size_type i, typename VHs::size_type j) const {
    VHs right, left;
    //whether the edges incident to a vertex block the left side and right side
    //of current vision ray.
    bool block_left(false), block_right(false);
    VH former = vs[i], nb;
    for (typename VHs::size_type l=i; l<j; l++) {
      bool left_v(false), right_v(false), has_predecessor(false);
        EHs& edges = incident_edges[vs[l]];
        for (typename EHs::size_type k=0; k<edges.size(); k++) {
          nb = get_neighbor(edges[k], vs[l]);
        if ( nb == former )  {
          has_predecessor = true;
          break;
        }
        if (CGAL::left_turn(q, vs[l]->point(), nb->point()))
          left_v = true;
        else
          right_v = CGAL::right_turn(q, vs[l]->point(), nb->point());
      }
      if (has_predecessor) {
        //if the current vertex connects to the vertex before by an edge,
        //the vertex before can help it to block.
        block_left = block_left || left_v;
        block_right = block_right || right_v;
      }
      else {
        block_left = left_v;
        block_right = right_v;
      }
      if (block_left && block_right) {
        //when both sides are blocked,
        //there is no need to change the vertex after.
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
    for (typename VHs::size_type l=0; l!=right.size(); l++)
      vs[i+l] = right[l];
    for (typename VHs::size_type l=0; l!=left.size(); l++)
      vs[i+l+right.size()] = left[left.size()-1-l];
  }



  void visibility_region_impl(const Face_const_handle f, const Point_2& q) const
  {
    vs.clear();
    polygon.clear();
    active_edges = std::set<EH, Closer_edge>(Closer_edge(geom_traits, q));
    incident_edges = std::map<VH, EHs, Less_vertex>(Less_vertex(geom_traits));
    edx  = std::map<EH, int, Less_edge>(Less_edge(geom_traits));

    EHs relevant_edges; //edges that can affect the visibility of query point.
    Arrangement_2 bbox;
    if (is_face_query)
      input_face(f);
    else
      input_face(f, relevant_edges, bbox);
    //the following code is the initiation of vision ray.
    //the direction of the initial ray is between the direction from q to last
    //vertex in vs and positive x-axis. By choosing this direction, we make sure
    //that all plane is swept and there is not needle at the beginning of
    //sweeping.
    Vector_2 dir;
    if (Direction_2(-1, 0) < Direction_2(Vector_2(q, vs.back()->point())))
      dir = Vector_2(1, 0) + Vector_2(q, vs.back()->point());
    else
      dir = Vector_2(0, -1);
    Point_2 dp = q + dir;

    //initiation of active_edges. for face queries,
    //all edges on the boundary can affect visibility.
    //for non-face queries, only relevant_edges has to be considered.
    if (is_face_query) {
      Ccb_halfedge_const_circulator curr = f->outer_ccb();
      Ccb_halfedge_const_circulator circ = curr;
      do {
        if (do_intersect_ray(
                    q, dp, curr->target()->point(), curr->source()->point()))
          active_edges.insert(curr);

      } while (++curr != circ);

      typename Arrangement_2::Hole_const_iterator hi;
      for (hi = f->holes_begin(); hi != f->holes_end(); ++hi) {
        Ccb_halfedge_const_circulator curr = circ = *hi;
        do {
          if (do_intersect_ray(
                      q, dp, curr->target()->point(), curr->source()->point()))
            active_edges.insert(curr);
        } while (++curr != circ);
      }
    }
    else {
      for (typename EHs::size_type i=0; i!=relevant_edges.size(); i++)
        if (do_intersect_ray(q, dp, relevant_edges[i]->source()->point(),
                                    relevant_edges[i]->target()->point()))
          active_edges.insert(relevant_edges[i]);
    }

    //angular sweep begins
//    std::cout<<active_edges.size()<<std::endl;
    for (typename VHs::size_type i=0; i!=vs.size(); i++) {
      VH vh = vs[i];
      EH closest_e = *active_edges.begin();
      EHs& edges = incident_edges[vh];
      EHs insert_ehs, remove_ehs;
      for (typename EHs::size_type j=0; j!=edges.size(); j++) {
        EH& e = edges[j];
        if (active_edges.find(e) == active_edges.end())
          insert_ehs.push_back(e);
        else
          remove_ehs.push_back(e);
      }
      typename EHs::size_type insert_cnt = insert_ehs.size();
      typename EHs::size_type remove_cnt = remove_ehs.size();
      if (insert_cnt == 1 && remove_cnt == 1) {
        const EH& ctemp_eh = *active_edges.find(remove_ehs.front());
        EH& temp_eh = const_cast<EH&>(ctemp_eh);
        temp_eh = insert_ehs.front();
      }
      else {
        for (typename EHs::size_type j=0; j!=remove_cnt; j++)
          active_edges.erase(remove_ehs[j]);
        for (typename EHs::size_type j=0; j!=insert_cnt; j++)
          active_edges.insert(insert_ehs[j]);
      }

      if (closest_e != *active_edges.begin()) {
        //when the closest edge changed
        if (is_face_query) {
          if (remove_cnt > 0 && insert_cnt > 0) {
            //some edges are added and some are deleted,
            //which means the vertex swept is part of visibility polygon.
            update_visibility(vh->point());
          }
          if (remove_cnt == 0 && insert_cnt > 0) {
            //only add some edges, means the view ray is blocked by new edges.
            //therefore first add the intersection of view ray and
            //former closet edge, then add the vertice swept.
            update_visibility(ray_seg_intersection(q,
                                                   vh->point(),
                                                   closest_e->target()->point(),
                                                   closest_e->source()->point())
                              );
            update_visibility(vh->point());
          }
          if (remove_cnt > 0 && insert_cnt == 0) {
            //only delete some edges, means some block is moved and the view ray
            //can reach the segments after the block.
            update_visibility(vh->point());
            update_visibility(
                ray_seg_intersection(q,
                                     vh->point(),
                                     (*active_edges.begin())->target()->point(),
                                     (*active_edges.begin())->source()->point()
                                     )
            );
          }
        }
        else {
          //extra work here for edge/vertex query is the index of cone_end1 and
          //cone_end2 will be recorded.
          if (remove_cnt > 0 && insert_cnt > 0) {
            //some edges are added and some are deleted,
            //which means the vertice swept is part of visibility polygon.
            if (update_visibility(vh->point())) {
              if (vh == cone_end1)
                cone_end1_idx = polygon.size()-1;
              else if (vh == cone_end2)
                cone_end2_idx = polygon.size()-1;
            }
          }
          if (remove_cnt == 0 && insert_cnt > 0) {
            //only add some edges, means the view ray is blocked by new edges.
            //therefore first add the intersection of view ray and former closet
            //edge, then add the vertice swept.
            update_visibility(ray_seg_intersection(q,
                                                   vh->point(),
                                                   closest_e->target()->point(),
                                                   closest_e->source()->point())
                              );
            if (update_visibility(vh->point())) {
              if (vh == cone_end1)
                cone_end1_idx = polygon.size()-1;
              else if (vh == cone_end2)
                cone_end2_idx = polygon.size()-1;
            }
          }
          if (remove_cnt > 0 && insert_cnt == 0) {
            //only delete some edges, means some block is removed and the vision
            //ray can reach the segments after the block.
            if (update_visibility(vh->point())) {
              if (vh == cone_end1)
                cone_end1_idx = polygon.size()-1;
              else if (vh == cone_end2)
                cone_end2_idx = polygon.size()-1;
            }
            update_visibility(
                ray_seg_intersection(q,
                                     vh->point(),
                                     (*active_edges.begin())->target()->point(),
                                     (*active_edges.begin())->source()->point())
            );
          }
        }
      }
    }
  }

  void print_edge(const EH e) const {
    std::cout << e->source()->point() <<"->"<< e->target()->point() <<std::endl;
  }

  //compute the intersection of ray(q->dp) and segment(s, t)
  //if they are collinear then return the endpoint which is closer to q.

  Point_2 ray_seg_intersection(
      const Point_2& q, const Point_2& dp,  // the ray
      const Point_2& s, const Point_2& t)   // the segment
  const
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
  }

  //check if p has been discovered before, if not update the visibility polygon
  bool update_visibility(const Point_2& p) const {
    if (polygon.empty()) {
      polygon.push_back(p);
      return true;
    }
    else if (Visibility_2::compare_xy_2(geom_traits, polygon.back(), p)
             != EQUAL)
    {
      polygon.push_back(p);
      return true;
    }
    return false;
  }

  //functor to decide which vertex is swept earlier by the rotational sweeping
  //ray
  class Is_swept_earlier:public CGAL::binary_function<VH, VH, bool> {
    const Point_2& q;
    const Geometry_traits_2* geom_traits;
  public:
    Is_swept_earlier(const Point_2& q, const Geometry_traits_2* traits) :
        q(q), geom_traits(traits) {}

    bool operator() (const VH v1, const VH v2) const {
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

    //return the quadrant of p with respect to o.
    int quadrant(const Point_2& o, const Point_2& p) const {
      typename Geometry_traits_2::Compare_x_2 compare_x =
                geom_traits->compare_x_2_object();

      typename Geometry_traits_2::Compare_y_2 compare_y =
                geom_traits->compare_y_2_object();

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

  //when the query point is in face, every edge is good.
  void input_neighbor_f( const Halfedge_const_handle e) const {
    VH v = e->target();
    if (!incident_edges.count(v))
      vs.push_back(v);
    incident_edges[v].push_back(e);
    incident_edges[v].push_back(e->next());
  }

  //check if p is in the visibility cone
  bool is_in_cone(const Point_2& p) const{
    if (is_big_cone)
      return (!CGAL::right_turn(cone_end1->point(), q, p)) ||
              (!CGAL::left_turn(cone_end2->point(), q, p));
    else
      return (!CGAL::right_turn(cone_end1->point(), q, p)) &&
              (!CGAL::left_turn(cone_end2->point(), q, p));
  }

  //for vertex and edge query: the visibility is limited in a cone.
  void input_edge(const Halfedge_const_handle e,
                  EHs& good_edges) const {
    for (typename EHs::size_type i = 0; i < bad_edges.size(); i++)
      if (e == bad_edges[i])
        return;
    VH v1 = e->target();
    VH v2 = e->source();
    //an edge will affect visibility only if it has an endpoint in the
    //visibility cone or it crosses the boundary of the cone.
    if (is_in_cone(v1->point()) || is_in_cone(v2->point()) ||
            do_intersect_ray(q, cone_end1->point(), v1->point(), v2->point()))
    {
      good_edges.push_back(e);
      if (!incident_edges.count(v1))
        vs.push_back(v1);
      incident_edges[v1].push_back(e);
      if (!incident_edges.count(v2))
        vs.push_back(v2);
      incident_edges[v2].push_back(e);
    }
  }

  //for face query: traverse the face to get all edges
  //and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh) const
  {
    Ccb_halfedge_const_circulator curr = fh->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      CGAL_assertion(curr->face() == fh);
      input_neighbor_f(curr);
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator curr = *hi, circ = *hi;
      do {
        CGAL_assertion(curr->face() == fh);
        input_neighbor_f(curr);
      } while (++curr != circ);
    }

    std::sort(vs.begin(), vs.end(), Is_swept_earlier(q, geom_traits));

    for (typename VHs::size_type i=0; i!=vs.size(); i++) {
      typename VHs::size_type j = i+1;
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
  //for vertex or edge query: traverse the face to get all edges
  //and sort vertices in counter-clockwise order.
  void input_face (Face_const_handle fh,
                   EHs& good_edges,
                   Arrangement_2& bbox) const
  {
    Ccb_halfedge_const_circulator curr = fh->outer_ccb();
    Ccb_halfedge_const_circulator circ = curr;
    do {
      CGAL_assertion(curr->face() == fh);
      input_edge(curr, good_edges);
    } while (++curr != circ);

    typename Arrangement_2::Hole_const_iterator hi;
    for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
      Ccb_halfedge_const_circulator curr = circ = *hi;
      do {
        CGAL_assertion(curr->face() == fh);
        input_edge(curr, good_edges);
      } while (++curr != circ);
    }

    //create a box that cover all vertices such that during the sweeping,
    //the vision ray will always intersect at least an edge.
    //this box doesn't intersect any relevant_edge.
    Points points;
    for (typename VHs::size_type i=0; i<vs.size(); i++) {
      points.push_back(vs[i]->point());
    }
    points.push_back(q);
    //first get the bounding box of all relevant points.
    typename Geometry_traits_2::Iso_rectangle_2 bb =
                                    bounding_box(points.begin(), points.end());

    Number_type xmin, xmax, ymin, ymax;
    typename Geometry_traits_2::Compute_x_2 compute_x =
                                            geom_traits->compute_x_2_object();

    typename Geometry_traits_2::Compute_y_2 compute_y =
                                            geom_traits->compute_y_2_object();

    //make the box a little bigger than bb so that it won't intersect any
    //relevant_edge.
    xmin = compute_x((bb.min)())-1;
    ymin = compute_y((bb.min)())-1;
    xmax = compute_x((bb.max)())+1;
    ymax = compute_y((bb.max)())+1;
    Point_2 box[4] = {Point_2(xmin, ymin), Point_2(xmax, ymin),
                      Point_2(xmax, ymax), Point_2(xmin, ymax)};

    Halfedge_handle e1 = bbox.insert_in_face_interior(Segment_2(box[0], box[1]),
                                                      bbox.unbounded_face());

    Halfedge_handle e2 = bbox.insert_from_left_vertex(Segment_2(box[1], box[2]),
                                                      e1->target());

    Halfedge_handle e3 = bbox.insert_from_right_vertex(Segment_2(box[2],box[3]),
                                                       e2->target());

    bbox.insert_at_vertices(Segment_2(box[0], box[3]),
                            e1->source(), e3->target());

    circ = curr = e1->face()->outer_ccb();
    do {
      VH v = curr->target();
      vs.push_back(v);
      incident_edges[v].push_back(curr);
      incident_edges[v].push_back(curr->next());
      good_edges.push_back(curr);
    } while(++curr != circ);

    std::sort(vs.begin(), vs.end(), Is_swept_earlier(q, geom_traits));

    for (typename VHs::size_type i=0; i!=vs.size(); i++) {
      typename VHs::size_type j = i+1;
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

};
} // end namespace CGAL



#endif


