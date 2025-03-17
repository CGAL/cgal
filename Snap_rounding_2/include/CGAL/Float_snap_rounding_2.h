// Copyright (c) 2025 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Léo Valque

#ifndef CGAL_FLOAT_SNAP_ROUNDING_2_H
#define CGAL_FLOAT_SNAP_ROUNDING_2_H


#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/enum.h>
// #include <CGAL/predicates_on_points_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/box_intersection_d.h>
#include <list>
#include <set>
#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/tss.h>

namespace CGAL {
namespace Box_intersection_d {

template<class NT_, int N>
class Box_with_index_d: public Box_d< NT_, N, ID_NONE> {
protected:
    size_t m_index;
public:
    typedef Box_d< NT_, N, ID_NONE> Base;
    typedef NT_                      NT;
    typedef size_t                   ID;

    Box_with_index_d() {}
    Box_with_index_d( ID i) : m_index(i) {}
    Box_with_index_d( bool complete, ID i): Base(complete), m_index(i) {}
    Box_with_index_d(NT l[N], NT h[N], ID i) : Base( l, h), m_index(i) {}
    Box_with_index_d( const Bbox_2& b, ID i) : Base( b), m_index(i) {}
    Box_with_index_d( const Bbox_3& b, ID i) : Base( b), m_index(i) {}
    ID  index() const { return m_index; }
    ID  id() const { return m_index; }
};

// Generic template signature of boxes, specialized for ID_FROM_HANDLE policy
#if 0
template<class NT_, int N, class Handle_, class IdPolicy = ID_FROM_INDEX>
class Box_with_index_d : public Box_d< NT_, N, IdPolicy> {
protected:
    size_t m_index;
public:
    typedef Box_d< NT_, N, IdPolicy> Base;
    typedef NT_                      NT;
    typedef size_t                   ID;

    Box_with_index_d() {}
    Box_with_index_d( ID i) : m_index(i) {}
    Box_with_index_d( bool complete, ID i): Base(complete), m_index(i) {}
    Box_with_index_d(NT l[N], NT h[N], ID i) : Base( l, h), m_index(i) {}
    Box_with_index_d( const Bbox_2& b, ID i) : Base( b), m_index(i) {}
    Box_with_index_d( const Bbox_3& b, ID i) : Base( b), m_index(i) {}
    ID  index() const { return m_index; }
};

// Specialization for ID_FROM_INDEX policy
template<class NT_, int N>
class Box_with_index_d<NT_, N, Handle_, ID_FROM_INDEX>
    : public Box_d< NT_, N, ID_NONE> {
protected:
    size_t m_index;
public:
    typedef Box_d< NT_, N, ID_NONE> Base;
    typedef NT_                      NT;
    typedef size_t                   ID;

    Box_with_index_d() {}
    Box_with_index_d( ID i) : m_index(i) {}
    Box_with_index_d( bool complete, ID i): Base(complete), m_index(i) {}
    Box_with_index_d(NT l[N], NT h[N], ID i) : Base( l, h), m_index(i) {}
    Box_with_index_d( const Bbox_2& b, ID i) : Base( b), m_index(i) {}
    Box_with_index_d( const Bbox_3& b, ID i) : Base( b), m_index(i) {}
    ID  index() const { return m_index; }
    ID  id() const { return m_index; }
};
#endif
}

template <class PointsRange , class PolylineRange>
void double_snap_rounding_2_disjoint(PointsRange &pts,
		                                 PolylineRange &polylines)
{
  using Point_2 = std::remove_cv_t<typename std::iterator_traits<typename PointsRange::iterator>::value_type>;
  using Kernel = typename Kernel_traits<Point_2>::Kernel;
  using Comparison_result = typename Kernel::Comparison_result;
  using Segment_2 =  typename Kernel::Segment_2;
  using Vector_2 =  typename Kernel::Vector_2;
  using PBox = CGAL::Box_intersection_d::Box_with_index_d<double,2>;
  using SBox = CGAL::Box_intersection_d::Box_with_index_d<double,2>;

  auto comp_by_x=[&](size_t ai, size_t bi){
    return compare(pts[ai].x(),pts[bi].x())==SMALLER;
  };
  auto comp_by_y=[&](size_t ai, size_t bi){
    return compare(pts[ai].y(),pts[bi].y())==SMALLER;
  };
  auto comp_by_x_first=[&](size_t ai, size_t bi){
    Comparison_result res=compare(pts[ai].x(),pts[bi].x());
    if(res==EQUAL)
      return compare(pts[ai].y(),pts[bi].y())==SMALLER;
    return res==SMALLER;
  };
  auto comp_by_y_first=[&](size_t ai, size_t bi){
    Comparison_result res=compare(pts[ai].y(),pts[bi].y());
    if(res==EQUAL)
      return compare(pts[ai].x(),pts[bi].x())==SMALLER;
    return res==SMALLER;
  };


  // Compute the order of the points along the 2 axis
  // Sorted the points may perform exact computations and thus refine the intervals of the coordinates values
  // This refine ensures that the order of the points will be preserved by the rounding
  using Iterator_set_x = typename std::set<size_t, decltype(comp_by_x_first)>::iterator;
  using Iterator_set_y = typename std::set<size_t, decltype(comp_by_y_first)>::iterator;
  std::set<size_t, decltype(comp_by_x_first)> p_sort_by_x(comp_by_x_first);
  std::set<size_t, decltype(comp_by_y_first)> p_sort_by_y(comp_by_y_first);
  for(size_t i=0; i!=pts.size(); ++i)
  {
    p_sort_by_x.insert(i);
    p_sort_by_y.insert(i);
  }
  // std::vector<size_t> p_sort_by_y(pts.size(),0);
  // std::iota(p_sort_by_x.begin(),p_sort_by_x.end(),0);
  // std::iota(p_sort_by_y.begin(),p_sort_by_y.end(),0);
  // std::sort(p_sort_by_x.begin(),p_sort_by_x.end(),comp_by_x);
  // std::sort(p_sort_by_y.begin(),p_sort_by_y.end(),comp_by_y);

  //Kd-Tree to exhibits pairs
  std::vector<PBox> points_boxes;
  std::vector<SBox> segs_boxes;
  for(size_t i=0; i<pts.size(); ++i)
    points_boxes.emplace_back(pts[i].bbox(),i);
  for(size_t i=0; i<polylines.size(); ++i)
    segs_boxes.emplace_back(pts[polylines[i][0]].bbox()+pts[polylines[i][1]].bbox(),i);

  auto round_bound=[](Point_2& p){
    return std::pow(p.x().approx().sup()-p.x().approx().inf(),2)+std::pow(p.y().approx().sup()-p.y().approx().inf(), 2);
  };
  auto callback=[&](PBox &bp, SBox &bseg){
    size_t pi=bp.index();
    size_t si=bseg.index();
    size_t si1=polylines[bseg.index()][0];
    size_t si2=polylines[bseg.index()][1];

    // the point is a vertex of the segment
    if((pi==si1) || (pi==si2))
      return;

    Point_2& p= pts[bp.index()];
    Segment_2 seg(pts[polylines[bseg.index()][0]], pts[polylines[bseg.index()][1]]);

    double round_bound_s=(std::max)(round_bound(pts[si1]), round_bound(pts[si2]));

    if(possibly(Kernel().compare_squared_distance_2_object()(p, seg, round_bound_s)!=CGAL::LARGER)){
      pts.emplace_back(p.x(), seg.supporting_line().y_at_x(p.x()));
      auto pair=p_sort_by_x.insert(pts.size()-1);
      if(pair.second)
        p_sort_by_y.insert(pts.size()-1);
      else
        pts.pop_back();
      polylines[si].emplace_back(*pair.first);
    }
  };

  CGAL::box_intersection_d(points_boxes.begin(), points_boxes.end(), segs_boxes.begin(), segs_boxes.end(), callback);

  //sort new vertices
  for(auto &polyline: polylines)
  {
    if(polyline.size()==2)
      continue;

    //Sort the points on the polyline along the original vector
    Vector_2 ref(pts[polyline[0]], pts[polyline[1]]);
    auto sort_along_ref=[&](size_t pi, size_t qi){
      Vector_2 v(pts[pi], pts[qi]);
      if(is_zero(ref.x()))
        return is_positive(v.y()*ref.y());
      return is_positive(v.x()*ref.x());
    };
    std::sort(polyline.begin(), polyline.end(), sort_along_ref);
  }

  //round
  for(auto &p: pts)
    p=Point_2(to_double(p.x()), to_double(p.y()));

  //The order of the points on an axis must be preserved for the correctness of the algorithm
  CGAL_assertion(std::is_sorted(p_sort_by_x.begin(),p_sort_by_x.end(),comp_by_x));
  CGAL_assertion(std::is_sorted(p_sort_by_y.begin(),p_sort_by_y.end(),comp_by_y));

  //remove duplicate_points
  std::vector< size_t > unique_points(p_sort_by_x.begin(),p_sort_by_x.end());
  std::sort(unique_points.begin(),unique_points.end(),comp_by_x_first);
  std::vector<Point_2> new_pts;
  std::vector<size_t> old_to_new_index(pts.size());
  for(size_t i=0; i!=pts.size(); ++i){
    if(i==0 || (pts[unique_points[i]]!=pts[unique_points[i-1]]))
      new_pts.push_back(pts[unique_points[i]]);
    old_to_new_index[unique_points[i]]=new_pts.size()-1;
  }

  std::swap(pts, new_pts);
  // Update the polylines by remapping the old indices to new indices
  for (auto& polyline : polylines) {
    std::vector<size_t> updated_polyline;
    for (size_t i=0; i<polyline.size(); ++i) {
        size_t new_pi=old_to_new_index[polyline[i]];
        if(i==0 || (new_pi!=updated_polyline[updated_polyline.size()-1]))
          updated_polyline.push_back(new_pi);
        assert(new_pi<pts.size());
    }
    std::swap(polyline, updated_polyline);
  }

  //Vertices can be on vertical or horizontal segment, we repair this
  p_sort_by_x.clear();
  p_sort_by_y.clear();
  for(size_t i=0; i!=pts.size(); ++i)
  {
    p_sort_by_x.insert(i);
    p_sort_by_y.insert(i);
  }
  // std::vector< size_t > vec_sort_by_x_first(p_sort_by_x.begin(),p_sort_by_x.end());
  // std::vector< size_t > vec_sort_by_y_first(p_sort_by_y.begin(),p_sort_by_y.end());
  // std::sort(vec_sort_by_x_first.begin(),vec_sort_by_x_first.end(),comp_by_x_first);
  // std::sort(vec_sort_by_y_first.begin(),vec_sort_by_y_first.end(),comp_by_y_first);
  for(auto &poly: polylines){
    std::vector<size_t> updated_polyline;
    updated_polyline.push_back(poly[0]);
    for(size_t i=1; i!=poly.size(); ++i){
      if(pts[poly[i-1]].x()==pts[poly[i]].x()){
        Iterator_set_x start, end;
        if(comp_by_x_first(poly[i-1],poly[i])){
          start=p_sort_by_x.upper_bound(poly[i-1]);
          end=p_sort_by_x.lower_bound(poly[i]);
        } else {
          start=p_sort_by_x.upper_bound(poly[i]);
          end=p_sort_by_x.lower_bound(poly[i-1]);
        }
        for(auto it=start; it!=end; ++it){
          updated_polyline.push_back(*it);
        }
      }
      if(pts[poly[i-1]].y()==pts[poly[i]].y()){
        Iterator_set_y start, end;
        if(comp_by_y_first(poly[i-1],poly[i])){
          start=p_sort_by_y.upper_bound(poly[i-1]);
          end=p_sort_by_y.lower_bound(poly[i]);
        } else {
          start=p_sort_by_y.upper_bound(poly[i]);
          end=p_sort_by_y.lower_bound(poly[i-1]);
        }
        for(auto it=start; it!=end; ++it){
          updated_polyline.push_back(*it);
        }
      }
      updated_polyline.push_back(poly[i]);
    }
    std::swap(poly, updated_polyline);
  }
  // compute_subcurves(input_begin, input_end, polylines);
}

template <class InputIterator , class OutputContainer>
typename OutputContainer::iterator double_snap_rounding_2(InputIterator  	input_begin,
		                                                      InputIterator  	input_end,
		                                                      OutputContainer&  output)
{
  using Segment_2 = std::remove_cv_t<typename std::iterator_traits<InputIterator>::value_type>;
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename OutputContainer::iterator>::value_type>;
  using Point_2 = typename Default_arr_traits<Segment_2>::Traits::Point_2;

  std::vector<Segment_2> segs;
  compute_subcurves(input_begin, input_end, std::back_inserter(segs));

  std::set<Point_2> unique_point_set;
  std::map<Point_2, int> point_to_index;
  std::vector<Point_2> pts;
  std::vector< std::vector< size_t> > polylines;

  // Transform range of the segments in the range of points and polyline of indexes
  for(typename std::vector<Segment_2>::iterator it=segs.begin(); it!=segs.end(); ++it)
  {
    const Point_2& p1 = it->source();
    const Point_2& p2 = it->target();

    // Check and insert the first endpoint if it's not already added
    if (unique_point_set.find(p1) == unique_point_set.end()) {
      unique_point_set.insert(p1);
      pts.push_back(p1);
      point_to_index[p1] = pts.size() - 1;
    }
    if (unique_point_set.find(p2) == unique_point_set.end()) {
      unique_point_set.insert(p2);
      pts.push_back(p2);
      point_to_index[p2] = pts.size() - 1;
    }
  }

  for(InputIterator it=input_begin; it!=input_end; ++it)
  {
    size_t index1 = point_to_index[it->source()];
    size_t index2 = point_to_index[it->target()];
    polylines.push_back({index1, index2});
  }

  double_snap_rounding_2_disjoint(pts, polylines);
  for(auto &poly: polylines){
    Polyline new_line;
    for(size_t pi: poly){
        new_line.push_back(pts[pi]);
    }
    output.push_back(new_line);
  }

  return output.begin();
}

} //namespace CGAL

#endif
