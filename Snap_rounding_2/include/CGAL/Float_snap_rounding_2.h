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
#include <CGAL/intersection_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/box_intersection_d.h>
#include <set>
#include <vector>

namespace CGAL {
namespace Box_intersection_d {

template<class NT_, int N>
class Box_with_index_d: public Box_d< NT_, N, ID_EXPLICIT> {
protected:
    size_t m_index;
public:
    typedef Box_d< NT_, N, ID_EXPLICIT> Base;
    typedef NT_                      NT;
    typedef size_t                   ID;

    Box_with_index_d() {}
    Box_with_index_d( ID i) : m_index(i) {}
    Box_with_index_d( bool complete, ID i): Base(complete), m_index(i) {}
    Box_with_index_d(NT l[N], NT h[N], ID i) : Base( l, h), m_index(i) {}
    Box_with_index_d( const Bbox_2& b, ID i) : Base( b), m_index(i) {}
    Box_with_index_d( const Bbox_3& b, ID i) : Base( b), m_index(i) {}
    ID  index() const { return m_index; }
    // ID  id() const { return m_index; }
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
  // Total order is required to using set
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

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Sort the input points" << std::endl;
#endif
  // Compute the order of the points along the 2 axis
  // Sorted the points may perform exact computations and thus refine the intervals of the coordinates values
  // This refine ensures that the order of the points will be preserved when rounded
  // However, except for this reason these sets are unused
  using Iterator_set_x = typename std::set<size_t, decltype(comp_by_x_first)>::iterator;
  using Iterator_set_y = typename std::set<size_t, decltype(comp_by_y_first)>::iterator;
  std::set<size_t, decltype(comp_by_x_first)> p_sort_by_x(comp_by_x_first);
  std::set<size_t, decltype(comp_by_y_first)> p_sort_by_y(comp_by_y_first);
  for(size_t i=0; i!=pts.size(); ++i)
  {
    p_sort_by_x.insert(i);
    p_sort_by_y.insert(i);
  }

  //Prepare boxes for box_intersection_d
  std::vector<PBox> points_boxes;
  std::vector<SBox> segs_boxes;
  for(size_t i=0; i<pts.size(); ++i)
    points_boxes.emplace_back(pts[i].bbox(),i);
  for(size_t i=0; i<polylines.size(); ++i)
    segs_boxes.emplace_back(pts[polylines[i][0]].bbox()+pts[polylines[i][1]].bbox(),i);

  // Bound the maximum squared distance between a point and its rounded value
  auto round_bound=[](Point_2& p){
    return std::pow(std::nextafter(p.x().approx().sup(),std::numeric_limits<double>::infinity())-p.x().approx().inf(),2)+
           std::pow(std::nextafter(p.y().approx().sup(),std::numeric_limits<double>::infinity())-p.y().approx().inf(),2);
  };

  // Callback used for box_intersection_d
  auto callback=[&](PBox &bp, SBox &bseg){
    size_t pi=bp.index();
    size_t si=bseg.index();
    size_t si1=polylines[bseg.index()][0];
    size_t si2=polylines[bseg.index()][1];

    if((pi==si1) || (pi==si2))
      return;

    Point_2& p= pts[bp.index()];
    Segment_2 seg(pts[polylines[bseg.index()][0]], pts[polylines[bseg.index()][1]]);

    // (A+B)^2 <= 4*max(A^2,B^2), we take some margin
    double bound=16*(std::max)(round_bound(pts[pi]), (std::max)(round_bound(pts[si1]), round_bound(pts[si2])));

    // If the segment is closed to the vertex, we subdivide it at same x coordinate that this vertex
    if(possibly(Kernel().compare_squared_distance_2_object()(p, seg, bound)!=CGAL::LARGER) &&
       compare(seg.source().x(),p.x())!=compare(seg.target().x(),p.x()))
    {
      pts.emplace_back(p.x(), seg.supporting_line().y_at_x(p.x()));
      auto pair=p_sort_by_x.insert(pts.size()-1);
      if(pair.second)
        p_sort_by_y.insert(pts.size()-1);
      else
        pts.pop_back(); // Remove the new point if it is already exist
      polylines[si].emplace_back(*pair.first); // Some duplicates maybe introduced, it will be removed later
    }
  };

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Exhibit pairs of possible intersections" << std::endl;
#endif

  do{
    size_t size_before=pts.size();
    CGAL::box_intersection_d(points_boxes.begin(), points_boxes.end(), segs_boxes.begin(), segs_boxes.end(), callback);
    points_boxes.clear();
    // The new vertices may intersect another segment when rounded, we repeat until they are not new vertices
#ifdef DOUBLE_2D_SNAP_VERBOSE
    std::cout << points_boxes.size()-size_before << " subdivisions performed" << std::endl;
#endif
    for(size_t i=size_before; i<pts.size(); ++i)
      points_boxes.emplace_back(pts[i].bbox(),i);
  } while(points_boxes.size()!=0);

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Form the polylines" << std::endl;
#endif
  for(auto &polyline: polylines)
  {
    if(polyline.size()==2)
      continue;

    // Sort the subdivision points on the polyline along the original vector
    Vector_2 ref(pts[polyline[0]], pts[polyline[1]]);
    auto sort_along_ref=[&](size_t pi, size_t qi){
      Vector_2 v(pts[pi], pts[qi]);
      if(is_zero(ref.x()))
        return is_positive(v.y()*ref.y());
      return is_positive(v.x()*ref.x());
    };
    size_t ps=polyline[0];
    size_t pt=polyline[1];
    std::sort(polyline.begin(), polyline.end(), sort_along_ref);
    CGAL_assertion((polyline[0]==ps) && (polyline[polyline.size()-1]==pt));
  }

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Round" << std::endl;
#endif
  for(auto &p: pts)
    p=Point_2(to_double(p.x()), to_double(p.y()));

  //The order of the points on an axis must be preserved for the correctness of the algorithm
  CGAL_assertion(std::is_sorted(p_sort_by_x.begin(),p_sort_by_x.end(),comp_by_x));
  CGAL_assertion(std::is_sorted(p_sort_by_y.begin(),p_sort_by_y.end(),comp_by_y));

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Remove duplicate points" << std::endl;
#endif
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

#ifdef DOUBLE_2D_SNAP_VERBOSE
  // The algorithm prevents the a vertex that goes through a segment but a vertex may lie on an horizontal/vertical segments after rounding
  std::cout << "Subdivide horizontal and vertical segments with vertices on them" << std::endl;
#endif
  //The order may have changed (Example: (1,1)<(1,2)<(1+e,1)), we recompute it
  p_sort_by_x.clear();
  p_sort_by_y.clear();
  for(size_t i=0; i!=pts.size(); ++i)
  {
    p_sort_by_x.insert(i);
    p_sort_by_y.insert(i);
  }

  for(auto &poly: polylines){
    std::vector<size_t> updated_polyline;
    updated_polyline.push_back(poly[0]);
    for(size_t i=1; i!=poly.size(); ++i){
      if(pts[poly[i-1]].x()==pts[poly[i]].x()){
        Iterator_set_x start, end;
        // Get all vertices between the two endpoints along x order
        if(comp_by_x_first(poly[i-1],poly[i])){
          start=p_sort_by_x.upper_bound(poly[i-1]);
          end=p_sort_by_x.lower_bound(poly[i]);
        } else {
          start=p_sort_by_x.upper_bound(poly[i]);
          end=p_sort_by_x.lower_bound(poly[i-1]);
        }
        // Add all endpoints between them to the polyline
        for(auto it=start; it!=end; ++it){
          updated_polyline.push_back(*it);
        }
      }
      if(pts[poly[i-1]].y()==pts[poly[i]].y()){
        Iterator_set_y start, end;
        // Get all vertices between the two endpoints along y order
        if(comp_by_y_first(poly[i-1],poly[i])){
          start=p_sort_by_y.upper_bound(poly[i-1]);
          end=p_sort_by_y.lower_bound(poly[i]);
        } else {
          start=p_sort_by_y.upper_bound(poly[i]);
          end=p_sort_by_y.lower_bound(poly[i-1]);
        }
        // Add all endpoints between them to the polyline
        for(auto it=start; it!=end; ++it){
          updated_polyline.push_back(*it);
        }
      }
      updated_polyline.push_back(poly[i]);
    }
    std::swap(poly, updated_polyline);
  }
}


/*
TODO doc
*/
template <class InputIterator , class OutputContainer>
typename OutputContainer::iterator double_snap_rounding_2(InputIterator  	input_begin,
		                                                      InputIterator  	input_end,
		                                                      OutputContainer&  output)
{
  using Segment_2 = std::remove_cv_t<typename std::iterator_traits<InputIterator>::value_type>;
  using Polyline = std::remove_cv_t<typename std::iterator_traits<typename OutputContainer::iterator>::value_type>;
  using Point_2 = typename Default_arr_traits<Segment_2>::Traits::Point_2;

  std::vector<Segment_2> segs(input_begin, input_end);
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Solved intersections" << std::endl;
#endif
  compute_subcurves(input_begin, input_end, std::back_inserter(segs));

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indexes" << std::endl;
#endif
  std::set<Point_2> unique_point_set;
  std::map<Point_2, int> point_to_index;
  std::vector<Point_2> pts;
  std::vector< std::vector< size_t> > polylines;

  // Transform range of the segments in the range of points and polyline of indexes
  for(typename std::vector<Segment_2>::iterator it=segs.begin(); it!=segs.end(); ++it)
  {
    const Point_2& p1 = it->source();
    const Point_2& p2 = it->target();

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

  for(InputIterator it=segs.begin(); it!=segs.end(); ++it)
  {
    size_t index1 = point_to_index[it->source()];
    size_t index2 = point_to_index[it->target()];
    polylines.push_back({index1, index2});
  }

  // Main algorithm
  double_snap_rounding_2_disjoint(pts, polylines);

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Output polylines
  output.clear();
  for(auto &poly: polylines){
    Polyline new_line;
    for(size_t pi: poly)
      new_line.push_back(pts[pi]);
    output.push_back(new_line);
  }

  return output.begin();
}

/*
TODO doc
*/
template <class InputIterator , class OutputContainer>
typename OutputContainer::iterator compute_snap_subcurves_2(InputIterator  	input_begin,
		                                                  InputIterator  	input_end,
		                                                  OutputContainer&  output)
{
  using Segment_2 = std::remove_cv_t<typename std::iterator_traits<InputIterator>::value_type>;
  using Point_2 = typename Default_arr_traits<Segment_2>::Traits::Point_2;

  std::vector<Segment_2> segs;
#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Solved intersections" << std::endl;
#endif

  compute_subcurves(input_begin, input_end, std::back_inserter(segs));

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Change format to range of points and indexes" << std::endl;
#endif
  std::set<Point_2> unique_point_set;
  std::map<Point_2, int> point_to_index;
  std::vector<Point_2> pts;
  std::vector< std::vector< size_t> > polylines;

  for(typename std::vector<Segment_2>::iterator it=segs.begin(); it!=segs.end(); ++it)
  {
    const Point_2& p1 = it->source();
    const Point_2& p2 = it->target();

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

  for(InputIterator it=segs.begin(); it!=segs.end(); ++it)
  {
    size_t index1 = point_to_index[it->source()];
    size_t index2 = point_to_index[it->target()];
    polylines.push_back({index1, index2});
  }

  // Main algorithm
  double_snap_rounding_2_disjoint(pts, polylines);

#ifdef DOUBLE_2D_SNAP_VERBOSE
  std::cout << "Build output" << std::endl;
#endif

  // Output a range of segments while removing duplicate ones
  std::set< std::pair<size_t,size_t> > set_out_segs;
  output.clear();
  for(auto &poly: polylines){
    for(size_t i=1; i<poly.size(); ++i)
      set_out_segs.emplace((std::min)(poly[i-1],poly[i]),(std::max)(poly[i-1],poly[i]));
  }
  for(auto &pair: set_out_segs)
    output.emplace_back(pts[pair.first], pts[pair.second]);

  return output.begin();
}

} //namespace CGAL

#endif
