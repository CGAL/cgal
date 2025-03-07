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

#include <CGAL/license/Snap_rounding_2.h>


#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
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

template <class PointRange, class PolylineRange>
void float_snap_rounding_2(PointRange& points,
                           PolylineRange& segments)
{
  using Point_2 = std::remove_cv_t<typename std::iterator_traits<typename PointRange::const_iterator>::value_type>;
  using Kernel = typename Kernel_traits<Point_2>::Kernel;
  using Segment_2 = typename Kernel::Segment_2;
  using Vector_2 = typename Kernel::Vector_2;
  using PBox = CGAL::Box_intersection_d::Box_with_index_d<double,2>;
  using SBox = CGAL::Box_intersection_d::Box_with_index_d<double,2>;

  auto round_bound=[](Point_2& p){
    return std::pow(p.x().approx().sup()-p.x().approx().inf(),2)+std::pow(p.y().approx().sup()-p.y().approx().inf(), 2);
  };


  //Kd-Tree to exhibits pairs
  std::vector<PBox> points_boxes;
  std::vector<SBox> segs_boxes;
  for(size_t i=0; i<points.size(); ++i)
    points_boxes.emplace_back(points[i].bbox(),i);
  for(size_t i=0; i<segments.size(); ++i)
    segs_boxes.emplace_back(points[segments[i][0]].bbox()+points[segments[i][1]].bbox(),i);

  auto callback=[&](PBox &bp, SBox &bseg){
    size_t pi=bp.index();
    size_t si=bseg.index();
    size_t si1=segments[bseg.index()][0];
    size_t si2=segments[bseg.index()][1];

    // the point is a vertex of the segment
    if((pi==si1) || (pi==si2))
      return;

    Point_2& p= points[bp.index()];
    Segment_2 seg(points[segments[bseg.index()][0]], points[segments[bseg.index()][1]]);

    // Compute the maximum distance between a point of s and its rounded version
    // TODO compute one a the beginning for performance
    // round_bound_si1= (std::max)(points_boxes[si1].bbox().xmax()-points_boxes[si1].bbox().xmin(),
    //                             points_boxes[si1].bbox().ymax()-points_boxes[si1].bbox().ymin());
    // round_bound_si2= (std::max)(points_boxes[si2].bbox().xmax()-points_boxes[si2].bbox().xmin(),
    //                             points_boxes[si2].bbox().ymax()-points_boxes[si2].bbox().ymin());
    // round_bound_s=(std::max)(round_bound_si1,round_bound_si2);
    double round_bound_s=(std::max)(round_bound(points[si1]), round_bound(points[si2]));

    if(possibly(Kernel().compare_squared_distance_2_object()(p, seg, round_bound_s)!=CGAL::LARGER)){
      points.emplace_back(p.x(), seg.supporting_line().y_at_x(p.x()));
      segments[si].emplace_back(points.size()-1);
    }
  };

  CGAL::box_intersection_d(points_boxes.begin(), points_boxes.end(), segs_boxes.begin(), segs_boxes.end(), callback);

  //sort new vertices
  for(auto &seg: segments)
  {
    if(seg.size()==2)
      continue;
    Vector_2 ref(points[seg[0]], points[seg[1]]);
    auto sort_along_ref=[&](size_t pi, size_t qi){
      Vector_2 v(points[pi], points[qi]);
      if(is_zero(ref.x()))
        return is_positive(v.y()*ref.y());
      return is_positive(v.x()*ref.x());
    };
    //Il faut conserver l'ordre des deux sommets d'entrées
    std::sort(seg.begin(), seg.end(), sort_along_ref);
  }

  //round
  // TODO bug potentiel, l'algo suppose que deux coordonnées en x ne peuvent pas s'inverser mais
  // l'intervalle d'un nombre pourrait être très petit et l'intervalle d'un autre très grand et provoquer ce genre
  // d'inversion
  for(auto &p: points)
    p=Point_2(to_double(p.x()), to_double(p.y()));

  //repair some vertices can be on segment

}

} //namespace CGAL

#endif
