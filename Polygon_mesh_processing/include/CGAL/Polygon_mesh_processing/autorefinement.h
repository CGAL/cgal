// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
#define CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H

#include <CGAL/license/Polygon_mesh_processing/autorefinement.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#ifdef CGAL_AUTOREF_USE_PROGRESS_DISPLAY
#include <boost/timer/progress_display.hpp>
#endif

// output
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>


#ifdef CGAL_PMP_AUTOREFINE_USE_DEFAULT_VERBOSE
#define CGAL_PMP_AUTOREFINE_VERBOSE(X) std::cout << X << "\n";
#endif

#ifndef CGAL_PMP_AUTOREFINE_VERBOSE
#define CGAL_PMP_AUTOREFINE_VERBOSE(MSG)
#endif

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/concurrent_vector.h>
#if TBB_INTERFACE_VERSION < 12010 && !defined(TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS)
#define CGAL_HAS_DEFINED_TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS
#define TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS 1
#endif
#include <tbb/concurrent_map.h>
#include <tbb/parallel_for.h>
#ifdef CGAL_AUTOREF_SET_POINT_IDS_USING_MUTEX
#include <CGAL/mutex.h>
#endif
#endif

#ifdef CGAL_AUTOREF_USE_FIXED_PROJECTION_TRAITS
#include <CGAL/Kernel_23/internal/Projection_traits_3.h>
#endif

#include <vector>

//#define CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
//#define CGAL_AUTOREFINE_DEBUG_COUNTERS
//#define CGAL_AUTOREF_DEBUG_DEPTH

#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
#include <CGAL/Real_timer.h>
#endif

#if defined(CGAL_AUTOREFINE_DEBUG_COUNTERS) || defined(CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS)
#include <CGAL/Real_timer.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

namespace Autorefinement {

/** \ingroup PMP_corefinement_grp
 *  %Default visitor model of `PMPAutorefinementVisitor`.
 *  All of its functions have an empty body. This class can be used as a
 *  base class if only some of the functions of the concept require to be
 *  overridden.
 */
struct Default_visitor
{
  inline void number_of_output_triangles(std::size_t /*nbt*/) {}
  inline void verbatim_triangle_copy(std::size_t /*tgt_id*/, std::size_t /*src_id*/) {}
  inline void new_subtriangle(std::size_t /*tgt_id*/, std::size_t /*src_id*/) {}
};

} // end of Autorefinement visitor


#ifndef DOXYGEN_RUNNING
namespace autorefine_impl {

enum Segment_inter_type { NO_INTERSECTION=0,
                          POINT_INTERSECTION,
                          POINT_P,
                          POINT_Q,
                          POINT_R,
                          POINT_S,
                          COPLANAR_SEGMENT_PQ,
                          COPLANAR_SEGMENT_RS,
                          COPLANAR_SEGMENT_PS,
                          COPLANAR_SEGMENT_QS,
                          COPLANAR_SEGMENT_PR,
                          COPLANAR_SEGMENT_QR,
                        };

// test intersection in the interior of segment pq and rs with pq and rs being coplanar segments
// note that for coplanar cases, we might report identical endpoints
template <class K>
Segment_inter_type
do_coplanar_segments_intersect(std::size_t pi, std::size_t qi,
                               std::size_t ri, std::size_t si,
                               const std::vector<typename K::Point_3>& points,
                               const K& k = K())
{
  typename K::Collinear_are_ordered_along_line_3 cln_order = k.collinear_are_ordered_along_line_3_object();
  typename K::Coplanar_orientation_3 cpl_orient=k.coplanar_orientation_3_object();

  const typename K::Point_3& p=points[pi];
  const typename K::Point_3& q=points[qi];
  const typename K::Point_3& r=points[ri];
  const typename K::Point_3& s=points[si];

  // first handle case of shared endpoints
  if (pi==ri)
  {
    if (si==qi || cpl_orient(p, q, s)!=COLLINEAR) return NO_INTERSECTION;
    // can be s, q or nothing
    if (cln_order(p,s,q))
      return POINT_S;
    if (cln_order(p,q,s))
      return POINT_Q;
    return NO_INTERSECTION;
  }
  else
  {
    if(pi==si)
    {
      if (qi==ri || cpl_orient(p, q, r)!=COLLINEAR) return NO_INTERSECTION;
      // can be r, q or nothing
      if (cln_order(p,r,q))
        return POINT_R;
      if (cln_order(p,q,r))
        return POINT_Q;
      return NO_INTERSECTION;
    }
    else
    {
      if (qi==ri)
      {
        if (pi==si || cpl_orient(p, q, s)!=COLLINEAR) return NO_INTERSECTION;
        // can be p, s or nothing
        if (cln_order(p,s,q))
          return POINT_S;
        if (cln_order(q,p,s))
          return POINT_P;
        return NO_INTERSECTION;
      }
      else
      {
        if (qi==si)
        {
          if (pi==ri || cpl_orient(p, q, r)!=COLLINEAR) return NO_INTERSECTION;
          // can be p, r or nothing
          if (cln_order(p,r,q))
            return POINT_R;
          if (cln_order(q,p,r))
            return POINT_P;
          return NO_INTERSECTION;
        }
      }
    }
  }

  // supporting_line intersects: points are coplanar
  ::CGAL::Orientation pqr = cpl_orient(p, q, r);
  ::CGAL::Orientation pqs = cpl_orient(p, q, s);

  if(pqr == COLLINEAR && pqs == COLLINEAR)
  {
    // segments are collinear
    bool r_in_pq = cln_order(p, r, q),
         s_in_pq = cln_order(p, s, q),
         p_in_rs = cln_order(r, p, s);

    if (r_in_pq)
    {
      // intersection could be rs, pr or qr
      if (s_in_pq)
        return COPLANAR_SEGMENT_RS;
      if (p_in_rs)
        return COPLANAR_SEGMENT_PR;
      CGAL_assertion(cln_order(r, q, s));
      return COPLANAR_SEGMENT_QR;
    }
    else
    {
      if (s_in_pq)
      {
        // intersection could be ps or qs
        if (p_in_rs)
          return COPLANAR_SEGMENT_PS;
        CGAL_assertion(cln_order(r, q, s));
        return COPLANAR_SEGMENT_QS;
      }
      else
        if (p_in_rs)
        {
          CGAL_assertion(cln_order(r, q, s));
          return COPLANAR_SEGMENT_PQ;
        }
    }
    return NO_INTERSECTION;
  }

  if(pqr != pqs)
  {
    ::CGAL::Orientation rsp = cpl_orient(r, s, p);

    if (rsp==COLLINEAR)
    {
      if (pqr==COLLINEAR || pqs==COLLINEAR)
      {
        throw std::runtime_error("no expected #1");
      }
      return POINT_P;
    }
    ::CGAL::Orientation rsq = cpl_orient(r, s, q);
    if (rsq==COLLINEAR)
    {
      if (pqr==COLLINEAR || pqs==COLLINEAR)
      {
        throw std::runtime_error("no expected #2");
      }
      return POINT_Q;
    }
    if (rsp!=rsq)
    {
      if (pqr==COLLINEAR) return POINT_R;
      if (pqs==COLLINEAR) return POINT_S;
      return POINT_INTERSECTION;
    }
  }

  return NO_INTERSECTION;
}


// imported from Intersections_3/include/CGAL/Intersections_3/internal/Triangle_3_Triangle_3_intersection.h
template<class K>
void coplanar_intersections(const std::array<typename K::Point_3, 3>& t1,
                            const std::array<typename K::Point_3, 3>& t2,
                            std::vector<std::tuple<typename K::Point_3, int, int>>& inter_pts)
{
  const typename K::Point_3& p1 = t1[0], &q1 = t1[1], &r1 = t1[2];
  const typename K::Point_3& p2 = t2[0], &q2 = t2[1], &r2 = t2[2];

  std::list<Intersections::internal::Point_on_triangle<K>> l_inter_pts;
  l_inter_pts.push_back(Intersections::internal::Point_on_triangle<K>(0,-1));
  l_inter_pts.push_back(Intersections::internal::Point_on_triangle<K>(0,-2));
  l_inter_pts.push_back(Intersections::internal::Point_on_triangle<K>(0,-3));
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  auto print_points = [&]()
  {
    std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
    for(auto p : l_inter_pts) {std::cout << "  (" << p.id1() << "," << p.id2() << ",[" << p.alpha << "]) ";} std::cout <<"\n";
  };
#endif

  //intersect t2 with the three half planes which intersection defines t1
  K k;
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  print_points();
#endif
  Intersections::internal::intersection_coplanar_triangles_cutoff(p1,q1,r1,1,p2,q2,r2,k,l_inter_pts); //line p1q1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  print_points();
#endif
  Intersections::internal::intersection_coplanar_triangles_cutoff(q1,r1,p1,2,p2,q2,r2,k,l_inter_pts); //line q1r1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  print_points();
#endif
  Intersections::internal::intersection_coplanar_triangles_cutoff(r1,p1,q1,3,p2,q2,r2,k,l_inter_pts); //line r1p1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  print_points();
#endif

  for (const Intersections::internal::Point_on_triangle<K>& pot : l_inter_pts)
    inter_pts.emplace_back(pot.point(p1,q1,r1,p2,q2,r2,k), pot.id1(), pot.id2());
}

// imported from Polygon_mesh_processing/internal/Corefinement/intersect_triangle_and_segment_3.h
template<class K>
std::optional<std::pair<typename K::Point_3, int>>
find_intersection(const typename K::Point_3& p, const typename K::Point_3& q,  //segment
                  const typename K::Point_3& a, const typename K::Point_3& b, const typename K::Point_3& c, //triangle
                  bool is_p_coplanar=false, bool is_q_coplanar=false) // note that in coref this was wrt a halfedge not p/q
{
  Orientation ab=orientation(p,q,a,b);
  Orientation bc=orientation(p,q,b,c);
  Orientation ca=orientation(p,q,c,a);

  if ( ab==POSITIVE || bc==POSITIVE || ca==POSITIVE )
    return std::nullopt;

  int nb_coplanar=(ab==COPLANAR?1:0) + (bc==COPLANAR?1:0) + (ca==COPLANAR?1:0);

  if (nb_coplanar==2)
  {
    // even if a common point is not new it is still needed to be reported so
    // that the intersection segment is known.
    if (ab!=COPLANAR)
      // intersection is c
      return std::make_pair(c, -3);
    else if (bc!=COPLANAR)
      // intersection is a
      return std::make_pair(a, -1);
    else
    {
      CGAL_assertion(ca!=COPLANAR);
      // intersection is b
      return std::make_pair(b, -2);
    }
  }

  typename K::Point_3 ipt = is_p_coplanar ? p :
                            is_q_coplanar ? q :
                            typename K::Construct_plane_line_intersection_point_3()
                              (a, b, c, p, q);

  if (nb_coplanar == 0)
    return std::make_pair(ipt, 0);


  CGAL_assertion(nb_coplanar==1);

  if (ab==COPLANAR)
    // intersection is ab
    return std::make_pair(ipt, 1);
  if (bc==COPLANAR)
    // intersection is bc
    return std::make_pair(ipt, 2);
  CGAL_assertion(ca==COPLANAR);
  // intersection is ca
  return std::make_pair(ipt, 3);
}

template <class K>
std::optional<std::pair<typename K::Point_3, int>>
test_edge(const typename K::Point_3& p, const typename K::Point_3& q,
          const typename K::Point_3& a, const typename K::Point_3& b, const typename K::Point_3& c,
          const Orientation abcp,
          const Orientation abcq)
{
  switch ( abcp ) {
  case POSITIVE:
    switch ( abcq ) {
    case POSITIVE:
      // the segment lies in the positive open halfspaces defined by the
      // triangle's supporting plane
    return std::nullopt;
    case NEGATIVE:
      // p sees the triangle in counterclockwise order
      return find_intersection<K>(p,q,a,b,c);
    //case COPLANAR:
    default:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in counterclockwise order
      return find_intersection<K>(p,q,a,b,c,false,true);
    }
  case NEGATIVE:
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return find_intersection<K>(q,p,a,b,c);
    case NEGATIVE:
      // the segment lies in the negative open halfspaces defined by the
      // triangle's supporting plane
    return std::nullopt;
    // case COPLANAR:
    default:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in clockwise order
      return find_intersection<K>(q,p,a,b,c,true,false);
    }
  default:
  //case COPLANAR: // p belongs to the triangle's supporting plane
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return find_intersection<K>(q,p,a,b,c,false, true);
    case NEGATIVE:
      // q sees the triangle in clockwise order
      return find_intersection<K>(p,q,a,b,c,true);
    //case COPLANAR:
    default:
      // the segment is coplanar with the triangle's supporting plane
      // we test whether the segment intersects the triangle in the common
      // supporting plane
      //if ( ::CGAL::Intersections::internal::do_intersect_coplanar(a,b,c,p,q,K()) )
      //{
        //handle coplanar intersection
        // nothing done as coplanar case handle in collect_intersections
        // and other intersection points will be collected with non-coplanar edges
      //}
      return std::nullopt;
    }
  }
}

template <class K>
bool collect_intersections(const std::array<typename K::Point_3, 3>& t1,
                           const std::array<typename K::Point_3, 3>& t2,
                           std::vector<std::tuple<typename K::Point_3, int, int>>& inter_pts)
{
  // test edges of t1 vs t2
  std::array<Orientation,3> ori;
  for (int i=0; i<3; ++i)
    ori[i] = orientation(t2[0],t2[1],t2[2],t1[i]);

  if (ori[0]== COPLANAR && ori[1]==COPLANAR && ori[2]==COPLANAR)
  {

    coplanar_intersections<K>(t1, t2, inter_pts);
#ifdef CGAL_AUTOREF_DEBUG_DEPTH
    for (auto p : inter_pts)
      if (depth(p)>2) throw std::runtime_error("Depth is not 4: "+std::to_string(depth(p)));
#endif

    return true;
  }

  for (int i=0; i<3; ++i)
  {
    int j=(i+1)%3;
    std::optional<std::pair<typename K::Point_3, int> > opt =
      test_edge<K>(t1[i], t1[j], t2[0], t2[1], t2[2], ori[i], ori[j]);
    if (opt)
    {
      if (ori[i]==COPLANAR)
        inter_pts.emplace_back(opt->first, -(i+1), opt->second);
      else if (ori[j]==COPLANAR)
        inter_pts.emplace_back(opt->first, -(j+1), opt->second);
      else
        inter_pts.emplace_back(opt->first,   i+1 , opt->second);
    }
  }

  // test edges of t2 vs t1
  for (int i=0; i<3; ++i)
    ori[i] = orientation(t1[0],t1[1],t1[2],t2[i]);
  for (int i=0; i<3; ++i)
  {
    int j=(i+1)%3;
    std::optional<std::pair<typename K::Point_3, int> > opt =
      test_edge<K>(t2[i], t2[j], t1[0], t1[1], t1[2], ori[i], ori[j]);
    if (opt)
    {
      if (ori[i]==COPLANAR)
        inter_pts.emplace_back(opt->first, opt->second, -(i+1));
      else if (ori[j]==COPLANAR)
        inter_pts.emplace_back(opt->first, opt->second, -(j+1));
      else
        inter_pts.emplace_back(opt->first, opt->second,   i+1 );
    }
  }

//  #warning TODO get rid of sort and unique calls
  // because we don't handle intersection type and can have edge-edge edge-vertex duplicates
  std::sort(inter_pts.begin(), inter_pts.end(), [](auto p, auto q){return std::get<0>(p)<std::get<0>(q);});
  auto last = std::unique(inter_pts.begin(), inter_pts.end(), [](auto p, auto q){return std::get<0>(p)==std::get<0>(q);});
  inter_pts.erase(last, inter_pts.end());

#ifdef CGAL_AUTOREF_DEBUG_DEPTH
  for (auto p : inter_pts)
    if (depth(p)>2) throw std::runtime_error("Depth is not 2: "+std::to_string(depth(p)));
#endif

  return false;
}

struct Triangle_data
{
  using Point_3 = Exact_predicates_exact_constructions_kernel::Point_3;
  std::vector<Point_3> points;
  std::vector<int> point_locations;
  std::vector<std::pair<std::size_t,std::size_t>> segments;
  std::vector<std::size_t> segment_input_triangle_ids;
  template <int i>
  std::size_t add_point(const std::tuple<Point_3,int, int>& tpl)
  {
    if (std::get<i>(tpl) < 0)
      return -std::get<i>(tpl)-1;
    points.push_back(std::get<0>(tpl));
    point_locations.push_back(std::get<i>(tpl));
    return points.size()-1;
  }

#ifndef NDEBUG
  bool on_same_edge(std::size_t i, std::size_t j)
#else
  bool on_same_edge_debug(std::size_t i, std::size_t j)
#endif
  {
    if (point_locations[i]==0 || point_locations[j]==0) return false;
    if (point_locations[i]>0)
    {
      if (point_locations[j]>0)
        return point_locations[i]==point_locations[j];
      return -point_locations[j]==point_locations[i] || point_locations[i]%3+1==-point_locations[j];
    }
    if (point_locations[j]<0)
      return true;
    return -point_locations[i]==point_locations[j] || point_locations[j]%3+1==-point_locations[i];
  }
#ifdef NDEBUG
  bool on_same_edge(std::size_t i, std::size_t j)
  {
    bool on = on_same_edge_debug(i,j);
    if (!on) return false;
    int eid = point_locations[i]>0 ? point_locations[i] : point_locations[j];
    if (eid<0) return true;
    if (!(collinear(points[eid-1], points[(eid)%3], points[i])))
    {
      std::cout << point_locations[i] << " " << point_locations[j] << " " << eid << "\n";
    }
    if (!(collinear(points[eid-1], points[(eid)%3], points[j])))
    {
      std::cout << point_locations[i] << " " << point_locations[j] << " " << eid << "\n";
    }
    CGAL_assertion(collinear(points[eid-1], points[(eid)%3], points[i]));
    CGAL_assertion(collinear(points[eid-1], points[(eid)%3], points[j]));
    return true;
  }
#endif
};

template <class EK,
#ifdef CGAL_AUTOREF_USE_FIXED_PROJECTION_TRAITS
          int dim,
#endif
          class PointVector>
void generate_subtriangles(std::size_t ti,
                           Triangle_data& triangle_data,
                           const std::set<std::pair<std::size_t, std::size_t> >& intersecting_triangles,
                           const std::set<std::pair<std::size_t, std::size_t> >& coplanar_triangles,
                           const std::vector<std::array<typename EK::Point_3,3>>& triangles,
                           PointVector& new_triangles
                           )
{
#ifdef CGAL_AUTOREFINE_DEBUG_COUNTERS
  struct Counter
  {
    int c1=0;
    int c2=0;
    int c3=0;
    int c4=0;
    int c5=0;
    int total=0;
    ~Counter()
    {
      std::cout << "intersection of 3 planes: " << c1 << "\n";
      std::cout << "coplanar segment intersection : " << c2 << "\n";
      std::cout << "coplanar segment overlap: " << c3 << "\n";
      std::cout << "no intersection: " << c4 << "\n";
      std::cout << "intersection filtered with bboxes: " << c5 << "\n";
      std::cout << "# pairs of segments : " << total << "\n";
      std::cout << "time computing segment intersections: " << timer1.time() << "\n";
      std::cout << "time sorting intersection points: " << timer2.time() << "\n";
      std::cout << "time for cdt of constraints: " << timer3.time() << "\n";
      std::cout << "time coplanar segment intersections: " << timer4.time() << "\n";
      std::cout << "time of do_coplanar_segments_intersect: " << timer5.time() << "\n";
      std::cout << "time of triplane intersection: " << timer6.time() << "\n";
    }
    CGAL::Real_timer timer1, timer2, timer3, timer4, timer5, timer6;
  };

  static Counter counter;
#define CGAL_AUTOREF_COUNTER_INSTRUCTION(X) X
#else
#define CGAL_AUTOREF_COUNTER_INSTRUCTION(X)
#endif

  std::vector<Exact_predicates_exact_constructions_kernel::Point_3>& points=triangle_data.points;
  std::vector<int>& point_locations=triangle_data.point_locations;
  std::vector<std::pair<std::size_t, std::size_t>>& segments=triangle_data.segments;
  std::vector<std::size_t>& in_triangle_ids=triangle_data.segment_input_triangle_ids;

// pre-compute segment intersections
  if (!segments.empty())
  {
    std::size_t nbs = segments.size();

    std::vector< std::vector<std::size_t> > points_on_segments(nbs);

    CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer1.start();)

    std::map<typename EK::Point_3, std::size_t> point_id_map;
//TODO: we already have sorted the points while deduplicating segments!
    for (std::size_t pid=0; pid<points.size(); ++pid)
    {
      CGAL_assertion_code(bool insertion_ok =)
      point_id_map.insert(std::make_pair(points[pid], pid))
      CGAL_assertion_code(.second);
      CGAL_assertion(insertion_ok);
    }

    auto get_point_id = [&](const typename EK::Point_3& pt)
    {
      auto insert_res = point_id_map.insert(std::make_pair(pt, points.size()));
      if (insert_res.second)
      {
        points.push_back(pt);
        point_locations.push_back(0);
      }
      return insert_res.first->second;
    };

    std::vector<Bbox_3> point_boxes(points.size());
    for (std::size_t i = 0; i<points.size(); ++i)
      point_boxes[i]=points[i].bbox();
    std::vector<Bbox_3> segment_boxes(nbs);
    for (std::size_t i = 0; i<nbs; ++i)
      segment_boxes[i]=point_boxes[segments[i].first]+point_boxes[segments[i].second];

    for (std::size_t i = 0; i<nbs-1; ++i)
    {
      for (std::size_t j = i+1; j<nbs; ++j)
      {
        if (intersecting_triangles.count(CGAL::make_sorted_pair(in_triangle_ids[i], in_triangle_ids[j]))!=0)
        {
          if ( !do_overlap(segment_boxes[i], segment_boxes[j]) )
          {
            CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c5;)
            CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.total;)
            continue;
          }

          CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer5.start();)
          Segment_inter_type seg_inter_type =
            do_coplanar_segments_intersect<EK>(segments[i].first, segments[i].second,
                                               segments[j].first, segments[j].second,
                                               points);
          CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer5.stop();)

          switch(seg_inter_type)
          {
            case POINT_P:
            {
              points_on_segments[j].push_back(segments[i].first);
              break;
            }
            case POINT_Q:
            {
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case POINT_R:
            {
              points_on_segments[i].push_back(segments[j].first);
              break;
            }
            case POINT_S:
            {
              points_on_segments[i].push_back(segments[j].second);
              break;
            }
            case POINT_INTERSECTION:
            {
              if (   coplanar_triangles.count(CGAL::make_sorted_pair(in_triangle_ids[i], in_triangle_ids[j])) == 0
                  && coplanar_triangles.count(CGAL::make_sorted_pair(ti, in_triangle_ids[j])) == 0
                  && coplanar_triangles.count(CGAL::make_sorted_pair(in_triangle_ids[i], ti)) == 0)
              {
                CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer6.start();)
                typename EK::Point_3 pt = typename EK::Construct_planes_intersection_point_3()(
                  triangles[in_triangle_ids[i]][0], triangles[in_triangle_ids[i]][1],triangles[in_triangle_ids[i]][2],
                  triangles[in_triangle_ids[j]][0], triangles[in_triangle_ids[j]][1],triangles[in_triangle_ids[j]][2],
                  triangles[ti][0], triangles[ti][1],triangles[ti][2]);
                CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer6.stop();)

                CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c1;)
                std::size_t pid = get_point_id(pt);
                points_on_segments[i].push_back(pid);
                points_on_segments[j].push_back(pid);
              }
              else
              {
                CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c2;)
                CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer4.start();)
                typename EK::Point_3 pt = typename EK::Construct_coplanar_segments_intersection_point_3()(
                  points[segments[i].first], points[segments[i].second],
                  points[segments[j].first], points[segments[j].second]);

                std::size_t pid = get_point_id(pt);
                points_on_segments[i].push_back(pid);
                points_on_segments[j].push_back(pid);
                CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer4.stop();)
              }
              break;
            }
            case COPLANAR_SEGMENT_PQ:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[j].push_back(segments[i].first);
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case COPLANAR_SEGMENT_RS:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].first);
              points_on_segments[i].push_back(segments[j].second);
              break;
            }
            case COPLANAR_SEGMENT_PR:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].first);
              points_on_segments[j].push_back(segments[i].first);
              break;
            }
            case COPLANAR_SEGMENT_QS:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].second);
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case COPLANAR_SEGMENT_PS:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].second);
              points_on_segments[j].push_back(segments[i].first);
              break;
            }
            case COPLANAR_SEGMENT_QR:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].first);
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case NO_INTERSECTION:
            {
              CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.c4;)
            }
          }
        }
        CGAL_AUTOREF_COUNTER_INSTRUCTION(++counter.total;)
      }
    }
    CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer1.stop();)
    CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer2.start();)
    std::size_t nb_new_segments=0;
    for (std::size_t i = 0; i<nbs; ++i)
    {
      if(!points_on_segments[i].empty())
      {
        int coord = 0;
        std::size_t src_id = segments[i].first, tgt_id = segments[i].second;
        typename EK::Point_3 src = points[src_id], tgt=points[tgt_id];
        if (src.x()==tgt.x())
        {
          coord=1;
          if (src.y()==tgt.y())
            coord=2;
        }
        if (src[coord]>tgt[coord])
        {
          std::swap(src_id, tgt_id);
          std::swap(src, tgt);
        }

        points_on_segments[i].push_back(src_id);
        std::swap(points_on_segments[i].front(), points_on_segments[i].back());
        std::sort(std::next(points_on_segments[i].begin()), points_on_segments[i].end(),
                  [&](std::size_t id1, std::size_t id2)
                  {
                    if (id1==id2) return false;
                    return points[id1][coord]<points[id2][coord];
                  });
        points_on_segments[i].push_back(tgt_id);
        auto last = std::unique(points_on_segments[i].begin(), points_on_segments[i].end());
        points_on_segments[i].erase(last, points_on_segments[i].end());
        nb_new_segments+=points_on_segments[i].size()-2;
      }
    }
    CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer2.stop();)

    // now fill segments with new segments
    segments.reserve(segments.size()+nb_new_segments);
    for (std::size_t i = 0; i<nbs; ++i)
    {
      if(!points_on_segments[i].empty())
      {
        segments[i]=std::make_pair(points_on_segments[i][0], points_on_segments[i][1]);
        for(std::size_t pos=1; pos<points_on_segments[i].size()-1; ++pos)
          segments.emplace_back(points_on_segments[i][pos], points_on_segments[i][pos+1]);
      }
    }

    CGAL_assertion(points.size()==std::set<typename EK::Point_3>(points.begin(), points.end()).size());
    CGAL_assertion(points.size()==point_id_map.size());
  }

  for (std::pair<std::size_t, size_t>& s : segments)
    if (s.second < s.first)
      std::swap(s.first,s.second);
  std::sort(segments.begin(), segments.end());
  auto last = std::unique(segments.begin(), segments.end());
  segments.erase(last, segments.end());

// init CDT + insert points and constraints
  CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer3.start();)
#ifdef CGAL_AUTOREF_USE_FIXED_PROJECTION_TRAITS
  typedef ::CGAL::internal::Projection_traits_3<EK, dim> P_traits;
#else
  typedef CGAL::Projection_traits_3<EK> P_traits;
#endif

  typedef CGAL::No_constraint_intersection_tag Itag;

  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits ,Default, Itag> CDT_2;
  //typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT;
  typedef CDT_2 CDT;

  const std::array<typename EK::Point_3,3>& t = triangles[ti];
  std::vector<typename CDT::Vertex_handle> vhandles(triangle_data.points.size());

#ifdef CGAL_AUTOREF_USE_FIXED_PROJECTION_TRAITS
  P_traits cdt_traits;
  bool orientation_flipped = false;
  CDT cdt(cdt_traits);
  // TODO: still need to figure out why I can't make the orientation_flipped correctly
  vhandles[0]=cdt.insert(t[0]);
  vhandles[1]=cdt.insert(t[1]);
  vhandles[2]=cdt.insert(t[2]);
#else
  // positive triangle normal
  typename EK::Vector_3 n = normal(t[0], t[1], t[2]);
  typename EK::Point_3 o(CGAL::ORIGIN);

  bool orientation_flipped = false;
  if ( typename EK::Less_xyz_3()(o+n,o) )
  {
    n=-n;
    orientation_flipped = true;
  }
  P_traits cdt_traits(n);

  CDT cdt(cdt_traits);
  vhandles[0]=cdt.insert_outside_affine_hull(t[0]);
  vhandles[1]=cdt.insert_outside_affine_hull(t[1]);
  vhandles[2] = cdt.tds().insert_dim_up(cdt.infinite_vertex(), orientation_flipped);
  vhandles[2]->set_point(t[2]);
#endif

  // insert points and fill vhandles
#if 0
  std::vector<std::size_t> indices(triangle_data.points.size()-3);
  std::iota(indices.begin(), indices.end(), 3);
#else
  //start by points on edges
  std::array<std::vector<std::size_t>, 3> indices_on_edges;
  std::vector<std::size_t> indices;
  for (std::size_t i=3; i<points.size(); ++i)
  {
    if (point_locations[i]==0)
    {
      CGAL_assertion(!collinear(points[0], points[1], points[i]));
      CGAL_assertion(!collinear(points[1], points[2], points[i]));
      CGAL_assertion(!collinear(points[2], points[0], points[i]));
      indices.push_back(i);
    }
    else
    {
      CGAL_assertion(point_locations[i]>0 && point_locations[i]<4);
      indices_on_edges[point_locations[i]-1].push_back(i);
    }
  }

  //sort points on edges and insert them
  typename CDT::Vertex_handle vinf=cdt.infinite_vertex();
  for (int i=0; i<3; ++i)
  {
    if (indices_on_edges[i].empty()) continue;
    std::size_t src_id=i, tgt_id=(i+1)%3;
    //look for a sort axis
    int coord = 0;
    if (points[src_id].x()==points[tgt_id].x())
    {
      coord=1;
      if (points[src_id].y()==points[tgt_id].y())
        coord=2;
    }
    if (points[src_id][coord]>points[tgt_id][coord])
      std::swap(src_id, tgt_id);

    std::sort(indices_on_edges[i].begin(), indices_on_edges[i].end(),
              [&](std::size_t id1, std::size_t id2)
              {
                return points[id1][coord]<points[id2][coord];
              });
    if (segments.empty()) // we might have intersection points duplicated
    {
      auto it=
        std::unique(indices_on_edges[i].begin(), indices_on_edges[i].end(),
                    [&](std::size_t id1, std::size_t id2)
                    {
                      return points[id1][coord]==points[id2][coord];
                    });
      indices_on_edges[i].erase(it, indices_on_edges[i].end());
    }

    std::size_t prev_id = src_id;
    for(std::size_t id : indices_on_edges[i])
    {
      CGAL_assertion(collinear(points[src_id], points[tgt_id], points[id]));
      typename CDT::Face_handle fh;
      CGAL_assertion_code(bool ok =)
      cdt.is_face(vhandles[prev_id], vhandles[tgt_id],vinf, fh);
      CGAL_assertion(ok);
      vhandles[id]=cdt.insert_in_edge(points[id], fh, fh->index(vinf));
      cdt.restore_Delaunay(vhandles[id]); // TODO maybe not each time but one global?
      CGAL_assertion(cdt.is_valid());
      prev_id=id;
    }
  }
#endif
  // then points in the interior
  typedef typename Pointer_property_map<typename EK::Point_3>::type Pmap;
  typedef Spatial_sort_traits_adapter_2<P_traits,Pmap> Search_traits;
  spatial_sort(indices.begin(), indices.end(),
               Search_traits(make_property_map(points), cdt_traits));

  typename CDT::Face_handle hint;
  for (std::size_t i : indices)
  {
    vhandles[i] = cdt.insert(points[i], hint);
    hint=vhandles[i]->face();
  }

  for (const std::pair<std::size_t, std::size_t>& ids : triangle_data.segments)
  {
    CGAL_assertion(ids.first < vhandles.size());
    CGAL_assertion(ids.second < vhandles.size());
    CGAL_assertion( vhandles[ids.first]!= typename CDT::Vertex_handle() );
    CGAL_assertion( vhandles[ids.second]!= typename CDT::Vertex_handle() );
    cdt.insert_constraint(vhandles[ids.first], vhandles[ids.second]);
  }
  CGAL_AUTOREF_COUNTER_INSTRUCTION(counter.timer3.stop();)

// extract new triangles
  for (typename CDT::Face_handle fh : cdt.finite_face_handles())
  {
    if (orientation_flipped)
      new_triangles.push_back( { CGAL::make_array(fh->vertex(0)->point(),
                                                  fh->vertex(cdt.cw(0))->point(),
                                                  fh->vertex(cdt.ccw(0))->point()), ti } );
    else
      new_triangles.push_back( { CGAL::make_array(fh->vertex(0)->point(),
                                                  fh->vertex(cdt.ccw(0))->point(),
                                                  fh->vertex(cdt.cw(0))->point()), ti } );
  }
}

} // end of autorefine_impl
#endif

/**
* \ingroup PMP_corefinement_grp
*
* refines a soup of triangles so that no pair of triangles intersects.
* Output triangles may share a common edge or a common vertex (but with the same indexed position in `points`).
* Note that points in `soup_points` can only be added (intersection points) at the end of the container, with the initial order preserved.
* Note that if `soup_points` contains two or more identical points then only the first copy (following the order in the `soup_points`)
* will be used in `soup_triangles`.
* `soup_triangles` will be updated to contain both the input triangles and the new subdivided triangles. Degenerate triangles will be removed.
* Also triangles in `soup_triangles` will be triangles without intersection first, followed by triangles coming from a subdivision induced
* by an intersection. The named parameter `visitor()` can be used to track
*
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam TriangleRange a model of the concepts `RandomAccessContainer`, `BackInsertionSequence` and `Swappable`, whose
* value type is a model of the concept `RandomAccessContainer` whose value type is convertible to `std::size_t` and that
* is constructible from an `std::initializer_list<std::size_t>` of size 3.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param soup_points points of the soup of polygons
* @param soup_triangles each element in the range describes a triangle using the indexed position of the points in `soup_points`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
*     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
*     \cgalParamDefault{`CGAL::Sequential_tag`}
*   \cgalParamNEnd
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `soup_points`}
*     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
*   \cgalParamNEnd
*   \cgalParamNBegin{visitor}
*     \cgalParamDescription{a visitor used to track the creation of new faces}
*     \cgalParamType{a class model of `PMPAutorefinementVisitor`}
*     \cgalParamDefault{`Autorefinement::Default_visitor`}
*     \cgalParamExtra{The visitor will be copied.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
void autorefine_triangle_soup(PointRange& soup_points,
                              TriangleRange& soup_triangles,
                              const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetPolygonSoupGeomTraits<PointRange, NamedParameters>::type GT;
  typedef typename GetPointMap<PointRange, NamedParameters>::const_type    Point_map;
  Point_map pm = choose_parameter<Point_map>(get_parameter(np, internal_np::point_map));

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  // visitor
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Autorefinement::Default_visitor//default
  > ::type Visitor;
  Visitor visitor(choose_parameter<Visitor>(get_parameter(np, internal_np::visitor)));


  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

#ifndef CGAL_LINKED_WITH_TBB
  static_assert (!parallel_execution,
                 "Parallel_tag is enabled but TBB is unavailable.");
#endif

  typedef std::size_t Input_TID;
  typedef std::pair<Input_TID, Input_TID> Pair_of_triangle_ids;

  // no need for a concurrent vector as the called function itself
  // takes care of sequentially writing into the output iterator
  std::vector<Pair_of_triangle_ids> si_pairs;

  // collect intersecting pairs of triangles
  CGAL_PMP_AUTOREFINE_VERBOSE("collect intersecting pairs");
  triangle_soup_self_intersections<Concurrency_tag>(soup_points, soup_triangles, std::back_inserter(si_pairs), np);

  if (si_pairs.empty())
  {
    if constexpr (!std::is_same_v<Autorefinement::Default_visitor, Visitor>)
    {
      visitor.number_of_output_triangles(soup_triangles.size());
      for(std::size_t i=0; i<soup_triangles.size(); ++i)
        visitor.verbatim_triangle_copy(i, i);
    }
    return;
  }

  // mark degenerate faces so that we can ignore them
  std::vector<bool> is_degen(soup_triangles.size(), false);

  for (const Pair_of_triangle_ids& p : si_pairs)
    if (p.first==p.second) // bbox inter reports (f,f) for degenerate faces
      is_degen[p.first] = true;

  // assign an id per triangle involved in an intersection
  // + the faces involved in the intersection
  std::vector<int> tri_inter_ids(soup_triangles.size(), -1);
  std::vector<Input_TID> intersected_faces;
  int tiid=-1;
  for (const Pair_of_triangle_ids& p : si_pairs)
  {
    if (tri_inter_ids[p.first]==-1 && !is_degen[p.first])
    {
      tri_inter_ids[p.first]=++tiid;
      intersected_faces.push_back(p.first);
    }
    if (tri_inter_ids[p.second]==-1 && !is_degen[p.second])
    {
      tri_inter_ids[p.second]=++tiid;
      intersected_faces.push_back(p.second);
    }
  }

  // init the vector of triangles used for the autorefinement of triangles
  typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
  // even if the info is duplicated with Triangle_data::point, we keep this container
  // so that we can use it in parallel calls of generate_subtriangles
  std::vector< std::array<EK::Point_3, 3> > triangles(tiid+1);
  // vector of data for refining triangles
  std::vector<autorefine_impl::Triangle_data> all_triangle_data(triangles.size());
  Cartesian_converter<GT, EK> to_exact;

  for(Input_TID f : intersected_faces)
  {
    std::size_t tid=tri_inter_ids[f];
    triangles[tid]= CGAL::make_array(
      to_exact( get(pm, soup_points[soup_triangles[f][0]]) ),
      to_exact( get(pm, soup_points[soup_triangles[f][1]]) ),
      to_exact( get(pm, soup_points[soup_triangles[f][2]]) ) );
    all_triangle_data[tid].points.resize(3);
    all_triangle_data[tid].points[0]=triangles[tri_inter_ids[f]][0];
    all_triangle_data[tid].points[1]=triangles[tri_inter_ids[f]][1];
    all_triangle_data[tid].points[2]=triangles[tri_inter_ids[f]][2];
    all_triangle_data[tid].point_locations.resize(3);
    all_triangle_data[tid].point_locations[0]=-1;
    all_triangle_data[tid].point_locations[1]=-2;
    all_triangle_data[tid].point_locations[2]=-3;
  }

  CGAL_PMP_AUTOREFINE_VERBOSE("compute intersections");
#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  Real_timer t;
  t.start();
#endif
  std::set<std::pair<std::size_t, std::size_t> > intersecting_triangles;
  std::set<std::pair<std::size_t, std::size_t> > coplanar_triangles;
  //TODO: PARALLEL_FOR #2
  for (const Pair_of_triangle_ids& p : si_pairs)
  {
    int i1 = tri_inter_ids[p.first],
        i2 = tri_inter_ids[p.second];

    if (i1==-1 || i2==-1) continue; //skip degenerate faces

    const std::array<EK::Point_3, 3>& t1 = triangles[i1];
    const std::array<EK::Point_3, 3>& t2 = triangles[i2];

    std::vector<std::tuple<EK::Point_3, int, int>> inter_pts;
    bool triangles_are_coplanar = autorefine_impl::collect_intersections<EK>(t1, t2, inter_pts);

    CGAL_assertion(
      CGAL::do_intersect(EK::Triangle_3(t1[0], t1[1], t1[2]), EK::Triangle_3(t2[0], t2[1], t2[2]))
      != inter_pts.empty());

    if (!inter_pts.empty())
    {
      std::size_t nbi = inter_pts.size();
      switch(nbi)
      {
        case 1:
          all_triangle_data[i1].add_point<1>(inter_pts[0]);
          all_triangle_data[i2].add_point<2>(inter_pts[0]);
        break;
        case 2:
        {
          std::size_t src_id=all_triangle_data[i1].add_point<1>(inter_pts[0]),
                      tgt_id=all_triangle_data[i1].add_point<1>(inter_pts[1]);
          if (!all_triangle_data[i1].on_same_edge(src_id, tgt_id))
          {
            all_triangle_data[i1].segments.emplace_back(src_id, tgt_id);
            all_triangle_data[i1].segment_input_triangle_ids.push_back(i2);
          }
          //
          src_id=all_triangle_data[i2].add_point<2>(inter_pts[0]);
          tgt_id=all_triangle_data[i2].add_point<2>(inter_pts[1]);
          if (!all_triangle_data[i2].on_same_edge(src_id, tgt_id))
          {
            all_triangle_data[i2].segments.emplace_back(src_id, tgt_id);
            all_triangle_data[i2].segment_input_triangle_ids.push_back(i1);
          }
        }
        break;
        default:
        {
          std::vector<std::size_t> ipt_ids1(nbi+1), ipt_ids2(nbi+1);

          for (std::size_t i=0;i<nbi; ++i)
          {
            ipt_ids1[i]=all_triangle_data[i1].add_point<1>(inter_pts[i]);
            ipt_ids2[i]=all_triangle_data[i2].add_point<2>(inter_pts[i]);
          }
          ipt_ids1.back()=ipt_ids1.front();
          ipt_ids2.back()=ipt_ids2.front();

          for (std::size_t i=0;i<nbi; ++i)
          {
            if (!all_triangle_data[i1].on_same_edge(ipt_ids1[i], ipt_ids1[i+1]))
            {
              all_triangle_data[i1].segments.emplace_back(ipt_ids1[i], ipt_ids1[i+1]);
              all_triangle_data[i1].segment_input_triangle_ids.push_back(i2);
            }
            if (!all_triangle_data[i2].on_same_edge(ipt_ids2[i], ipt_ids2[i+1]))
            {
              all_triangle_data[i2].segments.emplace_back(ipt_ids2[i], ipt_ids2[i+1]);
              all_triangle_data[i2].segment_input_triangle_ids.push_back(i1);
            }
          }
        }
      }
      intersecting_triangles.insert(CGAL::make_sorted_pair(i1, i2));
      if (triangles_are_coplanar)
        coplanar_triangles.insert(CGAL::make_sorted_pair(i1, i2));
    }
  }
#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #2" << std::endl;
  t.reset();
#endif

#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.start();
#endif

  // deduplicate inserted segments
  auto deduplicate_inserted_segments = [&](std::size_t ti)
  {
    if (!all_triangle_data[ti].segments.empty())
    {
      std::vector<EK::Point_3>& points=all_triangle_data[ti].points;
      std::vector<std::pair<std::size_t, std::size_t>>& segments=all_triangle_data[ti].segments;

      std::vector<std::size_t> indices(points.size()-3);
      std::iota(indices.begin(), indices.end(),3);

      std::sort(indices.begin(), indices.end(), [&points](std::size_t i, std::size_t j)
                                                { return points[i]<points[j]; });

      std::vector<std::size_t> id_map(points.size());
      id_map[0]=0;
      id_map[1]=1;
      id_map[2]=2;
      std::vector<std::size_t> unique_ids;
      unique_ids.reserve(indices.size());

      //make points unique + create mapping between indices
      for (std::size_t i=0; i<indices.size(); ++i)
      {
        std::size_t new_id=unique_ids.size()+3;
        unique_ids.push_back(indices[i]);
        id_map[indices[i]]=new_id;
        while(i+1!=indices.size() && points[indices[i]]==points[indices[i+1]])
        {
          id_map[indices[++i]]=new_id;
        }
      }

      CGAL_assertion(points.size() == all_triangle_data[ti].point_locations.size());
      if (unique_ids.size()==indices.size())
        // TODO: do we want to keep points sorted? if yes always swap the 3 containers
        return; // no duplicates

      // now make points unique
      using EPoint_3 = EK::Point_3; // workaround for MSVC 2022 bug
      std::vector<EPoint_3> tmp;
      std::vector<int> tmp_locations;
      tmp.reserve(unique_ids.size()+3);
      tmp_locations.reserve(unique_ids.size()+3);
      tmp.push_back(points[0]);
      tmp.push_back(points[1]);
      tmp.push_back(points[2]);
      tmp_locations.push_back(-1);
      tmp_locations.push_back(-2);
      tmp_locations.push_back(-3);
      for(std::size_t i : unique_ids)
      {
        tmp.push_back(points[i]);
        tmp_locations.push_back(all_triangle_data[ti].point_locations[i]);
      }
      tmp.swap(points);
      tmp_locations.swap(all_triangle_data[ti].point_locations);
      CGAL_assertion(points.size() == all_triangle_data[ti].point_locations.size());

      // now make segments unique
      std::size_t nbs = segments.size();
      std::vector<std::pair<std::size_t, std::size_t>> filtered_segments;
      std::vector<std::size_t> filtered_in_triangle_ids;
      filtered_segments.reserve(nbs);
      filtered_in_triangle_ids.reserve(nbs);
      std::set<std::pair<std::size_t, std::size_t>> segset;
      for (std::size_t si=0; si<nbs; ++si)
      {
        CGAL_assertion(id_map.size()>segments[si].first);
        CGAL_assertion(id_map.size()>segments[si].second);
        std::pair<std::size_t, std::size_t> seg_ids =
          CGAL::make_sorted_pair(id_map[segments[si].first], id_map[segments[si].second]);
        if (segset.insert(seg_ids).second)
        {
          filtered_segments.push_back(seg_ids);
          filtered_in_triangle_ids.push_back(all_triangle_data[ti].segment_input_triangle_ids[si]);
        }
      }
      filtered_in_triangle_ids.swap(all_triangle_data[ti].segment_input_triangle_ids);
      filtered_segments.swap(segments);

      CGAL_assertion(points.size()-3==std::set<EK::Point_3>(std::next(points.begin(),3), points.end()).size());
      CGAL_assertion(points.size()==std::set<EK::Point_3>(points.begin(), points.end()).size());
    }
  };

#ifdef CGAL_LINKED_WITH_TBB
  if (parallel_execution)
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, triangles.size()),
                      [&](const tbb::blocked_range<size_t>& r) {
                        for (size_t ti = r.begin(); ti != r.end(); ++ti)
                          deduplicate_inserted_segments(ti);
                      }
                      );
  }
  else
#endif
  for (std::size_t ti = 0; ti < triangles.size(); ++ti) {
    deduplicate_inserted_segments(ti);
  }

#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #3" << std::endl;
  t.reset();
#endif

  CGAL_PMP_AUTOREFINE_VERBOSE("triangulate faces");
  // now refine triangles
#ifdef CGAL_LINKED_WITH_TBB
  std::conditional_t<parallel_execution,
                     tbb::concurrent_vector<std::pair<std::array<EK::Point_3,3>, std::size_t>>,
                     std::vector<std::pair<std::array<EK::Point_3,3>, std::size_t>>> new_triangles;
#else
  std::vector<std::pair<std::array<EK::Point_3,3>, std::size_t>> new_triangles;
#endif

#ifdef CGAL_AUTOREF_USE_PROGRESS_DISPLAY
  boost::timer::progress_display pd(triangles.size());
#endif

  auto refine_triangles = [&](std::size_t ti)
  {
    if (all_triangle_data[ti].points.size()==3)
      new_triangles.push_back({triangles[ti], ti});
    else
    {
#ifdef CGAL_AUTOREF_USE_FIXED_PROJECTION_TRAITS
      const std::array<typename EK::Point_3, 3>& t = triangles[ti];
      typename EK::Vector_3 orth = CGAL::normal(t[0], t[1], t[2]); // TODO::avoid construction?
      int c = CGAL::abs(orth[0]) > CGAL::abs(orth[1]) ? 0 : 1;
      c = CGAL::abs(orth[2]) > CGAL::abs(orth[c]) ? 2 : c;

      if (c == 0) {
        autorefine_impl::generate_subtriangles<EK, 0>(ti, all_triangle_data[ti], intersecting_triangles, coplanar_triangles, triangles, new_triangles);
      }
      else if (c == 1) {
        autorefine_impl::generate_subtriangles<EK, 1>(ti, all_triangle_data[ti], intersecting_triangles, coplanar_triangles, triangles, new_triangles);
      }
      else if (c == 2) {
        autorefine_impl::generate_subtriangles<EK, 2>(ti, all_triangle_data[ti], intersecting_triangles, coplanar_triangles, triangles, new_triangles);
      }
#else
      autorefine_impl::generate_subtriangles<EK>(ti, all_triangle_data[ti], intersecting_triangles, coplanar_triangles, triangles, new_triangles);
#endif
    }

#ifdef CGAL_AUTOREF_USE_PROGRESS_DISPLAY
    ++pd;
#endif
  };

#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.start();
#endif
#ifdef CGAL_LINKED_WITH_TBB
  if (parallel_execution)
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, triangles.size()),
                      [&](const tbb::blocked_range<size_t>& r) {
                        for (size_t ti = r.begin(); ti != r.end(); ++ti)
                          refine_triangles(ti);
                      }
                      );
  }
  else
#endif
  for (std::size_t ti = 0; ti < triangles.size(); ++ti) {
    refine_triangles(ti);
  }

#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #1" << std::endl;
  t.reset();
#endif

  // brute force output: create a soup, orient and to-mesh
  CGAL_PMP_AUTOREFINE_VERBOSE("create output soup");

  Cartesian_converter<EK, GT> to_input;

#ifdef CGAL_LINKED_WITH_TBB
  typedef std::conditional_t<parallel_execution,
                             tbb::concurrent_map<EK::Point_3, std::size_t>,
                             std::map<EK::Point_3, std::size_t>> Point_id_map;
#else
  typedef std::map<EK::Point_3, std::size_t> Point_id_map;
#endif
  Point_id_map point_id_map;

#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
  std::vector<EK::Point_3> exact_soup_points;
#endif

  // TODO: parallel_for?
  // for input points, we on purpose keep duplicated points and isolated points
  for (std::size_t pid = 0; pid<soup_points.size(); ++pid)
  {
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
    auto insert_res =
#endif
    point_id_map.insert(
      std::make_pair(to_exact(get(pm,soup_points[pid])), pid));
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
      exact_soup_points.push_back(insert_res.first->first);
#endif
  }

  TriangleRange soup_triangles_out;
  soup_triangles_out.reserve(soup_triangles.size());

  if constexpr (!std::is_same_v<Autorefinement::Default_visitor, Visitor>)
  {
    std::size_t nbt=0;
    for (Input_TID f=0; f<soup_triangles.size(); ++f)
    {
      int tiid = tri_inter_ids[f];
      if (tiid == -1) ++nbt;
    }
    visitor.number_of_output_triangles(nbt+new_triangles.size());
  }

  // raw copy of input triangles with no intersection
  std::vector<std::size_t> tri_inter_ids_inverse(triangles.size());
  for (Input_TID f=0; f<soup_triangles.size(); ++f)
  {
    if (is_degen[f]) continue; //skip degenerate faces

    int tiid = tri_inter_ids[f];
    if (tiid == -1)
    {
      visitor.verbatim_triangle_copy(soup_triangles_out.size(), f);
      soup_triangles_out.push_back(
        {soup_triangles[f][0], soup_triangles[f][1], soup_triangles[f][2]}
      );
    }
    else
    {
      tri_inter_ids_inverse[tiid]=f;
    }
  }

  // import refined triangles

  // Lambda to retrieve the id of intersection points
  auto get_point_id = [&](const EK::Point_3& pt)
  {
    auto insert_res = point_id_map.insert(std::make_pair(pt, soup_points.size()));
    if (insert_res.second)
    {
      soup_points.push_back(to_input(pt));
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
      exact_soup_points.push_back(pt);
#endif
    }
    return insert_res.first->second;
  };

#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.start();
#endif

  std::size_t offset = soup_triangles_out.size();
#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  std::string mode = "parallel";
#endif

// It might be possible to optimise the hardcoded value below
// but the less triangles the faster will anyway be the operation.
// So it's probably not critical.
#ifdef CGAL_LINKED_WITH_TBB
  if(parallel_execution && new_triangles.size() > 50)
  {
#ifdef CGAL_AUTOREF_SET_POINT_IDS_USING_MUTEX
    //option 1 (using a mutex)
    CGAL_MUTEX point_container_mutex;
    /// Lambda concurrent_get_point_id()
    auto concurrent_get_point_id = [&](const EK::Point_3& pt)
    {
      auto insert_res = point_id_map.insert(std::make_pair(pt, -1));

      if (insert_res.second)
      {
        CGAL_SCOPED_LOCK(point_container_mutex);
        insert_res.first->second=soup_points.size();
        soup_points.push_back(to_input(pt));
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
        exact_soup_points.push_back(pt);
#endif
      }
      return insert_res.first;
    };

    soup_triangles_out.resize(offset + new_triangles.size());
    //use map iterator triple for triangles to create them concurrently and safely
    std::vector<std::array<typename Point_id_map::iterator, 3>> triangle_buffer(new_triangles.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, new_triangles.size()),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t ti = r.begin(); ti != r.end(); ++ti) {
          const std::array<EK::Point_3, 3>& t = new_triangles[ti].first;
          visitor.new_subtriangle(offset+ti, tri_inter_ids_inverse[new_triangles[ti].second]);
          triangle_buffer[ti] = CGAL::make_array(concurrent_get_point_id(t[0]), concurrent_get_point_id(t[1]), concurrent_get_point_id(t[2]));
        }
      }
    );
    tbb::parallel_for(tbb::blocked_range<size_t>(0, new_triangles.size()),
      [&](const tbb::blocked_range<size_t>& r) {
          for (size_t ti = r.begin(); ti != r.end(); ++ti)
          {
            soup_triangles_out[offset + ti] =
              { triangle_buffer[ti][0]->second,
                triangle_buffer[ti][1]->second,
                triangle_buffer[ti][2]->second };
          }
        }
    );
#else
    //option 2 (without mutex)
    /// Lambda concurrent_get_point_id()
    tbb::concurrent_vector<typename Point_id_map::iterator> iterators;
    auto concurrent_get_point_id = [&](const EK::Point_3& pt)
    {
      auto insert_res = point_id_map.insert(std::make_pair(pt, -1));
      if (insert_res.second)
        iterators.push_back(insert_res.first);
      return insert_res.first;
    };

    //use map iterator triple for triangles to create them concurrently and safely
    soup_triangles_out.resize(offset + new_triangles.size());
    std::vector<std::array<typename Point_id_map::iterator, 3>> triangle_buffer(new_triangles.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, new_triangles.size()),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t ti = r.begin(); ti != r.end(); ++ti)
        {
          const std::array<EK::Point_3, 3>& t = new_triangles[ti].first;
          visitor.new_subtriangle(offset+ti, tri_inter_ids_inverse[new_triangles[ti].second]);
          triangle_buffer[ti] = CGAL::make_array(concurrent_get_point_id(t[0]), concurrent_get_point_id(t[1]), concurrent_get_point_id(t[2]));
        }
      }
    );

    // the map is now filled we can safely set the point ids
    std::size_t pid_offset=soup_points.size();
    soup_points.resize(pid_offset+iterators.size());
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
    exact_soup_points.resize(soup_points.size());
#endif

    tbb::parallel_for(tbb::blocked_range<size_t>(0, iterators.size()),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t ti = r.begin(); ti != r.end(); ++ti)
        {
          soup_points[pid_offset+ti] = to_input(iterators[ti]->first);
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
          exact_soup_points[pid_offset+ti] = iterators[ti]->first;
#endif
          iterators[ti]->second=pid_offset+ti;
        }
      }
    );

    tbb::parallel_for(tbb::blocked_range<size_t>(0, new_triangles.size()),
      [&](const tbb::blocked_range<size_t>& r) {
          for (size_t ti = r.begin(); ti != r.end(); ++ti)
          {
            soup_triangles_out[offset + ti] =
              { triangle_buffer[ti][0]->second,
                triangle_buffer[ti][1]->second,
                triangle_buffer[ti][2]->second };
          }
        }
    );
#endif
  }
  else
#endif
  {
#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
    mode = "sequential";
#endif
    soup_triangles_out.reserve(offset + new_triangles.size());
    for (const std::pair<std::array<EK::Point_3,3>, std::size_t>& t_and_id : new_triangles)
    {
      visitor.new_subtriangle(soup_triangles_out.size(), tri_inter_ids_inverse[t_and_id.second]);
      soup_triangles_out.push_back({ get_point_id(t_and_id.first[0]),
                                     get_point_id(t_and_id.first[1]),
                                     get_point_id(t_and_id.first[2]) });
    }
  }



#ifdef CGAL_AUTOREF_USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #4 (" << mode << ")" << std::endl;
  t.reset();
#endif

#ifndef CGAL_NDEBUG
  CGAL_PMP_AUTOREFINE_VERBOSE("check soup");
  CGAL_assertion( !does_triangle_soup_self_intersect<Concurrency_tag>(exact_soup_points, soup_triangles_out) );
#else
#ifdef CGAL_DEBUG_PMP_AUTOREFINE
  CGAL_PMP_AUTOREFINE_VERBOSE("check soup");
  if (does_triangle_soup_self_intersect<Concurrency_tag>(exact_soup_points, soup_triangles_out))
    throw std::runtime_error("ERROR: invalid output, there is most probably a bug");
#endif
#endif
  using std::swap;
  swap(soup_triangles, soup_triangles_out);

  CGAL_PMP_AUTOREFINE_VERBOSE("done");
}

/**
 * \ingroup PMP_corefinement_grp
 * refines a triangle mesh so that no triangles intersects in their interior.
 *
 * Note that this function is only provided as a shortcut for calling `autorefine_triangle_soup()`
 * with a mesh. For any advance usage the aforementioned function should be called directly.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm input triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * @warning `clear(tm)` will be called before filling `tm` with the refined mesh.
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
void
autorefine(      TriangleMesh& tm,
           const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  // GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  std::vector<typename GT::Point_3> soup_points;
  std::vector<std::array<std::size_t, 3> > soup_triangles;

  polygon_mesh_to_polygon_soup(tm, soup_points, soup_triangles, np);

  autorefine_triangle_soup(soup_points, soup_triangles, np);

  clear(tm);
  repair_polygon_soup(soup_points, soup_triangles);

  duplicate_non_manifold_edges_in_polygon_soup(soup_points, soup_triangles);
  polygon_soup_to_polygon_mesh(soup_points, soup_triangles, tm);
}


} } // end of CGAL::Polygon_mesh_processing

#ifdef CGAL_LINKED_WITH_TBB
#ifdef CGAL_HAS_DEFINED_TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS
#undef TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS
#endif
#endif

#endif // CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
