//TODO: add for soup face the id of the input face. not sure it is easy to report intersection edge as a pair of vertex id
//TODO: only return intersection segments (pay attention to degenerate triangles that are currently ignored)
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

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#ifdef USE_PROGRESS_DISPLAY
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
#include <tbb/concurrent_map.h>
#include <tbb/parallel_for.h>
#ifdef SET_POINT_IDS_USING_MUTEX
#include <CGAL/mutex.h>
#endif
#endif

#include <vector>

#define TEST_RESOLVE_INTERSECTION
#define DEDUPLICATE_SEGMENTS
//#define USE_DEBUG_PARALLEL_TIMERS
//#define DEBUG_COUNTERS
//#define USE_FIXED_PROJECTION_TRAITS
//#define DEBUG_DEPTH

#ifdef USE_DEBUG_PARALLEL_TIMERS
#include <CGAL/Real_timer.h>
#endif

#ifdef USE_FIXED_PROJECTION_TRAITS
#include <CGAL/Kernel_23/internal/Projection_traits_3.h>
#endif

#if defined(DEBUG_COUNTERS) || defined(USE_DEBUG_PARALLEL_TIMERS)
#include <CGAL/Real_timer.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

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
                               const typename K::Vector_3& plane_normal,
                               const K& k = K())
{
  typename K::Collinear_are_ordered_along_line_3 cln_order = k.collinear_are_ordered_along_line_3_object();
#ifdef USE_PROJECTED_ORIENTATION_2_FOR_COPLANAR_ORIENTATION_TESTS
  auto cpl_orient =
    [&plane_normal](const typename K::Point_3& p, const typename K::Point_3& q, const typename K::Point_3& r)
  {
    return ::CGAL::orientation(q-p, r-p, plane_normal);
  };
#else
  typename K::Coplanar_orientation_3 cpl_orient=k.coplanar_orientation_3_object();
  CGAL_USE(plane_normal);
#endif

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

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
template <class Kernel>
void old_intersection_coplanar_triangles_cutoff(const typename Kernel::Point_3& p,
                                                const typename Kernel::Point_3& q,
                                                const typename Kernel::Point_3& r,
                                                const Kernel& k,
                                                std::list<typename Kernel::Point_3>& inter_pts)
{
  typedef typename std::list<typename Kernel::Point_3>::iterator Iterator;

  if(inter_pts.empty())
    return;

  typename Kernel::Coplanar_orientation_3 orient = k.coplanar_orientation_3_object();
  typename Kernel::Construct_line_3 line = k.construct_line_3_object();

  //orient(p,q,r,r) is POSITIVE
  std::map<const typename Kernel::Point_3*,Orientation> orientations;
  for (Iterator it=inter_pts.begin();it!=inter_pts.end();++it)
    orientations[ &(*it) ]=orient(p,q,r,*it);

  CGAL_kernel_assertion_code(int pt_added = 0;)

  const typename Kernel::Point_3* prev = &(*boost::prior(inter_pts.end()));
  Iterator stop = inter_pts.size() > 2 ? inter_pts.end() : boost::prior(inter_pts.end());
  for(Iterator it=inter_pts.begin(); it!=stop; ++it)
  {
    const typename Kernel::Point_3& curr = *it;
    Orientation or_prev = orientations[prev],
                or_curr = orientations[&curr];

    if((or_prev == POSITIVE && or_curr == NEGATIVE) ||
       (or_prev == NEGATIVE && or_curr == POSITIVE))
    {
      typename Intersection_traits<Kernel, typename Kernel::Line_3, typename Kernel::Line_3>::result_type
        obj = ::CGAL::Intersections::internal::intersection(line(p,q), line(*prev,curr), k);

      // assert "not empty"
      CGAL_kernel_assertion(bool(obj));

      const typename Kernel::Point_3* inter = ::CGAL::Intersections::internal::intersect_get<typename Kernel::Point_3>(obj);
      CGAL_kernel_assertion(inter != nullptr);

      prev = &(*inter_pts.insert(it,*inter));
      orientations[prev] = COLLINEAR;
      CGAL_kernel_assertion_code(++pt_added;)
    }

    prev = &(*it);
  }

  CGAL_kernel_assertion(pt_added<3);
  Iterator it = inter_pts.begin();
  while(it!=inter_pts.end())
  {
    if(orientations[&(*it)] == NEGATIVE)
      inter_pts.erase(it++);
    else
      ++it;
  }
}
#endif

// imported from Intersections_3/include/CGAL/Intersections_3/internal/Triangle_3_Triangle_3_intersection.h
template<class K>
void coplanar_intersections(const std::array<typename K::Point_3, 3>& t1,
                            const std::array<typename K::Point_3, 3>& t2,
                            std::vector<typename K::Point_3>& inter_pts)
{
  const typename K::Point_3& p1 = t1[0], q1 = t1[1], r1 = t1[2];
  const typename K::Point_3& p2 = t2[0], q2 = t2[1], r2 = t2[2];

  std::list<Intersections::internal::Point_on_triangle<K>> l_inter_pts;
  l_inter_pts.push_back(Intersections::internal::Point_on_triangle<K>(-1,0));
  l_inter_pts.push_back(Intersections::internal::Point_on_triangle<K>(-1,1));
  l_inter_pts.push_back(Intersections::internal::Point_on_triangle<K>(-1,2));

#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::list<typename K::Point_3> old_l_inter_pts;
  old_l_inter_pts.push_back(p2);
  old_l_inter_pts.push_back(q2);
  old_l_inter_pts.push_back(r2);


  auto enum_to_string = [](CGAL::Orientation o)
  {
    if (o==COLLINEAR) return std::string("COLLINEAR");
    if (o==POSITIVE) return std::string("POSITIVE");
    return std::string("NEGATIVE");
  };

  auto to_string = [](const auto& t)
  {
    std::stringstream sstr;
    sstr << "4 "
         << to_double(t[0].x()) << " "  << to_double(t[0].y()) << " "  << to_double(t[0].z()) << " "
         << to_double(t[1].x()) << " "  << to_double(t[1].y()) << " "  << to_double(t[1].z()) << " "
         << to_double(t[2].x()) << " "  << to_double(t[2].y()) << " "  << to_double(t[2].z()) << " "
         << to_double(t[0].x()) << " "  << to_double(t[0].y()) << " "  << to_double(t[0].z()) << "\n";
    return sstr.str();
  };

  std::cout << "intersection_coplanar_triangles\n";
  std::ofstream("/tmp/t1.polylines.txt") << to_string(t1) << "\n";
  std::ofstream("/tmp/t2.polylines.txt") << to_string(t2) << "\n";

  std::cout << "Position of vertices of t1: ";
  std::cout << enum_to_string( coplanar_orientation(p2,q2,r2,p1)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(p2,q2,r2,q1)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(p2,q2,r2,r1)) << "\n";
  std::cout << "                            ";
  std::cout << enum_to_string( coplanar_orientation(q2,r2,p2,p1)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(q2,r2,p2,q1)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(q2,r2,p2,r1)) << "\n";
  std::cout << "                            ";
  std::cout << enum_to_string( coplanar_orientation(r2,p2,q2,p1)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(r2,p2,q2,q1)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(r2,p2,q2,r1)) << "\n";
  std::cout << "Position of vertices of t2: ";
  std::cout << enum_to_string( coplanar_orientation(p1,q1,r1,p2)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(p1,q1,r1,q2)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(p1,q1,r1,r2)) << "\n";
  std::cout << "                            ";
  std::cout << enum_to_string( coplanar_orientation(q1,r1,p1,p2)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(q1,r1,p1,q2)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(q1,r1,p1,r2)) << "\n";
  std::cout << "                            ";
  std::cout << enum_to_string( coplanar_orientation(r1,p1,q1,p2)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(r1,p1,q1,q2)) << " - ";
  std::cout << enum_to_string( coplanar_orientation(r1,p1,q1,r2)) << "\n";

  auto print_points = [&]()
  {
    for(auto p : l_inter_pts) std::cout << "  (" << p.id1() << "," << p.id2() << ",[" << p.alpha << "]) "; std::cout <<"\n";
  };
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
#endif

  //intersect t2 with the three half planes which intersection defines t1
  K k;
  Intersections::internal::intersection_coplanar_triangles_cutoff(p1,q1,r1,0,p2,q2,r2,k,l_inter_pts); //line p1q1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
  old_intersection_coplanar_triangles_cutoff(p1,q1,r1,k,old_l_inter_pts);
  CGAL_assertion(l_inter_pts.size()==old_l_inter_pts.size());
  for (std::size_t i=0; i<l_inter_pts.size(); ++i)
  {
    if (*(std::next(old_l_inter_pts.begin(),i))!=std::next(l_inter_pts.begin(),i)->point(p1,q1,r1,p2,q2,r2,k))
    {
      std::cout <<"ERROR with point #" << i << "\n";
      throw std::runtime_error("invalid output 0");
    }
  }
#endif
  Intersections::internal::intersection_coplanar_triangles_cutoff(q1,r1,p1,1,p2,q2,r2,k,l_inter_pts); //line q1r1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
  old_intersection_coplanar_triangles_cutoff(q1,r1,p1,k,old_l_inter_pts);
  CGAL_assertion(l_inter_pts.size()==old_l_inter_pts.size());
  for (std::size_t i=0; i<l_inter_pts.size(); ++i)
  {
    if (*(std::next(old_l_inter_pts.begin(),i))!=std::next(l_inter_pts.begin(),i)->point(p1,q1,r1,p2,q2,r2,k))
    {
      std::cout <<"ERROR with point #" << i << "\n";
      throw std::runtime_error("invalid output 1");
    }
  }
#endif
  Intersections::internal::intersection_coplanar_triangles_cutoff(r1,p1,q1,2,p2,q2,r2,k,l_inter_pts); //line r1p1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
  old_intersection_coplanar_triangles_cutoff(r1,p1,q1,k,old_l_inter_pts);
  CGAL_assertion(l_inter_pts.size()==old_l_inter_pts.size());
  for (std::size_t i=0; i<l_inter_pts.size(); ++i)
  {
    if (*(std::next(old_l_inter_pts.begin(),i))!=std::next(l_inter_pts.begin(),i)->point(p1,q1,r1,p2,q2,r2,k))
    {
      std::cout <<"ERROR with point #" << i << "\n";
      throw std::runtime_error("invalid output 2");
    }
  }
#endif

#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::size_t start=inter_pts.size();
#endif

  for (const Intersections::internal::Point_on_triangle<K>& pot : l_inter_pts)
    inter_pts.push_back( pot.point(p1,q1,r1,p2,q2,r2,k) );

#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  for (std::size_t i=0; i<l_inter_pts.size(); ++i)
  {
    typename K::Point_3 inter = inter_pts[start+i];
    CGAL_assertion( coplanar_orientation(p1,q1,r1,inter)!=NEGATIVE );
    CGAL_assertion( coplanar_orientation(q1,r1,p1,inter)!=NEGATIVE );
    CGAL_assertion( coplanar_orientation(r1,p1,q1,inter)!=NEGATIVE );
    CGAL_assertion( coplanar_orientation(p2,q2,r2,inter)!=NEGATIVE );
    CGAL_assertion( coplanar_orientation(q2,r2,p2,inter)!=NEGATIVE );
    CGAL_assertion( coplanar_orientation(r2,p2,q2,inter)!=NEGATIVE );
  }

  std::ofstream debug("interpts.xyz");
  debug << std::setprecision(17);
  debug << l_inter_pts.size() << "\n";
  for (auto pot : l_inter_pts)
    debug << pot.point(p1,q1,r1,p2,q2,r2,k)  << "\n";
  debug.close();
  // std::cout <<"check!\n";
  // int i;
  // std::cin >> i;
#endif

}

// imported from Polygon_mesh_processing/internal/Corefinement/intersect_triangle_segment_3.h
template<class K>
void
find_intersection(const typename K::Point_3& p, const typename K::Point_3& q,  //segment
                  const typename K::Point_3& a, const typename K::Point_3& b, const typename K::Point_3& c, //triangle
                  std::vector<typename K::Point_3>& inter_pts,
                  bool is_p_coplanar=false, bool is_q_coplanar=false) // note that in coref this was wrt a halfedge not p/q
{
  Orientation ab=orientation(p,q,a,b);
  Orientation bc=orientation(p,q,b,c);
  Orientation ca=orientation(p,q,c,a);

  if ( ab==POSITIVE || bc==POSITIVE || ca==POSITIVE )
    return;

  int nb_coplanar=(ab==COPLANAR?1:0) + (bc==COPLANAR?1:0) + (ca==COPLANAR?1:0);

/*
  if ( nb_coplanar==0 )
    return result_type(ON_FACE,hd,is_src_coplanar,is_tgt_coplanar);

  if (nb_coplanar==1){
    if (ab==COPLANAR)
      // intersection is ab
      return result_type(ON_EDGE,next(hd,tm),is_src_coplanar,is_tgt_coplanar);
    if (bc==COPLANAR)
      // intersection is bc
      return result_type(ON_EDGE,prev(hd,tm),is_src_coplanar,is_tgt_coplanar);
    CGAL_assertion(ca==COPLANAR);
    // intersection is ca
    return result_type(ON_EDGE,hd,is_src_coplanar,is_tgt_coplanar);
  }
*/

  if (is_p_coplanar)
  {
    inter_pts.push_back(p);
    return;
  }
  if (is_q_coplanar)
  {
    inter_pts.push_back(q);
    return;
  }

  if (nb_coplanar!=2)
  {
    inter_pts.push_back(
      typename K::Construct_plane_line_intersection_point_3()(a, b, c, p, q)
    );
  }
  else
  {
    if (ab!=COPLANAR)
    {
      // intersection is c
      inter_pts.push_back(c);
      return;
    }

    if (bc!=COPLANAR)
    {
      // intersection is a
      inter_pts.push_back(a);
      return;
    }
    CGAL_assertion(ca!=COPLANAR);
    // intersection is b
    inter_pts.push_back(b);
  }
}

template <class K>
void test_edge(const typename K::Point_3& p, const typename K::Point_3& q,
               const typename K::Point_3& a, const typename K::Point_3& b, const typename K::Point_3& c,
               const Orientation abcp,
               const Orientation abcq,
               std::vector<typename K::Point_3>& inter_pts)
{
  switch ( abcp ) {
  case POSITIVE:
    switch ( abcq ) {
    case POSITIVE:
      // the segment lies in the positive open halfspaces defined by the
      // triangle's supporting plane
    break;
    case NEGATIVE:
      // p sees the triangle in counterclockwise order
      find_intersection<K>(p,q,a,b,c,inter_pts);
    break;
    //case COPLANAR:
    default:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in counterclockwise order
      find_intersection<K>(p,q,a,b,c,inter_pts,false,true);
    }
  break;
  case NEGATIVE:
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      find_intersection<K>(q,p,a,b,c,inter_pts);
    break;
    case NEGATIVE:
      // the segment lies in the negative open halfspaces defined by the
      // triangle's supporting plane
    break;
    // case COPLANAR:
    default:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in clockwise order
      find_intersection<K>(q,p,a,b,c,inter_pts,true,false);
    }
  break;
  default:
  //case COPLANAR: // p belongs to the triangle's supporting plane
    switch ( abcq ) {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      find_intersection<K>(q,p,a,b,c,inter_pts,false, true);
    break;
    case NEGATIVE:
      // q sees the triangle in clockwise order
      find_intersection<K>(p,q,a,b,c,inter_pts,true);
    break;
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
    break;
    }
  }
}

template <class K>
bool collect_intersections(const std::array<typename K::Point_3, 3>& t1,
                           const std::array<typename K::Point_3, 3>& t2,
                           std::vector<typename K::Point_3>& inter_pts)
{
  // test edges of t1 vs t2
  std::array<Orientation,3> ori;
  for (int i=0; i<3; ++i)
    ori[i] = orientation(t2[0],t2[1],t2[2],t1[i]);

  if (ori[0]== COPLANAR && ori[1]==COPLANAR && ori[2]==COPLANAR)
  {
    coplanar_intersections<K>(t1, t2, inter_pts);
#ifdef DEBUG_DEPTH
    for (auto p : inter_pts)
      if (depth(p)>2) throw std::runtime_error("Depth is not 4: "+std::to_string(depth(p)));
#endif

    return true;
  }

  for (int i=0; i<3; ++i)
  {
    int j=(i+1)%3;
    test_edge<K>(t1[i], t1[j], t2[0], t2[1], t2[2], ori[i], ori[j], inter_pts);
  }

  // test edges of t2 vs t1
  for (int i=0; i<3; ++i)
    ori[i] = orientation(t1[0],t1[1],t1[2],t2[i]);
  for (int i=0; i<3; ++i)
  {
    int j=(i+1)%3;
    test_edge<K>(t2[i], t2[j], t1[0], t1[1], t1[2], ori[i], ori[j], inter_pts);
  }

  // because we don't handle intersection type and can have edge-edge edge-vertex duplicates
  std::sort(inter_pts.begin(), inter_pts.end());
  auto last = std::unique(inter_pts.begin(), inter_pts.end());
  inter_pts.erase(last, inter_pts.end());

#ifdef DEBUG_DEPTH
  for (auto p : inter_pts)
    if (depth(p)>2) throw std::runtime_error("Depth is not 2: "+std::to_string(depth(p)));
#endif

  return false;
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <class EK,
#ifdef USE_FIXED_PROJECTION_TRAITS
          int dim,
#endif
          class PointVector
>
void generate_subtriangles(std::size_t ti,
                           std::vector<std::pair<std::size_t, std::size_t>>& segments,
                           std::vector<typename EK::Point_3>& points,
                           const std::vector<std::size_t>& in_triangle_ids,
                           const std::set<std::pair<std::size_t, std::size_t> >& intersecting_triangles,
                           const std::set<std::pair<std::size_t, std::size_t> >& coplanar_triangles,
                           const std::vector<std::array<typename EK::Point_3,3>>& triangles,
                           PointVector& new_triangles
                           )
{
  // std::cout << "generate_subtriangles()\n";
  // std::cout << std::setprecision(17);

#ifdef USE_FIXED_PROJECTION_TRAITS
  typedef ::CGAL::internal::Projection_traits_3<EK, dim> P_traits;
#else
  typedef CGAL::Projection_traits_3<EK> P_traits;
#endif

#ifndef TEST_RESOLVE_INTERSECTION
  typedef CGAL::Exact_intersections_tag Itag;
#else
  typedef CGAL::No_constraint_intersection_tag Itag;
#endif

  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits ,Default, Itag> CDT_2;
  //typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT;
  typedef CDT_2 CDT;

  const std::array<typename EK::Point_3,3>& t = triangles[ti];

  // positive triangle normal
  typename EK::Vector_3 n = normal(t[0], t[1], t[2]);
  typename EK::Point_3 o(CGAL::ORIGIN);

#ifdef USE_FIXED_PROJECTION_TRAITS
  P_traits cdt_traits;
  bool orientation_flipped = false;
  CDT cdt(cdt_traits);
  // TODO: still need to figure out why I can't make the orientation_flipped correctly
  cdt.insert(t[0]);
  cdt.insert(t[1]);
  cdt.insert(t[2]);
#else
  bool orientation_flipped = false;
  if ( typename EK::Less_xyz_3()(o+n,o) )
  {
    n=-n;
    orientation_flipped = true;
  }

  P_traits cdt_traits(n);
  CDT cdt(cdt_traits);
  cdt.insert_outside_affine_hull(t[0]);
  cdt.insert_outside_affine_hull(t[1]);
  typename CDT::Vertex_handle v = cdt.tds().insert_dim_up(cdt.infinite_vertex(), orientation_flipped);
  v->set_point(t[2]);
#endif

#ifdef DEBUG_COUNTERS
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
#define COUNTER_INSTRUCTION(X) X
#else
#define COUNTER_INSTRUCTION(X)
#endif


#ifdef TEST_RESOLVE_INTERSECTION
  //~ static std::ofstream debug("inter_segments.polylines.txt");
  //~ debug.precision(17);

  //~ std::cout << "points.size() " << points.size() << "\n";
  //~ std::set<std::size_t> all_triangles_indices(in_triangle_ids.begin(), in_triangle_ids.end());
  //~ all_triangles_indices.insert(ti);

  //~ std::ofstream debug("triangles.polylines.txt");
  //~ debug << std::setprecision(17);
  //~ for (std::size_t i : all_triangles_indices)
    //~ debug << "4 "
          //~ << triangles[i][0] << " "
          //~ << triangles[i][1] << " "
          //~ << triangles[i][2] << " "
          //~ << triangles[i][0] << "\n";
  //~ debug.close();
  //~ debug.open("triangle.off");
  //~ debug << std::setprecision(17);
  //~ debug << "OFF\n3 1 0\n";
  //~ debug << triangles[ti][0] << "\n"
        //~ << triangles[ti][1] << "\n"
        //~ << triangles[ti][2] << "\n 3 0 1 2\n";
  //~ debug.close();

  // pre-compute segment intersections
  if (!segments.empty())
  {
    std::size_t nbs = segments.size();

    std::vector< std::vector<std::size_t> > points_on_segments(nbs);

    COUNTER_INSTRUCTION(counter.timer1.start();)

    std::map<typename EK::Point_3, std::size_t> point_id_map;

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
        points.push_back(pt);
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
            COUNTER_INSTRUCTION(++counter.c5;)
            COUNTER_INSTRUCTION(++counter.total;)
            continue;
          }

          COUNTER_INSTRUCTION(counter.timer5.start();)
          Segment_inter_type seg_inter_type =
            do_coplanar_segments_intersect<EK>(segments[i].first, segments[i].second,
                                               segments[j].first, segments[j].second,
                                               points, n);
          COUNTER_INSTRUCTION(counter.timer5.stop();)

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
                  &&  coplanar_triangles.count(CGAL::make_sorted_pair(ti, in_triangle_ids[j])) == 0
                  &&  coplanar_triangles.count(CGAL::make_sorted_pair(in_triangle_ids[i], ti)) == 0)
              {
                COUNTER_INSTRUCTION(counter.timer6.start();)
                typename EK::Point_3 pt = typename EK::Construct_planes_intersection_point_3()(
                  triangles[in_triangle_ids[i]][0], triangles[in_triangle_ids[i]][1],triangles[in_triangle_ids[i]][2],
                  triangles[in_triangle_ids[j]][0], triangles[in_triangle_ids[j]][1],triangles[in_triangle_ids[j]][2],
                  triangles[ti][0], triangles[ti][1],triangles[ti][2]);
                COUNTER_INSTRUCTION(counter.timer6.stop();)

                COUNTER_INSTRUCTION(++counter.c1;)
                std::size_t pid = get_point_id(pt);
                points_on_segments[i].push_back(pid);
                points_on_segments[j].push_back(pid);
              }
              else
              {
                COUNTER_INSTRUCTION(++counter.c2;)
                COUNTER_INSTRUCTION(counter.timer4.start();)
                typename EK::Point_3 pt = typename EK::Construct_coplanar_segments_intersection_point_3()(
                  points[segments[i].first], points[segments[i].second],
                  points[segments[j].first], points[segments[j].second]);

                std::size_t pid = get_point_id(pt);
                points_on_segments[i].push_back(pid);
                points_on_segments[j].push_back(pid);
                COUNTER_INSTRUCTION(counter.timer4.stop();)
                //~ std::ofstream debug ("/tmp/triangles.polylines.txt");
                //~ debug << "4 " << triangles[ti][0] << " " << triangles[ti][1] << " " << triangles[ti][2] << " " << triangles[ti][0] << "\n";
                //~ debug << "4 " << triangles[in_triangle_ids[i]][0] << " " << triangles[in_triangle_ids[i]][1] << " " << triangles[in_triangle_ids[i]][2] << " " << triangles[in_triangle_ids[i]][0] << "\n";
                //~ debug << "4 " << triangles[in_triangle_ids[j]][0] << " " << triangles[in_triangle_ids[j]][1] << " " << triangles[in_triangle_ids[j]][2] << " " << triangles[in_triangle_ids[j]][0] << "\n";
                //~ debug.close();
                //~ throw std::runtime_error("Unexpected case 1");
              }
              break;
            }
            case COPLANAR_SEGMENT_PQ:
            {
              COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[j].push_back(segments[i].first);
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case COPLANAR_SEGMENT_RS:
            {
              COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].first);
              points_on_segments[i].push_back(segments[j].second);
              break;
            }
            case COPLANAR_SEGMENT_PR:
            {
              COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].first);
              points_on_segments[j].push_back(segments[i].first);
              break;
            }
            case COPLANAR_SEGMENT_QS:
            {
              COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].second);
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case COPLANAR_SEGMENT_PS:
            {
              COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].second);
              points_on_segments[j].push_back(segments[i].first);
              break;
            }
            case COPLANAR_SEGMENT_QR:
            {
              COUNTER_INSTRUCTION(++counter.c3;)
              points_on_segments[i].push_back(segments[j].first);
              points_on_segments[j].push_back(segments[i].second);
              break;
            }
            case NO_INTERSECTION:
            {
              COUNTER_INSTRUCTION(++counter.c4;)
            }
          }
        }
        COUNTER_INSTRUCTION(++counter.total;)
      }
    }
    COUNTER_INSTRUCTION(counter.timer1.stop();)
    COUNTER_INSTRUCTION(counter.timer2.start();)
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

        //~ {
        //~ std::cout << "cst_points.size() " << cst_points.size() << "\n";
        //~ std::ofstream debugbis("subcst.polylines.txt");
        //~ debugbis << std::setprecision(17);

        //~ std::cout << "splittng " << segments[i] << "\n";
        //~ static int kkk=-1;
        //~ ++kkk;
        //~ std::ofstream debugter("subcst_"+std::to_string(i)+".polylines.txt");
        //~ for (std::size_t k=0; k<=points_on_segments[i].size(); ++k)
        //~ {
          //~ exact(cst_points[src_id+k]);
          //~ exact(cst_points[src_id+k+1]);
          //~ debugbis << "2 " << cst_points[src_id+k] << " " <<  cst_points[src_id+k+1] << "\n";
          //~ debugter << "2 " << cst_points[src_id+k] << " " <<  cst_points[src_id+k+1] << "\n";
        //~ }
        //~ std::cout << "---\n";
        //~ }
      }
    }
    COUNTER_INSTRUCTION(counter.timer2.stop();)

    //~ int max_degree = 0;
    //~ for (const auto p : cst_points)
      //~ max_degree = std::max(max_degree, depth(p));
    //~ std::cout << "max_degree " << max_degree << "\n";

    //~ if (max_degree > 10){
      //~ for (const auto p : cst_points)
        //~ std::cout << " -- " << p << "(" << depth(p) << ")\n";
      //~ std::cout << "segments:\n";
      //~ for (auto s : segments)
        //~ std::cout << "  " << depth(s[0]) << " " << depth(s[1]) << "\n";
      //~ exit(1);
    //~ }

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

  // TODO: sorted pair to be constructed when pushing_back
  // TODO: only needed in case of coplanar segments?
  for (std::pair<std::size_t, size_t>& s : segments)
    if (s.second < s.first)
      std::swap(s.first,s.second);
  std::sort(segments.begin(), segments.end());
  auto last = std::unique(segments.begin(), segments.end());
  segments.erase(last, segments.end());
#endif

  //~ std::ofstream("/tmp/tri.xyz") << std::setprecision(17) << triangles[ti][0] << "\n"
                                                         //~ << triangles[ti][1] << "\n"
                                                         //~ << triangles[ti][2] << "\n";
  //~ std::ofstream debug("/tmp/cst.polylines.txt");
  //~ debug << std::setprecision(17);
  //~ for(auto s : segments)
    //~ debug << "2 " << points[s.first] << " " << points[s.second] << "\n";
  //~ debug.close();

  COUNTER_INSTRUCTION(counter.timer3.start();)
  if (segments.empty())
    cdt.insert(points.begin(), points.end());
  else
    cdt.insert_constraints(points.begin(), points.end(), segments.begin(), segments.end());
  COUNTER_INSTRUCTION(counter.timer3.stop();)

#ifdef CGAL_DEBUG_PMP_AUTOREFINE_DUMP_TRIANGULATIONS
    static int k = 0;
    std::stringstream buffer;
    buffer.precision(17);
    int nbt=0;
#endif
    for (typename CDT::Face_handle fh : cdt.finite_face_handles())
    {
      if (orientation_flipped)
        new_triangles.push_back( CGAL::make_array(fh->vertex(0)->point(),
                                                  fh->vertex(cdt.cw(0))->point(),
                                                  fh->vertex(cdt.ccw(0))->point()) );
      else
        new_triangles.push_back( CGAL::make_array(fh->vertex(0)->point(),
                                                  fh->vertex(cdt.ccw(0))->point(),
                                                  fh->vertex(cdt.cw(0))->point()) );
#ifdef CGAL_DEBUG_PMP_AUTOREFINE_DUMP_TRIANGULATIONS
      ++nbt;
      buffer << fh->vertex(0)->point() << "\n";
      buffer << fh->vertex(cdt.ccw(0))->point() << "\n";
      buffer << fh->vertex(cdt.cw(0))->point() << "\n";
#endif
    }

#ifdef CGAL_DEBUG_PMP_AUTOREFINE_DUMP_TRIANGULATIONS
    std::ofstream dump("triangulation_"+std::to_string(k)+".off");
    dump << "OFF\n" << 3*nbt << " " << nbt << " 0\n";
    dump << buffer.str();
    for (int i=0; i<nbt; ++i)
      dump << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
    ++k;
#endif
}

} // end of autorefine_impl
#endif

/**
* \ingroup PMP_corefinement_grp
*
* refines a soup of triangles so that no pair of triangles intersects in their interior.
* Note that points in `soup_points` can only be added (intersection points) a the end of the container, with the initial order preserved.
* Note that if `soup_points` contains two or more identical points and only the first copy (following the order in the `soup_points`)
* will be used in `soup_triangles`.
* `soup_triangles` will be updated to contain both the input triangles and the new subdivides triangles. Degenerate triangles will be removed.
* Also triangles in `soup_triangles` will be triangles without intersection first, followed by triangles coming from a subdivision induced
* by an intersection. The named parameter `visitor()` can be used to track
*
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam TriIdsRange a model of the concepts `RandomAccessContainer`, `BackInsertionSequence` and `Swappable`, whose
* value type is a model of the concept `RandomAccessContainer` whose value type is convertible to `std::size_t` and that
* is constructible from an `std::initializer_list<std::size_t>` of size 3.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param soup_points points of the soup of polygons
* @param soup_triangles each element in the range describes a triangle using the indexed position of the points in `soup_points`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
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
*     \cgalParamType{a class model of `PMPFooBar`}
*     \cgalParamDefault{`Autorefinement::Default_visitor<Bar, Foo>`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class PointRange, class TriIdsRange, class NamedParameters = parameters::Default_named_parameters>
void autorefine_triangle_soup(PointRange& soup_points,
                              TriIdsRange& soup_triangles,
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

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

#ifndef CGAL_LINKED_WITH_TBB
  static_assert (!parallel_execution,
                 "Parallel_tag is enabled but TBB is unavailable.");
#endif

  typedef std::size_t Input_TID;
  typedef std::pair<Input_TID, Input_TID> Pair_of_triangle_ids;

  std::vector<Pair_of_triangle_ids> si_pairs; // TODO: check std::vector is fine with Parallel_tag

  // collect intersecting pairs of triangles
  CGAL_PMP_AUTOREFINE_VERBOSE("collect intersecting pairs");
  triangle_soup_self_intersections<Concurrency_tag>(soup_points, soup_triangles, std::back_inserter(si_pairs), np);

  if (si_pairs.empty()) return;

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
  std::vector< std::array<typename EK::Point_3, 3> > triangles(tiid+1);
  Cartesian_converter<GT, EK> to_exact;

  for(Input_TID f : intersected_faces)
  {
    triangles[tri_inter_ids[f]]= CGAL::make_array(
      to_exact( get(pm, soup_points[soup_triangles[f][0]]) ),
      to_exact( get(pm, soup_points[soup_triangles[f][1]]) ),
      to_exact( get(pm, soup_points[soup_triangles[f][2]]) ) );
  }

  std::vector< std::vector<std::array<EK::Point_3,2> > > all_segments(triangles.size());
  std::vector< std::vector<EK::Point_3> > all_points(triangles.size());
  std::vector< std::vector<std::size_t> > all_in_triangle_ids(triangles.size());

  CGAL_PMP_AUTOREFINE_VERBOSE("compute intersections");
#ifdef USE_DEBUG_PARALLEL_TIMERS
  Real_timer t;
  t.start();
#endif
  std::set<std::pair<std::size_t, std::size_t> > intersecting_triangles; // TODO replace with vector<boost::flat_set<>>>
  std::set<std::pair<std::size_t, std::size_t> > coplanar_triangles; // TODO replace with vector<boost::flat_set<>>>
  //TODO: PARALLEL_FOR #2
  for (const Pair_of_triangle_ids& p : si_pairs)
  {
    int i1 = tri_inter_ids[p.first],
        i2 = tri_inter_ids[p.second];

    if (i1==-1 || i2==-1) continue; //skip degenerate faces

    const std::array<typename EK::Point_3, 3>& t1 = triangles[i1];
    const std::array<typename EK::Point_3, 3>& t2 = triangles[i2];

    std::vector<typename EK::Point_3> inter_pts;
    bool triangles_are_coplanar = autorefine_impl::collect_intersections<EK>(t1, t2, inter_pts);

    CGAL_assertion(
      CGAL::do_intersect(typename EK::Triangle_3(t1[0], t1[1], t1[2]), typename EK::Triangle_3(t2[0], t2[1], t2[2]))
      != inter_pts.empty());

    if (!inter_pts.empty())
    {
      std::size_t nbi = inter_pts.size();
      switch(nbi)
      {
        case 1:
          all_points[i1].push_back(inter_pts[0]);
          all_points[i2].push_back(inter_pts[0]);
        break;
        case 2:
          all_segments[i1].push_back(CGAL::make_array(inter_pts[0], inter_pts[1]));
          all_segments[i2].push_back(CGAL::make_array(inter_pts[0], inter_pts[1]));
          all_in_triangle_ids[i1].push_back(i2);
          all_in_triangle_ids[i2].push_back(i1);
        break;
        default:
          for (std::size_t i=0;i<nbi; ++i)
          {
            all_segments[i1].push_back(CGAL::make_array(inter_pts[i], inter_pts[(i+1)%nbi]));
            all_segments[i2].push_back(all_segments[i1].back());
            all_in_triangle_ids[i1].push_back(i2);
            all_in_triangle_ids[i2].push_back(i1);
          }
      }
      intersecting_triangles.insert(CGAL::make_sorted_pair(i1, i2));
      if (triangles_are_coplanar)
        coplanar_triangles.insert(CGAL::make_sorted_pair(i1, i2));
    }
  }
#ifdef USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #2" << std::endl;
  t.reset();
#endif

#ifdef DEDUPLICATE_SEGMENTS
#ifdef USE_DEBUG_PARALLEL_TIMERS
  t.start();
#endif

  // deduplicate inserted segments
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> all_segments_ids(all_segments.size());

  auto deduplicate_inserted_segments = [&](std::size_t ti)
  {
    if (!all_segments[ti].empty())
    {
      std::map<EK::Point_3, std::size_t> point_id_map;


      auto get_point_id = [&](const typename EK::Point_3& pt)
      {
        auto insert_res = point_id_map.insert(std::make_pair(pt, all_points[ti].size()));
        if (insert_res.second)
          all_points[ti].push_back(pt);
        return insert_res.first->second;
      };


      if (!all_points[ti].empty())
      {
        std::vector<typename EK::Point_3> tmp;
        tmp.swap(all_points[ti]);
        for (const typename EK::Point_3& pt : tmp)
          get_point_id(pt);
      }

      std::size_t nbs = all_segments[ti].size();
      std::vector<std::array<EK::Point_3,2>> filtered_segments;
      std::vector<std::size_t> filtered_in_triangle_ids;
      filtered_segments.reserve(nbs);
      std::set<std::pair<std::size_t, std::size_t>> segset;
      for (std::size_t si=0; si<nbs; ++si)
      {
        EK::Point_3 src = all_segments[ti][si][0],
                    tgt = all_segments[ti][si][1];
        std::size_t src_id = get_point_id(src), tgt_id=get_point_id(tgt);
        if (segset.insert(
              CGAL::make_sorted_pair(src_id, tgt_id)).second)
        {
          all_segments_ids[ti].emplace_back(src_id, tgt_id);
          filtered_in_triangle_ids.push_back(all_in_triangle_ids[ti][si]);
        }
      }
      if (all_segments_ids[ti].size()!=nbs)
        filtered_in_triangle_ids.swap(all_in_triangle_ids[ti]);
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

#ifdef USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #3" << std::endl;
  t.reset();
#endif
#endif

  CGAL_PMP_AUTOREFINE_VERBOSE("triangulate faces");
  // now refine triangles
#ifdef CGAL_LINKED_WITH_TBB
  std::conditional_t<parallel_execution,
                     tbb::concurrent_vector<std::array<EK::Point_3, 3>>,
                     std::vector<std::array<EK::Point_3,3>>> new_triangles;
#else
  std::vector<std::array<EK::Point_3,3>> new_triangles;
#endif

#ifdef USE_PROGRESS_DISPLAY
  boost::timer::progress_display pd(triangles.size());
#endif

  auto refine_triangles = [&](std::size_t ti)
  {
    if (all_segments[ti].empty() && all_points[ti].empty())
      new_triangles.push_back(triangles[ti]);
    else
    {
#ifdef USE_FIXED_PROJECTION_TRAITS
      const std::array<typename EK::Point_3, 3>& t = triangles[ti];
      auto is_constant_in_dim = [](const std::array<typename EK::Point_3, 3>& t, int dim)
      {
        return t[0][dim] == t[1][dim] && t[0][dim] != t[2][dim];
      };

      typename EK::Vector_3 orth = CGAL::normal(t[0], t[1], t[2]); // TODO::avoid construction?
      int c = CGAL::abs(orth[0]) > CGAL::abs(orth[1]) ? 0 : 1;
      c = CGAL::abs(orth[2]) > CGAL::abs(orth[c]) ? 2 : c;

      if (c == 0) {
        autorefine_impl::generate_subtriangles<EK, 0>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      }
      else if (c == 1) {
        autorefine_impl::generate_subtriangles<EK, 1>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      }
      else if (c == 2) {
        autorefine_impl::generate_subtriangles<EK, 2>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      }
#else
      autorefine_impl::generate_subtriangles<EK>(ti, all_segments_ids[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, coplanar_triangles, triangles, new_triangles);
#endif
    }

#ifdef USE_PROGRESS_DISPLAY
    ++pd;
#endif
  };

#ifdef USE_DEBUG_PARALLEL_TIMERS
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

#ifdef USE_DEBUG_PARALLEL_TIMERS
  t.stop();
  std::cout << t.time() << " sec. for #1" << std::endl;
  t.reset();
#endif

  // brute force output: create a soup, orient and to-mesh
  CGAL_PMP_AUTOREFINE_VERBOSE("create output soup");

  Cartesian_converter<EK, GT> to_input;
  // TODO: reuse the fact that maps per triangle are already sorted

#ifdef CGAL_LINKED_WITH_TBB
  std::conditional_t<parallel_execution,
                     tbb::concurrent_map<EK::Point_3, std::size_t>,
                     std::map<EK::Point_3, std::size_t>> point_id_map;
#else
  std::map<EK::Point_3, std::size_t> point_id_map;
#endif

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

  TriIdsRange soup_triangles_out;
  soup_triangles_out.reserve(soup_triangles.size()); // TODO: remove #deg tri?

  // raw copy of input triangles with no intersection
  for (Input_TID f=0; f<soup_triangles.size(); ++f)
  {
    if (is_degen[f]) continue; //skip degenerate faces

    int tiid = tri_inter_ids[f];
    if (tiid == -1)
    {
      soup_triangles_out.push_back(
        {soup_triangles[f][0], soup_triangles[f][1], soup_triangles[f][2]}
      );
    }
  }

  // import refined triangles

  // Lambda to retrieve the id of intersection points
  auto get_point_id = [&](const typename EK::Point_3& pt)
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

#ifdef USE_DEBUG_PARALLEL_TIMERS
  t.start();
#endif

  std::size_t offset = soup_triangles_out.size();
#ifdef USE_DEBUG_PARALLEL_TIMERS
  std::string mode = "parallel";
#endif

  //TODO: 100 should be fined tune and depends on #threads
#ifdef CGAL_LINKED_WITH_TBB
  if(parallel_execution && new_triangles.size() > 100)
  {
#ifdef SET_POINT_IDS_USING_MUTEX
    //option 1 (using a mutex)
    CGAL_MUTEX point_container_mutex;
    /// Lambda concurrent_get_point_id()
    auto concurrent_get_point_id = [&](const typename EK::Point_3 pt)
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
    std::vector<std::array<tbb::concurrent_map<EK::Point_3, std::size_t>::iterator, 3>> triangle_buffer(new_triangles.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, new_triangles.size()),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t ti = r.begin(); ti != r.end(); ++ti) {
          const std::array<EK::Point_3, 3>& t = new_triangles[ti];
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
    tbb::concurrent_vector<tbb::concurrent_map<EK::Point_3, std::size_t>::iterator> iterators;
    auto concurrent_get_point_id = [&](const typename EK::Point_3 pt)
    {
      auto insert_res = point_id_map.insert(std::make_pair(pt, -1));
      if (insert_res.second)
        iterators.push_back(insert_res.first);
      return insert_res.first;
    };

    //use map iterator triple for triangles to create them concurrently and safely
    soup_triangles_out.resize(offset + new_triangles.size());
    std::vector<std::array<tbb::concurrent_map<EK::Point_3, std::size_t>::iterator, 3>> triangle_buffer(new_triangles.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, new_triangles.size()),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t ti = r.begin(); ti != r.end(); ++ti) {
          if (offset + ti > soup_triangles_out.size()) {
              std::cout << "ti = " << ti << std::endl;
          }
          const std::array<EK::Point_3, 3>& t = new_triangles[ti];
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
#ifdef USE_DEBUG_PARALLEL_TIMERS
    mode = "sequential";
#endif
    soup_triangles_out.reserve(offset + new_triangles.size());
    for (const std::array<EK::Point_3,3>& t : new_triangles)
      soup_triangles_out.push_back({ get_point_id(t[0]), get_point_id(t[1]), get_point_id(t[2])});
  }



#ifdef USE_DEBUG_PARALLEL_TIMERS
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
 * @tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm input triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
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

  autorefine_triangle_soup(soup_points, soup_triangles);

  clear(tm);
  repair_polygon_soup(soup_points, soup_triangles);

  duplicate_non_manifold_edges_in_polygon_soup(soup_points, soup_triangles);
  polygon_soup_to_polygon_mesh(soup_points, soup_triangles, tm);
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
