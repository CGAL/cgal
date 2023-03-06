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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// output
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#ifndef CGAL_PMP_AUTOREFINE_VERBOSE
#define CGAL_PMP_AUTOREFINE_VERBOSE(MSG)
#endif

#include <vector>

//#define TEST_RESOLVE_INTERSECTION
//#define DEDUPLICATE_SEGMENTS
//#define USE_FIXED_PROJECTION_TRAITS
//#define DEBUG_DEPTH

#ifdef USE_FIXED_PROJECTION_TRAITS
#include <CGAL/Kernel_23/internal/Projection_traits_3.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

#ifndef DOXYGEN_RUNNING
namespace autorefine_impl {

enum Segment_inter_type { NO_INTERSECTION=0, COPLANAR_SEGMENTS, POINT_INTERSECTION };

template <class K>
Segment_inter_type
do_coplanar_segments_intersect(const std::array<typename K::Point_3, 2>& s1,
                               const std::array<typename K::Point_3, 2>& s2,
                               const K& k = K())
{
  // supporting_line intersects: points are coplanar
  typename K::Coplanar_orientation_3 cpl_orient=k.coplanar_orientation_3_object();
  ::CGAL::Orientation or1 = cpl_orient(s1[0], s1[1], s2[0]);
  ::CGAL::Orientation or2 = cpl_orient(s1[0], s1[1], s2[1]);

  if(or1 == COLLINEAR && or2 == COLLINEAR)
  {
    // segments are collinear
    typename K::Collinear_are_ordered_along_line_3 cln_order = k.collinear_are_ordered_along_line_3_object();
    return (cln_order(s1[0], s2[0], s1[1]) ||
            cln_order(s1[0], s2[1], s1[1]) ||
            cln_order(s2[0], s1[0], s2[1])) ? COPLANAR_SEGMENTS : NO_INTERSECTION;
  }

  if(or1 != or2)
  {
    or1 = cpl_orient(s2[0], s2[1], s1[0]);
    return (or1 == COLLINEAR || or1 != cpl_orient(s2[0], s2[1], s1[1])) ? POINT_INTERSECTION : NO_INTERSECTION;
  }

  return NO_INTERSECTION;
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

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
  auto print_points = [&]()
  {
    for(auto p : l_inter_pts) std::cout << "  (" << p.id1() << "," << p.id2() << ",[" << p.alpha << "]) "; std::cout <<"\n";
  };
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
#endif

  //intersect t2 with the three half planes which intersection defines t1
  K k;
  intersection_coplanar_triangles_cutoff(p1,q1,r1,0,p2,q2,r2,k,l_inter_pts); //line p1q1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
#endif
  intersection_coplanar_triangles_cutoff(q1,r1,p1,1,p2,q2,r2,k,l_inter_pts); //line q1r1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
#endif
  intersection_coplanar_triangles_cutoff(r1,p1,q1,2,p2,q2,r2,k,l_inter_pts); //line r1p1
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << l_inter_pts.size() << "\n";
  print_points();
#endif

  for (const Intersections::internal::Point_on_triangle<K>& pot : l_inter_pts)
    inter_pts.push_back( pot.point(p1,q1,r1,p2,q2,r2,k) );
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
void collect_intersections(const std::array<typename K::Point_3, 3>& t1,
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

    return;
  }

  for (int i=0; i<3; ++i)
  {
    int j=(i+1)%3;
    test_edge<K>(t1[i], t1[j], t2[0], t2[1], t2[2], ori[i], ori[j], inter_pts);
    //~ if (inter_pts.size()>1) return;
  }

  // test edges of t2 vs t1
  for (int i=0; i<3; ++i)
    ori[i] = orientation(t1[0],t1[1],t1[2],t2[i]);
  for (int i=0; i<3; ++i)
  {
    int j=(i+1)%3;
    test_edge<K>(t2[i], t2[j], t1[0], t1[1], t1[2], ori[i], ori[j], inter_pts);
    //~ if (inter_pts.size()>1) return;
  }

  // because we don't handle intersection type and can have edge-edge edge-vertex duplicates
  std::sort(inter_pts.begin(), inter_pts.end());
  auto last = std::unique(inter_pts.begin(), inter_pts.end());
  inter_pts.erase(last, inter_pts.end());


  for (auto p : inter_pts)
    if (depth(p)>2) throw std::runtime_error("Depth is not 2: "+std::to_string(depth(p)));
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <class EK
#ifdef USE_FIXED_PROJECTION_TRAITS
, int dim
#endif
>
void generate_subtriangles(std::size_t ti,
                           std::vector<std::array<typename EK::Point_3, 2>>& segments,
                           const std::vector<typename EK::Point_3>& points,
                           const std::vector<std::size_t>& in_triangle_ids,
                           const std::set<std::pair<std::size_t, std::size_t> >& intersecting_triangles,
                           const std::vector<std::array<typename EK::Point_3,3>>& triangles,
                           std::vector<std::array<typename EK::Point_3,3>>& new_triangles)
{
  //~ std::cout << "generate_subtriangles()\n";
  std::cout << std::setprecision(17);

#ifdef USE_FIXED_PROJECTION_TRAITS
  typedef ::CGAL::internal::Projection_traits_3<EK, dim> P_traits;
#else
  typedef CGAL::Projection_traits_3<EK> P_traits;
#endif
  typedef CGAL::Exact_intersections_tag Itag;

  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits
#ifndef TEST_RESOLVE_INTERSECTION
  ,Default, Itag
#endif
> CDT_2;
  //typedef CGAL::Constrained_triangulation_plus_2<CDT_base> CDT;
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
  P_traits cdt_traits(n);
  bool orientation_flipped = false;
  if ( typename EK::Less_xyz_3()(o+n,o) )
  {
    n=-n;
    orientation_flipped = true;
  }
  CDT cdt(cdt_traits);
  cdt.insert_outside_affine_hull(t[0]);
  cdt.insert_outside_affine_hull(t[1]);
  typename CDT::Vertex_handle v = cdt.tds().insert_dim_up(cdt.infinite_vertex(), orientation_flipped);
  v->set_point(t[2]);
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
    //~ std::cout << "nbs " << nbs << "\n";

    //~ if (nbs==8)
    //~ {
      //~ for (std::size_t i = 0; i<nbs; ++i)
        //~ std::ofstream("cst_"+std::to_string(i)+".polylines.txt") << std::setprecision(17) << "2 " << segments[i][0] << " " << segments[i][1] << "\n";
    //~ }

    auto supporting_plane = [](const std::array<typename EK::Point_3, 3>& t)
    {
      return typename EK::Plane_3(t[0], t[1], t[2]);
    };

    std::vector< std::vector<typename EK::Point_3> > points_on_segments(nbs);
    for (std::size_t i = 0; i<nbs-1; ++i)
    {
      for (std::size_t j = i+1; j<nbs; ++j)
      {
        if (intersecting_triangles.count(CGAL::make_sorted_pair(in_triangle_ids[i], in_triangle_ids[j]))!=0)
        {
          Segment_inter_type seg_inter_type = do_coplanar_segments_intersect<EK>(segments[i], segments[j]);
          switch(seg_inter_type)
          {
            case POINT_INTERSECTION:
            {
              auto res = CGAL::intersection(supporting_plane(triangles[in_triangle_ids[i]]),
                                            supporting_plane(triangles[in_triangle_ids[j]]),
                                            supporting_plane(triangles[ti]));

              if (const typename EK::Point_3* pt_ptr = boost::get<typename EK::Point_3>(&(*res)))
              {
                points_on_segments[i].push_back(*pt_ptr);
                points_on_segments[j].push_back(*pt_ptr);

                //~ std::cout << "new inter " << *pt_ptr << " (" << depth(points_on_segments[i].back()) << ")" << "\n";

              }
            }
            // break; No break because of the coplanar case
            case COPLANAR_SEGMENTS:
            {
              // We can have hard cases if two triangles are coplanar....

              //~ std::cout << "coplanar inter: " << i << " " << j << "\n";

              typename EK::Segment_3 s1(segments[i][0], segments[i][1]);
              typename EK::Segment_3 s2(segments[j][0], segments[j][1]);// TODO: avoid this construction
              auto inter = CGAL::intersection(s1, s2);

              if (inter == boost::none) throw std::runtime_error("Unexpected case #2");

              if (const typename EK::Point_3* pt_ptr = boost::get<typename EK::Point_3>(&(*inter)))
              {
                points_on_segments[i].push_back(*pt_ptr);
                points_on_segments[j].push_back(*pt_ptr);

                //~ std::cout << "new inter bis " << *pt_ptr << " (" << depth(points_on_segments[i].back()) << ")" <<  "\n";
              }
              else
              {
                if (const typename EK::Segment_3* seg_ptr = boost::get<typename EK::Segment_3>(&(*inter)))
                {
                  points_on_segments[i].push_back(seg_ptr->source());
                  points_on_segments[j].push_back(seg_ptr->source());
                  points_on_segments[i].push_back(seg_ptr->target());
                  points_on_segments[j].push_back(seg_ptr->target());

                  //~ std::cout << "new inter seg " << *seg_ptr << " (" << depth(*seg_ptr) << ")" <<  "\n";

                }
                else
                  throw std::runtime_error("BOOM\n");
              }

#if 0
              //this code works if triangles are not coplanar
              // coplanar intersection that is not a point
              int coord = 0;
              const typename EK::Segment_3& s = segments[i];
              typename EK::Point_3 src = s[0], tgt=s[1];
              if (src.x()==tgt.x())
              {
                coord=1;
                if (src.y()==tgt.y())
                  coord==2;
              }

              std::vector<typename EK::Point_3> tmp_pts = {
                src, tgt, segments[j][0], segments[j][1] };

              std::sort(tmp_pts.begin(), tmp_pts.end(),
                        [coord](const typename EK::Point_3& p, const typename EK::Point_3& q)
                        {return p[coord]<q[coord];});

              points_on_segments[i].push_back(tmp_pts[1]);
              points_on_segments[i].push_back(tmp_pts[2]);
              points_on_segments[j].push_back(tmp_pts[1]);
              points_on_segments[j].push_back(tmp_pts[2]);
#endif
              //~ std::cout << "new inter coli " << segments[j][0] << "\n";
              //~ std::cout << "new inter coli " << segments[j][1] << "\n";
              //~ std::cout << "new inter coli " << segments[i][0] << "\n";
              //~ std::cout << "new inter coli " << segments[i][1] << "\n";

              //~ points_on_segments[j].push_back(*pt_ptr);



              //~ std::cerr << "ERROR: intersection is a segment\n";
              //~ std::cout << std::setprecision(17);
              //~ exact(segments[i]);
              //~ exact(segments[j]);
              //~ std::cout << segments[i] << "\n";
              //~ std::cout << segments[j] << "\n";
              //~ debug << "4 " << triangles[in_triangle_ids[i]] << " " << triangles[in_triangle_ids[i]][0] << "\n";
              //~ debug << "4 " << triangles[in_triangle_ids[j]] << " " << triangles[in_triangle_ids[j]][0] << "\n";
              //~ debug << "4 " << triangles[ti] << " " << triangles[ti][0] << "\n";
              //~ exit(1);
            }
            break;
            default:
            break;
          }
        }
      }
    }

    std::vector<typename EK::Point_3> cst_points;
    std::vector<std::pair<std::size_t, std::size_t>> csts;
    for (std::size_t i = 0; i<nbs; ++i)
    {
      if(!points_on_segments[i].empty())
      {
        // TODO: predicate on input triangles
        int coord = 0;
        const std::array<typename EK::Point_3, 2>& s = segments[i];
        typename EK::Point_3 src = s[0], tgt=s[1];
        if (src.x()==tgt.x())
        {
          coord=1;
          if (src.y()==tgt.y())
            coord==2;
        }
        if (src[coord]>tgt[coord])
          std::swap(src, tgt);

        std::sort(points_on_segments[i].begin(), points_on_segments[i].end(), [coord](const typename EK::Point_3& p, const typename EK::Point_3& q){return p[coord]<q[coord];});
//TODO: use reserve on cst_points?
        std::size_t src_id=cst_points.size();
        cst_points.push_back(src);
        cst_points.insert(cst_points.end(), points_on_segments[i].begin(), points_on_segments[i].end());
        cst_points.push_back(tgt);


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

        for (std::size_t k=0; k<=points_on_segments[i].size(); ++k)
          csts.emplace_back(src_id+k, src_id+k+1);
      }
    }

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

    cdt.insert_constraints(cst_points.begin(), cst_points.end(), csts.begin(), csts.end());

    std::vector<std::array<typename EK::Point_3,2>> no_inter_segments;
    no_inter_segments.reserve(nbs);
    for (std::size_t i = 0; i<nbs; ++i)
      if(points_on_segments[i].empty())
        no_inter_segments.push_back(segments[i]);
    no_inter_segments.swap(segments);
  }
  //~ std::cout << "done\n";
#endif

  cdt.insert_constraints(segments.begin(), segments.end());
  cdt.insert(points.begin(), points.end());

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

template <class PointRange, class TriIdsRange, class Point_3, class NamedParameters = parameters::Default_named_parameters>
void autorefine_soup_output(const PointRange& input_points,
                            const TriIdsRange& id_triples,
                            std::vector<Point_3>& soup_points,
                            std::vector<std::array<std::size_t, 3> >& soup_triangles,
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

  typedef std::size_t Input_TID;
  typedef std::pair<Input_TID, Input_TID> Pair_of_triangle_ids;

  std::vector<Pair_of_triangle_ids> si_pairs;

  // collect intersecting pairs of triangles
  CGAL_PMP_AUTOREFINE_VERBOSE("collect intersecting pairs");
  triangle_soup_self_intersections<Concurrency_tag>(input_points, id_triples, std::back_inserter(si_pairs), np);

  if (si_pairs.empty()) return;

  // mark degenerate faces so that we can ignore them
  std::vector<bool> is_degen(id_triples.size(), false);

  for (const Pair_of_triangle_ids& p : si_pairs)
    if (p.first==p.second) // bbox inter reports (f,f) for degenerate faces
      is_degen[p.first] = true;

  // assign an id per triangle involved in an intersection
  // + the faces involved in the intersection
  std::vector<int> tri_inter_ids(id_triples.size(), -1);
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
      to_exact( get(pm, input_points[id_triples[f][0]]) ),
      to_exact( get(pm, input_points[id_triples[f][1]]) ),
      to_exact( get(pm, input_points[id_triples[f][2]]) ) );
  }

  std::vector< std::vector<std::array<EK::Point_3,2> > > all_segments(triangles.size());
  std::vector< std::vector<EK::Point_3> > all_points(triangles.size());
  std::vector< std::vector<std::size_t> > all_in_triangle_ids(triangles.size());

  CGAL_PMP_AUTOREFINE_VERBOSE("compute intersections");

  std::set<std::pair<std::size_t, std::size_t> > intersecting_triangles;
  for (const Pair_of_triangle_ids& p : si_pairs)
  {
    int i1 = tri_inter_ids[p.first],
        i2 = tri_inter_ids[p.second];

    if (i1==-1 || i2==-1) continue; //skip degenerate faces

    const std::array<typename EK::Point_3, 3>& t1 = triangles[i1];
    const std::array<typename EK::Point_3, 3>& t2 = triangles[i2];

    std::vector<typename EK::Point_3> inter_pts;
    autorefine_impl::collect_intersections<EK>(t1, t2, inter_pts);

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
    }
  }

  // deduplicate inserted segments
  Cartesian_converter<EK, GT> to_input;
  std::map<EK::Point_3, std::size_t> point_id_map;
#if ! defined(CGAL_NDEBUG) || defined(CGAL_DEBUG_PMP_AUTOREFINE)
  std::vector<EK::Point_3> exact_soup_points;
#endif

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

  // filter duplicated segments
#ifdef DEDUPLICATE_SEGMENTS
  for(std::size_t ti=0; ti<triangles.size(); ++ti)
  {
    if (!all_segments[ti].empty())
    {
      std::size_t nbs = all_segments[ti].size();
      std::vector<std::array<EK::Point_3,2>> filtered_segments;
      std::vector<std::size_t> filtered_in_triangle_ids;
      filtered_segments.reserve(nbs);
      std::set<std::pair<std::size_t, std::size_t>> segset;
      for (std::size_t si=0; si<nbs; ++si)
      {
        EK::Point_3 src = all_segments[ti][si][0],
                    tgt = all_segments[ti][si][1];
        if (segset.insert(
              CGAL::make_sorted_pair( get_point_id(src),
                                      get_point_id(tgt))).second)
        {
          filtered_segments.push_back(all_segments[ti][si]);
          filtered_in_triangle_ids.push_back(all_in_triangle_ids[ti][si]);
        }
      }
      if (filtered_segments.size()!=nbs)
      {
        filtered_segments.swap(all_segments[ti]);
        filtered_in_triangle_ids.swap(all_in_triangle_ids[ti]);
      }
    }
  }
#endif

  CGAL_PMP_AUTOREFINE_VERBOSE("triangulate faces");
  // now refine triangles
  std::vector<std::array<EK::Point_3,3>> new_triangles;
  for(std::size_t ti=0; ti<triangles.size(); ++ti)
  {
    if (all_segments[ti].empty() && all_points[ti].empty())
      new_triangles.push_back(triangles[ti]);
    else
    {
      #ifdef USE_FIXED_PROJECTION_TRAITS
      const std::array<typename EK::Point_3, 3>& t = triangles[ti];
      auto is_constant_in_dim = [](const std::array<typename EK::Point_3, 3>& t, int dim)
      {
        return t[0][dim]==t[1][dim] && t[0][dim]!=t[2][dim];
      };

      typename EK::Vector_3 orth = CGAL::normal(t[0], t[1], t[2]); // TODO::avoid construction?
      int c = CGAL::abs(orth[0]) > CGAL::abs(orth[1]) ? 0 : 1;
      c = CGAL::abs(orth[2]) > CGAL::abs(orth[c]) ? 2 : c;

      if(c == 0) {
        autorefine_impl::generate_subtriangles<EK, 0>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      } else if(c == 1) {
        autorefine_impl::generate_subtriangles<EK, 1>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      } else if(c == 2) {
        autorefine_impl::generate_subtriangles<EK, 2>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      }
      #else
      autorefine_impl::generate_subtriangles<EK>(ti, all_segments[ti], all_points[ti], all_in_triangle_ids[ti], intersecting_triangles, triangles, new_triangles);
      #endif
    }


  }

  // brute force output: create a soup, orient and to-mesh
  CGAL_PMP_AUTOREFINE_VERBOSE("create output soup");

  std::vector <std::size_t> input_point_ids;
  input_point_ids.reserve(input_points.size());
  for (const auto& p : input_points)
    input_point_ids.push_back(get_point_id(to_exact(get(pm,p))));

  for (Input_TID f=0; f<id_triples.size(); ++f)
  {
    if (is_degen[f]) continue; //skip degenerate faces

    int tiid = tri_inter_ids[f];
    if (tiid == -1)
    {
      soup_triangles.emplace_back(
        CGAL::make_array(input_point_ids[id_triples[f][0]],
                         input_point_ids[id_triples[f][1]],
                         input_point_ids[id_triples[f][2]])
      );
    }
  }
  for (const std::array<EK::Point_3,3>& t : new_triangles)
  {
    soup_triangles.emplace_back(CGAL::make_array(get_point_id(t[0]), get_point_id(t[1]), get_point_id(t[2])));
  }

#ifndef CGAL_NDEBUG
  CGAL_PMP_AUTOREFINE_VERBOSE("check soup");
  CGAL_assertion( !does_triangle_soup_self_intersect(exact_soup_points, soup_triangles) );
#else
#ifdef CGAL_DEBUG_PMP_AUTOREFINE
  CGAL_PMP_AUTOREFINE_VERBOSE("check soup");
  if (does_triangle_soup_self_intersect(exact_soup_points, soup_triangles))
    throw std::runtime_error("invalid output");
#endif
#endif
  CGAL_PMP_AUTOREFINE_VERBOSE("done");
}
#endif

/**
 * \ingroup PMP_corefinement_grp
 * refines a triangle mesh so that no triangles intersects in their interior.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm input triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalParamNBegin{geom_traits}
 *   \cgalParamDescription{an instance of a geometric traits class}
 *   \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *   \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *   \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 * \cgalParamNEnd
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *  \cgalParamNEnd
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
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  std::vector<typename GT::Point_3> in_soup_points;
  std::vector<std::array<std::size_t, 3> > in_soup_triangles;
  std::vector<typename GT::Point_3> out_soup_points;
  std::vector<std::array<std::size_t, 3> > out_soup_triangles;

  polygon_mesh_to_polygon_soup(tm, in_soup_points, in_soup_triangles);

  autorefine_soup_output(in_soup_points, in_soup_triangles,
                         out_soup_points, out_soup_triangles);

  clear(tm);
  orient_polygon_soup(out_soup_points, out_soup_triangles);
  polygon_soup_to_polygon_mesh(out_soup_points, out_soup_triangles, tm);
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
