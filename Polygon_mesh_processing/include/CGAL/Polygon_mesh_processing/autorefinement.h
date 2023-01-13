//TODO: add for soup face the id of the input face. not sure it is easy to report intersection edge as a pair of vertex id

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

namespace CGAL {
namespace Polygon_mesh_processing {

#ifndef DOXYGEN_RUNNING
namespace autorefine_impl {

template <class EK>
void generate_subtriangles(const typename EK::Triangle_3& t,
                           const std::vector<typename EK::Segment_3>& segments,
                           const std::vector<typename EK::Point_3>& points,
                           std::vector<typename EK::Triangle_3>& new_triangles)
{
  typedef CGAL::Projection_traits_3<EK> P_traits;
  typedef CGAL::Exact_intersections_tag Itag;

  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, Default, Itag> CDT_2;
  //typedef CGAL::Constrained_triangulation_plus_2<CDT_base> CDT;
  typedef CDT_2 CDT;

  // positive triangle normal
  typename EK::Vector_3 n = normal(t[0], t[1], t[2]);
  typename EK::Point_3 o(0,0,0);

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
        new_triangles.emplace_back(fh->vertex(0)->point(),
                                   fh->vertex(cdt.cw(0))->point(),
                                   fh->vertex(cdt.ccw(0))->point());
      else
        new_triangles.emplace_back(fh->vertex(0)->point(),
                                   fh->vertex(cdt.ccw(0))->point(),
                                   fh->vertex(cdt.cw(0))->point());
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

template <class EK>
struct Intersection_visitor
{
  std::vector< std::vector<typename EK::Segment_3> >& all_segments;
  std::vector< std::vector<typename EK::Point_3> >& all_points;
  std::pair<int, int> ids;

  Intersection_visitor(std::vector< std::vector<typename EK::Segment_3> >& all_segments,
                       std::vector< std::vector<typename EK::Point_3> >& all_points)
    : all_segments (all_segments)
    , all_points(all_points)
  {}

  void set_triangle_ids(int i1, int i2)
  {
    ids = {i1, i2};
  }

  typedef void result_type;
  void operator()(const typename EK::Point_3& p)
  {
    all_points[ids.first].push_back(p);
    all_points[ids.second].push_back(p);
  }

  void operator()(const typename EK::Segment_3& s)
  {
    all_segments[ids.first].push_back(s);
    all_segments[ids.second].push_back(s);
  }

  void operator()(const typename EK::Triangle_3& t)
  {
    for (std::size_t i=0; i<3; ++i)
    {
      typename EK::Segment_3 s(t[i], t[(i+1)%3]);
      all_segments[ids.first].push_back(s);
      all_segments[ids.second].push_back(s);
    }

  }

  void operator()(const std::vector<typename EK::Point_3>& poly)
  {
    std::size_t nbp = poly.size();
    for (std::size_t i=0; i<nbp; ++i)
    {
      typename EK::Segment_3 s(poly[i], poly[(i+1)%nbp]);
      all_segments[ids.first].push_back(s);
      all_segments[ids.second].push_back(s);
    }
  }
};

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
  std::vector< EK::Triangle_3 > triangles(tiid+1);
  Cartesian_converter<GT, EK> to_exact;

  for(Input_TID f : intersected_faces)
  {
    triangles[tri_inter_ids[f]]= EK::Triangle_3(
      to_exact( get(pm, input_points[id_triples[f][0]]) ),
      to_exact( get(pm, input_points[id_triples[f][1]]) ),
      to_exact( get(pm, input_points[id_triples[f][2]]) ) );
  }

  std::vector< std::vector<EK::Segment_3> > all_segments(triangles.size());
  std::vector< std::vector<EK::Point_3> > all_points(triangles.size());

  CGAL_PMP_AUTOREFINE_VERBOSE("compute intersections");
  typename EK::Intersect_3 intersection = EK().intersect_3_object();
  autorefine_impl::Intersection_visitor<EK> intersection_visitor(all_segments, all_points);

  for (const Pair_of_triangle_ids& p : si_pairs)
  {
    int i1 = tri_inter_ids[p.first],
        i2 = tri_inter_ids[p.second];

    if (i1==-1 || i2==-1) continue; //skip degenerate faces

    const EK::Triangle_3& t1 = triangles[i1];
    const EK::Triangle_3& t2 = triangles[i2];

    auto inter = intersection(t1, t2);

    if (inter != boost::none)
    {
      intersection_visitor.set_triangle_ids(i1, i2);
      boost::apply_visitor(intersection_visitor, *inter);
    }
  }

  CGAL_PMP_AUTOREFINE_VERBOSE("triangulate faces");
  // now refine triangles
  std::vector<EK::Triangle_3> new_triangles;
  for(std::size_t ti=0; ti<triangles.size(); ++ti)
  {
    if (all_segments[ti].empty() && all_points[ti].empty())
      new_triangles.push_back(triangles[ti]);
    else
      autorefine_impl::generate_subtriangles<EK>(triangles[ti], all_segments[ti], all_points[ti], new_triangles);
  }


  // brute force output: create a soup, orient and to-mesh
  CGAL_PMP_AUTOREFINE_VERBOSE("create output soup");
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
  for (const typename EK::Triangle_3& t : new_triangles)
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
