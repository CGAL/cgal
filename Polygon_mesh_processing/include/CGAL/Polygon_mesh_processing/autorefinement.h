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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// output
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#ifndef NDEBUG
// debug
#include <CGAL/Surface_mesh.h>
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
  // typedef CGAL::Exact_intersections_tag Itag;
  typedef CGAL::No_constraint_intersection_requiring_constructions_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, Default, Itag> CDT;


  P_traits cdt_traits(normal(t[0], t[1], t[2]));
  CDT cdt(cdt_traits);

  cdt.insert_outside_affine_hull(t[0]);
  cdt.insert_outside_affine_hull(t[1]);
  typename CDT::Vertex_handle v = cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
  v->set_point(t[2]);

  for (const typename EK::Segment_3& s : segments)
    cdt.insert_constraint(s[0], s[1]);

  for (const typename EK::Point_3& p : points)
    cdt.insert(p);

  for (typename CDT::Face_handle fh : cdt.finite_face_handles())
  {
    new_triangles.emplace_back(fh->vertex(0)->point(),
                               fh->vertex(cdt.ccw(0))->point(),
                               fh->vertex(cdt.cw(0))->point());
  }
}

template <class EK>
struct Intersection_visitor
{
  std::vector<typename EK::Segment_3>& all_segments_1;
  std::vector<typename EK::Segment_3>& all_segments_2;
  std::vector<typename EK::Point_3>& all_points_1;
  std::vector<typename EK::Point_3>& all_points_2;

  Intersection_visitor(std::vector<typename EK::Segment_3>& all_segments_1,
                       std::vector<typename EK::Segment_3>& all_segments_2,
                       std::vector<typename EK::Point_3>& all_points_1,
                       std::vector<typename EK::Point_3>& all_points_2)
    : all_segments_1(all_segments_1)
    , all_segments_2(all_segments_2)
    , all_points_1(all_points_1)
    , all_points_2(all_points_2)
  {}

  typedef void result_type;
  void operator()(const typename EK::Point_3& p)
  {
    all_points_1.push_back(p);
    all_points_2.push_back(p);
  }

  void operator()(const typename EK::Segment_3& s)
  {
    all_segments_1.push_back(s);
    all_segments_2.push_back(s);
  }

  void operator()(const typename EK::Triangle_3& t)
  {
    for (std::size_t i=1; i<3; ++i)
    {
      typename EK::Segment_3 s(t[i-1], t[i]);
      all_segments_1.push_back(s);
      all_segments_2.push_back(s);
    }
  }

  void operator()(const std::vector<typename EK::Point_3>& poly)
  {
    std::size_t nbp = poly.size();
    for (std::size_t i=1; i<nbp; ++i)
    {
      typename EK::Segment_3 s(poly[i-1], poly[i]);
      all_segments_1.push_back(s);
      all_segments_2.push_back(s);
    }
  }
};

template <class EK, class TriangleMesh, class VPM, class TID_Map>
bool is_output_valid(TriangleMesh& tm , VPM vpm, TID_Map tid_map, const std::vector< std::vector<Triangle_3<EK>>>& triangles)
{
  typedef typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type IK;
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  std::vector<typename EK::Point_3> soup_points;
  std::vector<std::array<std::size_t, 3> > soup_triangles;
  Cartesian_converter<IK, EK> to_exact;
  std::map<typename EK::Point_3, std::size_t> point_id_map;

  auto get_point_id = [&](const typename EK::Point_3& pt)
  {
    auto insert_res = point_id_map.insert(std::make_pair(pt, soup_points.size()));
    if (insert_res.second)
      soup_points.push_back(pt);
    return insert_res.first->second;
  };

  for (face_descriptor f : faces(tm))
  {
    int tid = get(tid_map, f);
    if (tid == -1)
    {
      halfedge_descriptor h = halfedge(f, tm);
      soup_triangles.emplace_back(
        CGAL::make_array(get_point_id(to_exact(get(vpm,source(h, tm)))),
                         get_point_id(to_exact(get(vpm,target(h, tm)))),
                         get_point_id(to_exact(get(vpm,target(next(h, tm), tm)))))
      );
    }
    else
    {
      for (const typename EK::Triangle_3& t : triangles[tid])
      {
        soup_triangles.emplace_back(CGAL::make_array(get_point_id(t[0]), get_point_id(t[1]), get_point_id(t[2])));
      }
    }
  }

  typedef Surface_mesh<typename EK::Point_3> Exact_mesh;
  Exact_mesh etm;
  orient_polygon_soup(soup_points, soup_triangles);
  polygon_soup_to_polygon_mesh(soup_points, soup_triangles, etm);
  typename Exact_mesh::Property_map<typename Exact_mesh::Vertex_index, bool> is_border_map =
    etm.template add_property_map<typename Exact_mesh::Vertex_index, bool>("v:is_border", false).first;
  for(typename Exact_mesh::Halfedge_index h : etm.halfedges())
  {
    if (CGAL::is_border(h, etm))
      is_border_map[target(h, etm)] = true;
  }

  //TODO: double check me
  auto skip_faces = [&](const std::pair<typename Exact_mesh::Face_index, typename Exact_mesh::Face_index>& p)
  {
    typename Exact_mesh::Halfedge_index h1 = etm.halfedge(p.first), h2=etm.halfedge(p.second);

    boost::container::small_vector<typename Exact_mesh::Halfedge_index, 3> bv1;
    if (is_border_map[source(h1, etm)]) bv1.push_back(prev(h1, etm));
    if (is_border_map[target(h1, etm)]) bv1.push_back(h1);
    if (is_border_map[target(next(h1, etm), etm)]) bv1.push_back(next(h1, etm));
    if (bv1.empty()) return false;

    boost::container::small_vector<typename Exact_mesh::Halfedge_index, 3> bv2;
    if (is_border_map[source(h2, etm)]) bv2.push_back(prev(h2, etm));
    if (is_border_map[target(h2, etm)]) bv2.push_back(h2);
    if (is_border_map[target(next(h2, etm), etm)]) bv2.push_back(next(h2, etm));
    if (bv2.empty()) return false;

    //collect identical border vertices
    boost::container::small_vector<std::pair<typename Exact_mesh::Halfedge_index, typename Exact_mesh::Halfedge_index>, 3> common;
    for(typename Exact_mesh::Halfedge_index h1 : bv1)
      for(typename Exact_mesh::Halfedge_index h2 : bv2)
        if (etm.point(target(h1, etm))==etm.point(target(h2,etm)))
          common.push_back(std::make_pair(h1,h2));

    if (common.empty()) return false;

    switch (common.size())
    {
      case 1:
      {
        // geometric check if the opposite segments intersect the triangles
        const typename EK::Triangle_3 t1(etm.point(source(h1,etm)),
                                         etm.point(target(h1,etm)),
                                         etm.point(target(next(h1,etm),etm)));
        const typename EK::Triangle_3 t2(etm.point(source(h2,etm)),
                                         etm.point(target(h2,etm)),
                                         etm.point(target(next(h2,etm),etm)));

        const typename EK::Segment_3 s1(etm.point(source(common[0].first,etm)), etm.point(target(next(common[0].first,etm),etm)));
        const typename EK::Segment_3 s2(etm.point(source(common[0].second,etm)), etm.point(target(next(common[0].second,etm),etm)));

        if(do_intersect(t1, s2) || do_intersect(t2, s1))
          return false;
        return true;
      }
      case 2:
      {
        // shared edge
        h1 = next(common[0].first, etm) == common[1].first ? common[1].first : common[0].first;
        h2 = next(common[0].second, etm) == common[1].second ? common[1].second : common[0].second;

        if ( is_border(etm.opposite(h1), etm) &&
             is_border(etm.opposite(h2), etm) )
        {
          if( CGAL::coplanar(etm.point(source(h1,etm)),
                             etm.point(target(h1,etm)),
                             etm.point(target(etm.next(h1),etm)),
                             etm.point(source(h1,etm))) &&
              CGAL::coplanar_orientation(etm.point(source(h1,etm)),
                             etm.point(target(h1,etm)),
                             etm.point(target(etm.next(h1),etm)),
                             etm.point(source(h1,etm))) == CGAL::POSITIVE)
          {
            return false;
          }
          return true;
        }
        else
        {
          // TODO: 2 identical border vertices, no common edge on the boundary. Not sure what to do
          return false;
        }
      }
      default: // size == 3
        return true;
    }
  };

  std::vector< std::pair<typename Exact_mesh::Face_index, typename Exact_mesh::Face_index> > si_faces;
  self_intersections(etm,
                     CGAL::filter_output_iterator(std::back_inserter(si_faces),
                     skip_faces));

  return si_faces.empty();
}

} // end of autorefine_impl

template <class TriangleMesh, class Point_3, class NamedParameters = parameters::Default_named_parameters>
void autorefine_soup_output(const TriangleMesh& tm,
                            std::vector<Point_3>& soup_points,
                            std::vector<std::array<std::size_t, 3> >& soup_triangles,
                            const NamedParameters& np = parameters::default_values())
{
  //TODO: what about degenerate faces?

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, tm));

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef std::pair<face_descriptor, face_descriptor> Pair_of_faces;

  std::vector<Pair_of_faces> si_pairs;

  // collect intersecting pairs of triangles
  self_intersections<Concurrency_tag>(tm, std::back_inserter(si_pairs), np);

  if (si_pairs.empty()) return;

  // assign an id per triangle involved in an intersection
  // + the faces involved in the intersection
  typedef CGAL::dynamic_face_property_t<int> Face_property_tag;
  typedef typename boost::property_map<TriangleMesh, Face_property_tag>::const_type Triangle_id_map;

  Triangle_id_map tid_map = get(Face_property_tag(), tm);
  for (face_descriptor f : faces(tm))
    put(tid_map, f, -1);

  std::vector<face_descriptor> intersected_faces;
  int tid=-1;
  for (const Pair_of_faces& p : si_pairs)
  {
    if (get(tid_map, p.first)==-1)
    {
      put(tid_map, p.first, ++tid);
      intersected_faces.push_back(p.first);
    }
    if (get(tid_map, p.second)==-1)
    {
      put(tid_map, p.second, ++tid);
      intersected_faces.push_back(p.second);
    }
  }

  // init the vector of triangles used for the autorefinement of triangles
  typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
  std::vector< std::vector<EK::Triangle_3> > triangles(tid+1);
  Cartesian_converter<GT, EK> to_exact;

  for(face_descriptor f : intersected_faces)
  {
    halfedge_descriptor h = halfedge(f, tm);
    triangles[get(tid_map, f)].emplace_back(
      to_exact( get(vpm, source(h, tm)) ),
      to_exact( get(vpm, target(h, tm)) ),
      to_exact( get(vpm, target(next(h, tm), tm)) ) );
  }

  typename EK::Intersect_3 intersection = EK().intersect_3_object();
  for (const Pair_of_faces& p : si_pairs)
  {
    int i1 = get(tid_map, p.first),
        i2 = get(tid_map, p.second);


    std::size_t nbt_1 = triangles[i1].size(),
                nbt_2 = triangles[i2].size();

    std::vector< std::vector<EK::Segment_3> > all_segments_1(nbt_1);
    std::vector< std::vector<EK::Segment_3> > all_segments_2(nbt_2);
    std::vector< std::vector<EK::Point_3> > all_points_1(nbt_1);
    std::vector< std::vector<EK::Point_3> > all_points_2(nbt_2);

    std::vector <EK::Triangle_3> t1_subtriangles, t2_subtriangles;
    for (std::size_t it1=0; it1<nbt_1; ++it1)
    {
      for (std::size_t it2=0; it2<nbt_2; ++it2)
      {
        const EK::Triangle_3& t1 = triangles[i1][it1];
        const EK::Triangle_3& t2 = triangles[i2][it2];

        auto inter = intersection(t1, t2);

        if (inter != boost::none)
        {
          autorefine_impl::Intersection_visitor<EK> intersection_visitor(all_segments_1[it1],  all_segments_2[it2],
                                                                         all_points_1[it1], all_points_2[it2]);

          boost::apply_visitor(intersection_visitor, *inter);
        }
      }
    }

    // now refine triangles
    std::vector<EK::Triangle_3> new_triangles;
    for(std::size_t it1=0; it1<nbt_1; ++it1)
    {
      if (all_segments_1[it1].empty() && all_points_1[it1].empty())
        new_triangles.push_back(triangles[i1][it1]);
      else
        autorefine_impl::generate_subtriangles<EK>(triangles[i1][it1], all_segments_1[it1], all_points_1[it1], new_triangles);
    }
    triangles[i1].swap(new_triangles);
    new_triangles.clear();
    for(std::size_t it2=0; it2<nbt_2; ++it2)
    {
      if (all_segments_2[it2].empty() && all_points_2[it2].empty())
        new_triangles.push_back(triangles[i2][it2]);
      else
        autorefine_impl::generate_subtriangles<EK>(triangles[i2][it2], all_segments_2[it2], all_points_2[it2], new_triangles);
    }
    triangles[i2].swap(new_triangles);
  }

  CGAL_assertion( autorefine_impl::is_output_valid(tm, vpm, tid_map, triangles) );

  // brute force output: create a soup, orient and to-mesh
  // WARNING: there is no reason when using double that identical exact points are identical in double
  Cartesian_converter<EK, GT> to_input;
  std::map<EK::Point_3, std::size_t> point_id_map;

  for (vertex_descriptor v : vertices(tm))
  {
    if (point_id_map.insert(std::make_pair(to_exact(get(vpm,v)), soup_points.size())).second)
      soup_points.push_back(get(vpm,v));
  }

  auto get_point_id = [&](const typename EK::Point_3& pt)
  {
    auto insert_res = point_id_map.insert(std::make_pair(pt, soup_points.size()));
    if (insert_res.second)
      soup_points.push_back(to_input(pt));
    return insert_res.first->second;
  };

  for (face_descriptor f : faces(tm))
  {
    int tid = get(tid_map, f);
    if (tid == -1)
    {
      halfedge_descriptor h = halfedge(f, tm);
      soup_triangles.emplace_back(
        CGAL::make_array(get_point_id(to_exact(get(vpm,source(h, tm)))),
                         get_point_id(to_exact(get(vpm,target(h, tm)))),
                         get_point_id(to_exact(get(vpm,target(next(h, tm), tm)))))
      );
    }
    else
    {
      for (const typename EK::Triangle_3& t : triangles[tid])
      {
        soup_triangles.emplace_back(CGAL::make_array(get_point_id(t[0]), get_point_id(t[1]), get_point_id(t[2])));
      }
    }
  }

}
#endif

/**
 * \ingroup PMP_corefinement_grp
 * \link coref_def_subsec autorefines \endlink `tm`. Refines a triangle mesh
 * so that no triangles intersects in their interior.
 * Self-intersection edges will be marked as constrained. If an edge that was marked as
 * constrained is split, its sub-edges will be marked as constrained as well.
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
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *     \cgalParamExtra{If the property map is writable, the indices of the faces of `tm1` and `tm2`
 *                     will be set after the corefinement is done.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{visitor}
 *     \cgalParamDescription{a visitor used to track the creation of new faces}
 *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
 *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
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
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  std::vector<typename GT::Point_3> soup_points;
  std::vector<std::array<std::size_t, 3> > soup_triangles;

  autorefine_soup_output(tm, soup_points, soup_triangles, np);

  clear(tm);
  orient_polygon_soup(soup_points, soup_triangles);
  polygon_soup_to_polygon_mesh(soup_points, soup_triangles, tm);
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_AUTOREFINEMENT_H
