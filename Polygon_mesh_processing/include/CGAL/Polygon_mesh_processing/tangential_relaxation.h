// Copyright (c) 2015-2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Ivan Paden


#ifndef CGAL_POLYGON_MESH_PROCESSING_TANGENTIAL_RELAXATION_H
#define CGAL_POLYGON_MESH_PROCESSING_TANGENTIAL_RELAXATION_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/property_map.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <tuple>
#include <type_traits>
#include <unordered_map>

namespace CGAL {
namespace Polygon_mesh_processing {

namespace internal {
struct Allow_all_moves{
  template <class vertex_descriptor, class Point_3>
  constexpr inline bool operator()(vertex_descriptor, Point_3, Point_3) const
  {
    return true;
  }
};
} // internal namespace


/*!
* \ingroup PMP_meshing_grp
* applies an iterative area-based tangential smoothing to the given range of vertices.
* Each vertex `v` is relocated to its gravity-weighted centroid, and the relocation vector
* is projected back to the tangent plane to the surface at `v`, iteratively.
* The connectivity remains unchanged.
*
* @tparam TriangleMesh model of `FaceGraph` and `VertexListGraph`.
*         The descriptor types `boost::graph_traits<TriangleMesh>::%face_descriptor`
*         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
* @tparam VertexRange range of `boost::graph_traits<TriangleMesh>::%vertex_descriptor`,
*         model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param vertices the range of vertices which will be relocated by relaxation
* @param tm the triangle mesh to which `vertices` belong
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `tm`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the `Point_3` type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex `Point_3` type.}
*     \cgalParamExtra{Exact constructions kernels are not supported by this function.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of smoothing iterations}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{edge_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm`.
*                           The endpoints of a constrained edge cannot be moved by relaxation.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no edges are constrained}
*     \cgalParamExtra{Boundary edges are always considered as constrained edges.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of `tm`.
*                           A constrained vertex cannot be modified during relaxation.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no vertices are constrained}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{relax_constraints}
*     \cgalParamDescription{If `true`, the end vertices of the edges set as constrained
*                           in `edge_is_constrained_map` and boundary edges move along the
*                           constrained polylines they belong to.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{allow_move_functor}
*     \cgalParamDescription{A function object used to determinate if a vertex move should be allowed or not}
*     \cgalParamType{Unary functor that provides `bool operator()(vertex_descriptor v, Point_3 src, Point_3 tgt)` returning `true`
*                    if the vertex `v` can be moved from `src` to `tgt`; `Point_3` being the value type of the vertex point map }
*     \cgalParamDefault{If not provided, all moves are allowed.}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* \todo check if it should really be a triangle mesh or if a polygon mesh is fine
*/
template <typename VertexRange, class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
void tangential_relaxation(const VertexRange& vertices,
                           TriangleMesh& tm,
                           const NamedParameters& np = parameters::default_values())
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT gt = choose_parameter(get_parameter(np, internal_np::geom_traits), GT());

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPMap;
  VPMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_property_map(vertex_point, tm));

  typedef Static_boolean_property_map<edge_descriptor, false> Default_ECM;
  typedef typename internal_np::Lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      Static_boolean_property_map<edge_descriptor, false> // default (no constraint)
    > ::type ECM;
  ECM ecm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                             Default_ECM());

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      Static_boolean_property_map<vertex_descriptor, false> // default (no constraint)
    > ::type VCM;
  VCM vcm = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                             Static_boolean_property_map<vertex_descriptor, false>());

  const bool relax_constraints = choose_parameter(get_parameter(np, internal_np::relax_constraints), false);
  const unsigned int nb_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);

  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::Point_3 Point_3;

  auto check_normals = [&](vertex_descriptor v)
  {
    bool first_run = true;
    Vector_3 prev = NULL_VECTOR, first = NULL_VECTOR;
    halfedge_descriptor first_h = boost::graph_traits<TriangleMesh>::null_halfedge();
    for (halfedge_descriptor hd : CGAL::halfedges_around_target(v, tm))
    {
      if (is_border(hd, tm)) continue;

      Vector_3 n = compute_face_normal(face(hd, tm), tm, np);
      if (n == CGAL::NULL_VECTOR) //for degenerate faces
        continue;

      if (first_run)
      {
        first_run = false;
        first = n;
        first_h = hd;
      }
      else
      {
        if (!get(ecm, edge(hd, tm)))
          if (to_double(n * prev) <= 0)
            return false;
      }
      prev = n;
    }

    if (first_run)
      return true; //vertex incident only to degenerate faces

    if (!get(ecm, edge(first_h, tm)))
      if (to_double(first * prev) <= 0)
        return false;

    return true;
  };

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::allow_move_functor_t,
      NamedParameters,
      internal::Allow_all_moves// default
    > ::type Shall_move;
  Shall_move shall_move = choose_parameter(get_parameter(np, internal_np::allow_move_functor),
                                           internal::Allow_all_moves());

  for (unsigned int nit = 0; nit < nb_iterations; ++nit)
  {
#ifdef CGAL_PMP_TANGENTIAL_RELAXATION_VERBOSE
    std::cout << "\r\t(Tangential relaxation iteration " << (nit + 1) << " / ";
    std::cout << nb_iterations << ") ";
    std::cout.flush();
#endif

    typedef std::tuple<vertex_descriptor, Vector_3, Point_3> VNP;
    std::vector< VNP > barycenters;
    auto gt_barycenter = gt.construct_barycenter_3_object();
    auto gt_project = gt.construct_projected_point_3_object();

    // at each vertex, compute vertex normal
    std::unordered_map<vertex_descriptor, Vector_3> vnormals;
    compute_vertex_normals(tm, boost::make_assoc_property_map(vnormals), np);

    // at each vertex, compute barycenter of neighbors
    for(vertex_descriptor v : vertices)
    {
      if (get(vcm, v) || CGAL::internal::is_isolated(v, tm))
        continue;

      // collect hedges to detect if we have to handle boundary cases
      std::vector<halfedge_descriptor> interior_hedges, border_halfedges;
      for(halfedge_descriptor h : halfedges_around_target(v, tm))
      {
        if (is_border_edge(h, tm) || get(ecm, edge(h, tm)))
          border_halfedges.push_back(h);
        else
          interior_hedges.push_back(h);
      }

      if (border_halfedges.empty())
      {
        const Vector_3& vn = vnormals.at(v);
        Vector_3 move = CGAL::NULL_VECTOR;
        unsigned int star_size = 0;
        for(halfedge_descriptor h :interior_hedges)
        {
          move = move + Vector_3(get(vpm, v), get(vpm, source(h, tm)));
          ++star_size;
        }
        CGAL_assertion(star_size > 0); //isolated vertices have already been discarded
        move = (1. / static_cast<double>(star_size)) * move;

        barycenters.emplace_back(v, vn, get(vpm, v) + move);
      }
      else
      {
        if (!relax_constraints) continue;
        Vector_3 vn(NULL_VECTOR);

        if (border_halfedges.size() == 2)// corners are constrained
        {
          vertex_descriptor ph0 = source(border_halfedges[0], tm);
          vertex_descriptor ph1 = source(border_halfedges[1], tm);
          double dot = to_double(Vector_3(get(vpm, v), get(vpm, ph0))
                                 * Vector_3(get(vpm, v), get(vpm, ph1)));
          // \todo shouldn't it be an input parameter?
          //check squared cosine is < 0.25 (~120 degrees)
          if (0.25 < dot*dot / ( squared_distance(get(vpm,ph0), get(vpm, v)) *
                                 squared_distance(get(vpm,ph1), get(vpm, v))) )
          {
            typename GT::Point_3 bary = gt_barycenter(get(vpm, ph0), 0.25, get(vpm, ph1), 0.25, get(vpm, v), 0.5);
            // to avoid shrinking of borders, we project back onto the incident segments
            typename GT::Segment_3 s1(get(vpm, ph0), get(vpm,v)),
                                   s2(get(vpm, ph1), get(vpm,v));

            typename GT::Point_3 p1 = gt_project(s1, bary), p2 = gt_project(s2, bary);

            bary = squared_distance(p1, bary)<squared_distance(p2,bary)? p1:p2;
            barycenters.emplace_back(v, vn, bary);
          }


        }
      }
    }

    // compute moves
    typedef std::pair<vertex_descriptor, Point_3> VP_pair;
    std::vector< std::pair<vertex_descriptor, Point_3> > new_locations;
    new_locations.reserve(barycenters.size());
    for(const VNP& vnp : barycenters)
    {
      vertex_descriptor v = std::get<0>(vnp);
      const Point_3& pv = get(vpm, v);
      const Vector_3& nv = std::get<1>(vnp);
      const Point_3& qv = std::get<2>(vnp); //barycenter at v

      new_locations.emplace_back(v, qv + (nv * Vector_3(qv, pv)) * nv);
    }

    // perform moves
    for(const VP_pair& vp : new_locations)
    {
      const Point_3 initial_pos = get(vpm, vp.first); // make a copy on purpose
      const Vector_3 move(initial_pos, vp.second);

      put(vpm, vp.first, vp.second);

      //check that no inversion happened
      double frac = 1.;
      while (frac > 0.03 //5 attempts maximum
          && (   !check_normals(vp.first)
              || !shall_move(vp.first, initial_pos, get(vpm, vp.first)))) //if a face has been inverted
      {
        frac = 0.5 * frac;
        put(vpm, vp.first, initial_pos + frac * move);//shorten the move by 2
      }
      if (frac <= 0.02)
        put(vpm, vp.first, initial_pos);//cancel move
    }
  }//end for loop (nit == nb_iterations)

#ifdef CGAL_PMP_TANGENTIAL_RELAXATION_VERBOSE
  std::cout << "\rTangential relaxation : "
    << nb_iterations << " iterations done." << std::endl;
#endif
}

/*!
* \ingroup PMP_meshing_grp
* applies an iterative area-based tangential smoothing to the given range of vertices based on the
* underlying sizing field.
* Each vertex `v` is relocated to its weighted centroid, where weights depend on the area of the
* adjacent triangle and its averaged vertex sizing values.
* The relocation vector is projected back to the tangent plane to the surface at `v`, iteratively.
* The connectivity remains unchanged.
*
* @tparam TriangleMesh model of `FaceGraph` and `VertexListGraph`.
*         The descriptor types `boost::graph_traits<TriangleMesh>::%face_descriptor`
*         and `boost::graph_traits<TriangleMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
* @tparam VertexRange range of `boost::graph_traits<TriangleMesh>::%vertex_descriptor`,
*         model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam SizingFunction model of `PMPSizingField`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param vertices the range of vertices which will be relocated by relaxation
* @param tm the triangle mesh to which `vertices` belong
* @param sizing a map containing sizing field for individual vertices.
*        Used to derive smoothing weights.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `tm`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the `Point_3` type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex `Point_3` type.}
*     \cgalParamExtra{Exact constructions kernels are not supported by this function.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of smoothing iterations}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{edge_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm`.
*                           The endpoints of a constrained edge cannot be moved by relaxation, unless `relax_constraints` is `true`.
*                           The endpoints of a constrained polyline can never be moved by relaxation.
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no edges are constrained}
*     \cgalParamExtra{Boundary edges are always considered as constrained edges.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of `tm`.
*                           A constrained vertex cannot be modified during relaxation.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no vertices are constrained}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{relax_constraints}
*     \cgalParamDescription{If `true`, the end vertices of the edges set as constrained
*                           in `edge_is_constrained_map` and boundary edges move along the
*                           constrained polylines they belong to.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{allow_move_functor}
*     \cgalParamDescription{A function object used to determinate if a vertex move should be allowed or not}
*     \cgalParamType{Unary functor that provides `bool operator()(vertex_descriptor v, Point_3 src, Point_3 tgt)` returning `true`
*                    if the vertex `v` can be moved from `src` to `tgt`; `Point_3` being the value type of the vertex point map }
*     \cgalParamDefault{If not provided, all moves are allowed.}
*   \cgalParamNEnd
*
* \cgalNamedParamsEnd
*
* \todo check if it should really be a triangle mesh or if a polygon mesh is fine
*/
template <typename VertexRange,
          class TriangleMesh,
          class SizingFunction,
          class NamedParameters = parameters::Default_named_parameters>
void tangential_relaxation_with_sizing(const VertexRange& vertices,
                                       TriangleMesh& tm,
                                       const SizingFunction& sizing,
                                       const NamedParameters& np = parameters::default_values())
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT gt = choose_parameter(get_parameter(np, internal_np::geom_traits), GT());

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPMap;
  VPMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
    get_property_map(vertex_point, tm));

  typedef Static_boolean_property_map<edge_descriptor, false> Default_ECM;
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Static_boolean_property_map<edge_descriptor, false> // default (no constraint)
  > ::type ECM;
  ECM ecm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
    Default_ECM());

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::vertex_is_constrained_t,
    NamedParameters,
    Static_boolean_property_map<vertex_descriptor, false> // default (no constraint)
  > ::type VCM;
  VCM vcm = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
    Static_boolean_property_map<vertex_descriptor, false>());

  const bool relax_constraints = choose_parameter(get_parameter(np, internal_np::relax_constraints), false);
  const unsigned int nb_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);

  typedef typename GT::Vector_3 Vector_3;
  typedef typename GT::Point_3 Point_3;

  auto check_normals = [&](vertex_descriptor v)
  {
    bool first_run = true;
    Vector_3 prev = NULL_VECTOR, first = NULL_VECTOR;
    halfedge_descriptor first_h = boost::graph_traits<TriangleMesh>::null_halfedge();
    for (halfedge_descriptor hd : CGAL::halfedges_around_target(v, tm))
    {
      if (is_border(hd, tm)) continue;

      Vector_3 n = compute_face_normal(face(hd, tm), tm, np);
      if (n == CGAL::NULL_VECTOR) //for degenerate faces
        continue;

      if (first_run)
      {
        first_run = false;
        first = n;
        first_h = hd;
      }
      else
      {
        if (!get(ecm, edge(hd, tm)))
          if (to_double(n * prev) <= 0)
            return false;
      }
      prev = n;
    }

    if (first_run)
      return true; //vertex incident only to degenerate faces

    if (!get(ecm, edge(first_h, tm)))
      if (to_double(first * prev) <= 0)
        return false;

    return true;
  };

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::allow_move_functor_t,
    NamedParameters,
    internal::Allow_all_moves// default
  > ::type Shall_move;
  Shall_move shall_move = choose_parameter(get_parameter(np, internal_np::allow_move_functor),
    internal::Allow_all_moves());

  for (unsigned int nit = 0; nit < nb_iterations; ++nit)
  {
#ifdef CGAL_PMP_TANGENTIAL_RELAXATION_VERBOSE
    std::cout << "\r\t(Tangential relaxation iteration " << (nit + 1) << " / ";
    std::cout << nb_iterations << ") ";
    std::cout.flush();
#endif

    typedef std::tuple<vertex_descriptor, Vector_3, Point_3> VNP;
    std::vector< VNP > barycenters;
    auto gt_barycenter = gt.construct_barycenter_3_object();
    auto gt_centroid = gt.construct_centroid_3_object();
    auto gt_area = gt.compute_area_3_object();

    // at each vertex, compute vertex normal
    std::unordered_map<vertex_descriptor, Vector_3> vnormals;
    compute_vertex_normals(tm, boost::make_assoc_property_map(vnormals), np);

    // at each vertex, compute centroids of neighbouring faces weighted by
    // area and sizing field
    for(vertex_descriptor v : vertices)
    {
      if (get(vcm, v) || CGAL::internal::is_isolated(v, tm))
        continue;

      // collect hedges to detect if we have to handle boundary cases
      std::vector<halfedge_descriptor> interior_hedges, border_halfedges;
      for(halfedge_descriptor h : halfedges_around_target(v, tm))
      {
        if (is_border_edge(h, tm) || get(ecm, edge(h, tm)))
          border_halfedges.push_back(h);
        else
          interior_hedges.push_back(h);
      }

      if (border_halfedges.empty())
      {
        const Vector_3& vn = vnormals.at(v);
        Vector_3 move = CGAL::NULL_VECTOR;
        double weight = 0;
        for(halfedge_descriptor h :interior_hedges)
        {
          // calculate weight
          // need v, v1 and v2
          const vertex_descriptor v1 = target(next(h, tm), tm);
          const vertex_descriptor v2 = source(h, tm);

          const double tri_area = gt_area(get(vpm, v), get(vpm, v1), get(vpm, v2));
          const double face_weight = tri_area
                                   / (1. / 3. * (sizing.get_sizing(v) + sizing.get_sizing(v1) + sizing.get_sizing(v2)));
          weight += face_weight;

          const Point_3 centroid = gt_centroid(get(vpm, v), get(vpm, v1), get(vpm, v2));
          move = move + Vector_3(get(vpm, v), centroid) * face_weight;
        }
        move = move / weight; //todo ip: what if weight ends up being close to 0?

        barycenters.emplace_back(v, vn, get(vpm, v) + move);
      }
      else
      {
        if (!relax_constraints) continue;
        Vector_3 vn(NULL_VECTOR);

        if (border_halfedges.size() == 2)// corners are constrained
        {
          vertex_descriptor ph0 = source(border_halfedges[0], tm);
          vertex_descriptor ph1 = source(border_halfedges[1], tm);
          double dot = to_double(Vector_3(get(vpm, v), get(vpm, ph0))
                                 * Vector_3(get(vpm, v), get(vpm, ph1)));
          // \todo shouldn't it be an input parameter?
          //check squared cosine is < 0.25 (~120 degrees)
          if (0.25 < dot*dot / ( squared_distance(get(vpm,ph0), get(vpm, v)) *
                                 squared_distance(get(vpm,ph1), get(vpm, v))) )
            barycenters.emplace_back(v, vn,
              gt_barycenter(get(vpm, ph0), 0.25, get(vpm, ph1), 0.25, get(vpm, v), 0.5));
        }
      }
    }

    // compute moves
    typedef std::pair<vertex_descriptor, Point_3> VP_pair;
    std::vector< std::pair<vertex_descriptor, Point_3> > new_locations;
    new_locations.reserve(barycenters.size());
    for(const VNP& vnp : barycenters)
    {
      vertex_descriptor v = std::get<0>(vnp);
      const Point_3& pv = get(vpm, v);
      const Vector_3& nv = std::get<1>(vnp);
      const Point_3& qv = std::get<2>(vnp); //barycenter at v

      new_locations.emplace_back(v, qv + (nv * Vector_3(qv, pv)) * nv);
    }

    // perform moves
    for(const VP_pair& vp : new_locations)
    {
      const Point_3 initial_pos = get(vpm, vp.first); // make a copy on purpose
      const Vector_3 move(initial_pos, vp.second);

      put(vpm, vp.first, vp.second);

      //check that no inversion happened
      double frac = 1.;
      while (frac > 0.03 //5 attempts maximum
             && (   !check_normals(vp.first)
                    || !shall_move(vp.first, initial_pos, get(vpm, vp.first)))) //if a face has been inverted
      {
        frac = 0.5 * frac;
        put(vpm, vp.first, initial_pos + frac * move);//shorten the move by 2
      }
      if (frac <= 0.02)
        put(vpm, vp.first, initial_pos);//cancel move
    }
  }//end for loop (nit == nb_iterations)

#ifdef CGAL_PMP_TANGENTIAL_RELAXATION_VERBOSE
  std::cout << "\rTangential relaxation : "
    << nb_iterations << " iterations done." << std::endl;
#endif
}

/*!
* \ingroup PMP_meshing_grp
* applies `tangential_relaxation()` to all the vertices of `tm`.
*/
template <class TriangleMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void tangential_relaxation(TriangleMesh& tm, const CGAL_NP_CLASS& np = parameters::default_values())
{
  tangential_relaxation(vertices(tm), tm, np);
}

template <class TriangleMesh,
          typename SizingFunction,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
void tangential_relaxation_with_sizing(TriangleMesh& tm,
                                       const SizingFunction& sizing,
                                       const CGAL_NP_CLASS& np = parameters::default_values())
{
  tangential_relaxation_with_sizing(vertices(tm), tm, sizing, np);
}

} } // CGAL::Polygon_mesh_processing

#endif //CGAL_POLYGON_MESH_PROCESSING_TANGENTIAL_RELAXATION_H
