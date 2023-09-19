// Copyright (c) 2010-2011  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/array.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

#include <boost/range/value_type.hpp>

#include <algorithm>
#include <iterator>
#include <map>
#include <queue>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace Triangulate_faces {

/** \ingroup PMP_meshing_grp
*   %Default new face visitor model of `PMPTriangulateFaceVisitor`.
*   All its functions have an empty body. This class can be used as a
*   base class if only some of the functions of the concept require to be
*   overridden.
*/
template<class PolygonMesh>
struct Default_visitor
  : public Hole_filling::Default_visitor
{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;

  void before_subface_creations(face_descriptor /*f_old*/) {}
  void after_subface_creations() {}
  void after_subface_created(face_descriptor /*f_new*/) {}
};

} // namespace Triangulate_faces

namespace internal {

template <typename PolygonMesh>
class Triangulate_polygon_mesh_modifier
{
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

private:
  template <typename VPM,
            typename Visitor,
            typename NamedParameters>
  bool triangulate_face_with_hole_filling(face_descriptor f,
                                          PolygonMesh& pmesh,
                                          const VPM vpm,
                                          Visitor visitor,
                                          const NamedParameters& np)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    using Point = typename boost::property_traits<VPM>::value_type;

    // gather halfedges around the face
    std::vector<Point> hole_points;
    std::vector<vertex_descriptor> border_vertices;
    CGAL_assertion(CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size() > 0);
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      vertex_descriptor v = source(h, pmesh);
      hole_points.push_back(get(vpm, v));
      border_vertices.push_back(v);
    }

    // use hole filling
    typedef CGAL::Triple<int, int, int> Face_indices;
    std::vector<Face_indices> patch;
    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch),
                                   np.use_2d_constrained_delaunay_triangulation(true));

    if(patch.empty())
      return false;

    // triangulate the hole
    std::map<std::pair<int, int>, halfedge_descriptor > halfedge_map;
    int i = 0;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      int j = std::size_t(i+1) == hole_points.size() ? 0 : i+1;
      halfedge_map[std::make_pair(i, j)] = h;
      ++i;
    }

    visitor.before_subface_creations(f);

    bool first = true;
    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(4);
    for(const Face_indices& triangle : patch)
    {
      if(first)
        first = false;
      else
        f = add_face(pmesh);

      visitor.after_subface_created(f);

      std::array<int, 4> indices = make_array(triangle.first,
                                              triangle.second,
                                              triangle.third,
                                              triangle.first);
      for (int i=0; i<3; ++i)
      {
        typename std::map< std::pair<int, int> , halfedge_descriptor >::iterator insert_res =
          halfedge_map.emplace(std::make_pair(indices[i], indices[i+1]),
                               boost::graph_traits<PolygonMesh>::null_halfedge()).first;
        if(insert_res->second == boost::graph_traits<PolygonMesh>::null_halfedge())
        {
          halfedge_descriptor nh = halfedge(add_edge(pmesh), pmesh);
          insert_res->second = nh;
          halfedge_map[std::make_pair(indices[i+1], indices[i])] = opposite(nh, pmesh);
        }
        hedges.push_back(insert_res->second);
      }

      hedges.push_back(hedges.front());
      for(int i=0; i<3;++i)
      {
        set_next(hedges[i], hedges[i+1], pmesh);
        set_face(hedges[i], f, pmesh);
        set_target(hedges[i], border_vertices[indices[i+1]], pmesh);
      }

      set_halfedge(f, hedges[0], pmesh);
      hedges.clear();
    }

    visitor.after_subface_creations();

    return true;
  }

public:
  template <typename NamedParameters>
  bool operator()(face_descriptor f,
                  PolygonMesh& pmesh,
                  const NamedParameters& np)
  {
    using Traits = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
    using VPM = typename GetVertexPointMap<PolygonMesh, NamedParameters>::type;

    using FT = typename Traits::FT;
    using Point_ref = typename boost::property_traits<VPM>::reference;

    using Visitor = typename internal_np::Lookup_named_param_def<
                                            internal_np::visitor_t,
                                            NamedParameters,
                                            Triangulate_faces::Default_visitor<PolygonMesh> // default
                                          >::type;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    CGAL_precondition(is_valid_face_descriptor(f, pmesh));

    Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));
    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_property_map(vertex_point, pmesh));
    Visitor visitor = choose_parameter<Visitor>(get_parameter(np, internal_np::visitor),
                                                Triangulate_faces::Default_visitor<PolygonMesh>());

    typename Traits::Construct_cross_product_vector_3 cross_product =
      traits.construct_cross_product_vector_3_object();

    typename boost::graph_traits<PolygonMesh>::degree_size_type  original_size = degree(f, pmesh);
    if(original_size <= 3)
      return true;

    if(original_size == 4)
    {
      halfedge_descriptor v0, v1, v2, v3;
      v0 = halfedge(f, pmesh);
      Point_ref p0 = get(vpm, target(v0, pmesh));
      v1 = next(v0, pmesh);
      Point_ref p1 = get(vpm, target(v1, pmesh));
      v2 = next(v1, pmesh);
      Point_ref p2 = get(vpm, target(v2, pmesh));
      v3 = next(v2, pmesh);
      Point_ref p3 = get(vpm, target(v3, pmesh));

      /* Chooses the diagonal that will split the quad in two triangles that maximize
       * the scalar product of the un-normalized normals of the two triangles.
       * The lengths of the un-normalized normals (computed using cross-products of two vectors)
       *  are proportional to the area of the triangles.
       * Maximize the scalar product of the two normals will avoid skinny triangles,
       * and will also taken into account the cosine of the angle between the two normals.
       * In particular, if the two triangles are oriented in different directions,
       * the scalar product will be negative.
       */
      visitor.before_subface_creations(f);

      const FT p1p3 = cross_product(p2-p1, p3-p2) * cross_product(p0-p3, p1-p0);
      const FT p0p2 = cross_product(p1-p0, p1-p2) * cross_product(p3-p2, p3-p0);
      halfedge_descriptor res = (p0p2>p1p3) ?  CGAL::Euler::split_face(v0, v2, pmesh)
                                            :  CGAL::Euler::split_face(v1, v3, pmesh);

      visitor.after_subface_created(face(res, pmesh));
      visitor.after_subface_created(face(opposite(res, pmesh), pmesh));

      visitor.after_subface_creations();

      return true;
    }

    return triangulate_face_with_hole_filling(f, pmesh, vpm, visitor, np);
  }
}; // class Triangulate_polygon_mesh_modifier

} // namespace internal

/**
* \ingroup PMP_meshing_grp
*
* triangulates a single face of a polygon mesh. This function depends on the package \ref PkgTriangulation2.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param f face to be triangulated
* @param pmesh the polygon mesh to which the face to be triangulated belongs
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{visitor}
*     \cgalParamDescription{a visitor that enables to track how faces are triangulated into subfaces}
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`}
*     \cgalParamDefault{`Triangulate_faces::Default_visitor<PolygonMesh>`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* This function calls `CGAL::Polygon_mesh_processing::triangulate_hole_polyline()`.
* Refer to its documentation for its named parameters.
*
* @pre The face `f` is not degenerate.
*
* @return `true` if the face has been triangulated.
*
* @see `triangulate_faces()`
*/
template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                      PolygonMesh& pmesh,
                      const NamedParameters& np = parameters::default_values())
{
  internal::Triangulate_polygon_mesh_modifier<PolygonMesh> modifier;
  return modifier(f, pmesh, np);
}

/**
* \ingroup PMP_meshing_grp
*
* triangulates given faces of a polygon mesh. This function depends on the package \ref PkgTriangulation2.
*
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`, model of `Range`.
*         Its iterator type is `InputIterator`.
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param face_range the range of faces to be triangulated
* @param pmesh the polygon mesh to be triangulated
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{visitor}
*     \cgalParamDescription{a visitor that enables to track how faces are triangulated into subfaces}
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`}
*     \cgalParamDefault{`Triangulate_faces::Default_visitor<PolygonMesh>`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* This function calls `CGAL::Polygon_mesh_processing::triangulate_hole_polyline()` for each face.
* Refer to its documentation for its named parameters.
*
* @pre No face within `face_range` is degenerate.
*
* @return `true` if all the faces have been triangulated.
*
* @see `triangulate_face()`
* @see `triangulate_polygons()`
*/
template <typename FaceRange,
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_faces(FaceRange face_range,
                       PolygonMesh& pmesh,
                       const NamedParameters& np = parameters::default_values())
{
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  bool result = true;

  // One needs to store the facets into a vector, because the list of
  // facets of the polyhedron will be modified during the loop, and
  // that invalidates the range [facets_begin(), facets_end()[.
  std::vector<face_descriptor> facets(std::begin(face_range), std::end(face_range));

  internal::Triangulate_polygon_mesh_modifier<PolygonMesh> modifier;
  for(face_descriptor f : facets)
  {
    if(!modifier(f, pmesh, np))
      result = false;
  }

  return result;
}

/**
* \ingroup PMP_meshing_grp
*
* triangulates all faces of a polygon mesh. This function depends on the package \ref PkgTriangulation2.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param pmesh the polygon mesh to be triangulated
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{visitor}
*     \cgalParamDescription{a visitor that enables to track how faces are triangulated into subfaces}
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`}
*     \cgalParamDefault{`Triangulate_faces::Default_visitor<PolygonMesh>`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* This function calls `CGAL::Polygon_mesh_processing::triangulate_hole_polyline()` on all the faces of the polygon mesh.
* Refer to its documentation for its named parameters.
*
* @pre No face of `pmesh` is degenerate.
*
* @return `true` if all the faces have been triangulated.
*
* @see `triangulate_face()`
* @see `triangulate_polygons()`
*/
template <typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_faces(PolygonMesh& pmesh,
                       const NamedParameters& np = parameters::default_values())
{
  return triangulate_faces(faces(pmesh), pmesh, np);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Polygon Soup

namespace Triangulate_polygons {

/** \ingroup PMP_meshing_grp
*   %Default new polygon visitor model of `PMPTriangulateFaceVisitor`.
*   All its functions have an empty body. This class can be used as a
*   base class if only some of the functions of the concept require to be
*   overridden.
*/
struct Default_visitor
  : public Hole_filling::Default_visitor
{
  template <typename Polygon>
  void before_subface_creations(const Polygon& /*f_old*/) {}

  template <typename Polygon>
  void after_subface_created(const Polygon& /*f_new*/) {}

  void after_subface_creations() {}
};

} // namespace Triangulate_polygons

namespace internal {

class Triangulate_polygon_soup_modifier
{
private:
  template<typename Polygon,
           typename PointRange,
           typename PolygonRange,
           typename PMap,
           typename Visitor,
           typename NamedParameters>
  bool triangulate_polygon_with_hole_filling(const Polygon& polygon,
                                             const PointRange& points,
                                             PolygonRange& triangulated_polygons, // output
                                             PMap pm,
                                             Visitor visitor,
                                             const NamedParameters& np)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    using Point = typename boost::property_traits<PMap>::value_type;

    // gather halfedges around the face
    std::vector<Point> hole_points;
    std::vector<std::size_t> hole_points_indices;

    for(std::size_t i : polygon)
    {
      hole_points.push_back(get(pm, points[i]));
      hole_points_indices.push_back(i);
    }

    // use hole filling
    typedef CGAL::Triple<int, int, int> Face_indices;
    std::vector<Face_indices> patch;
    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch), np);

    if(patch.empty())
      return false;

    visitor.before_subface_creations(polygon);

    for(const Face_indices& triangle : patch)
    {
      triangulated_polygons.push_back({hole_points_indices[triangle.first],
                                       hole_points_indices[triangle.second],
                                       hole_points_indices[triangle.third]});
      visitor.after_subface_created(triangulated_polygons.back());
    }

    visitor.after_subface_creations();
    return true;
  }

public:
  template <typename Polygon,
            typename PointRange,
            typename PolygonRange,
            typename NamedParameters>
  bool operator()(const Polygon& polygon,
                  const PointRange& points,
                  PolygonRange& triangulated_polygons,
                  const NamedParameters& np)
  {
    // PointMap
    using PMap = typename GetPointMap<PointRange, NamedParameters>::const_type;
    using Point_ref = typename boost::property_traits<PMap>::reference;

    // Kernel
    using Point = typename boost::property_traits<PMap>::value_type;
    using Def_Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Traits = typename internal_np::Lookup_named_param_def<
                              internal_np::geom_traits_t,
                              NamedParameters,
                              Def_Kernel>::type;
    using FT = typename Traits::FT;

    // Visitor
    using Visitor = typename internal_np::Lookup_named_param_def<
                                            internal_np::visitor_t,
                                            NamedParameters,
                                            Triangulate_polygons::Default_visitor // default
                                          >::type;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    PMap pm = choose_parameter<PMap>(get_parameter(np, internal_np::point_map));
    Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));
    Visitor visitor = choose_parameter<Visitor>(get_parameter(np, internal_np::visitor),
                                                Triangulate_polygons::Default_visitor());

    typename Traits::Construct_cross_product_vector_3 cross_product =
      traits.construct_cross_product_vector_3_object();

    const std::size_t original_size = polygon.size();
    if(original_size == 4)
    {
      Point_ref p0 = get(pm, points[polygon[0]]);
      Point_ref p1 = get(pm, points[polygon[1]]);
      Point_ref p2 = get(pm, points[polygon[2]]);
      Point_ref p3 = get(pm, points[polygon[3]]);

      /* Chooses the diagonal that will split the quad in two triangles that maximize
       * the scalar product of the un-normalized normals of the two triangles.
       * The lengths of the un-normalized normals (computed using cross-products of two vectors)
       *  are proportional to the area of the triangles.
       * Maximize the scalar product of the two normals will avoid skinny triangles,
       * and will also taken into account the cosine of the angle between the two normals.
       * In particular, if the two triangles are oriented in different directions,
       * the scalar product will be negative.
       */
      visitor.before_subface_creations(polygon);

      const FT p1p3 = cross_product(p2-p1, p3-p2) * cross_product(p0-p3, p1-p0);
      const FT p0p2 = cross_product(p1-p0, p1-p2) * cross_product(p3-p2, p3-p0);
      if(p0p2 > p1p3)
      {
        triangulated_polygons.push_back({polygon[0], polygon[1], polygon[2]});
        triangulated_polygons.push_back({polygon[0], polygon[2], polygon[3]});
      }
      else
      {
        triangulated_polygons.push_back({polygon[0], polygon[1], polygon[3]});
        triangulated_polygons.push_back({polygon[1], polygon[2], polygon[3]});
      }

      visitor.after_subface_created(triangulated_polygons[triangulated_polygons.size()-2]);
      visitor.after_subface_created(triangulated_polygons[triangulated_polygons.size()-1]);

      visitor.after_subface_creations();

      return true;
    }

    return triangulate_polygon_with_hole_filling(polygon, points, triangulated_polygons, pm, visitor, np);
  }
}; // class Triangulate_polygon_soup_modifier

} // namespace internal

/**
* \ingroup PMP_meshing_grp
*
* triangulates all polygons of a polygon soup. This function depends on the package \ref PkgTriangulation2.
*
* @tparam PointRange a model of `ConstRange`. The value type of its iterator is the point type.
* @tparam PolygonRange a model of the concepts `SequenceContainer` and `Swappable`,
*                      whose `value_type` is itself a model of the concept `SequenceContainer`
*                      whose `value_type` is `std::size_t`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param points the point geometry of the soup to be triangulated
* @param polygons the polygons to be triangulated
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the point set `points`}
*     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
*                    of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
*     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{visitor}
*     \cgalParamDescription{a visitor that enables to track how polygons are divided into triangles}
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`}
*     \cgalParamDefault{`Triangulate_polygons::Default_visitor`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* This function calls `CGAL::Polygon_mesh_processing::triangulate_hole_polyline()` for each polygon.
* Refer to its documentation for its named parameters.
*
* @pre No polygon within `polygons` is degenerate.
*
* @return `true` if all the polygons have been triangulated.
*
* @see `triangulate_faces()`
*/
template <typename PointRange,
          typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_polygons(const PointRange& points,
                          PolygonRange& polygons,
                          const NamedParameters& np = parameters::default_values())
{
  using Polygon = typename boost::range_value<PolygonRange>::type;

  PolygonRange triangulated_polygons;
  triangulated_polygons.reserve(polygons.size());

  bool success = true;

  internal::Triangulate_polygon_soup_modifier modifier;
  for(const Polygon& polygon : polygons)
  {
    if(polygon.size() <= 3)
    {
      triangulated_polygons.push_back(polygon);
      continue;
    }

    if(!modifier(polygon, points, triangulated_polygons, np))
      success = false;
  }

  std::swap(polygons, triangulated_polygons);

  return success;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
