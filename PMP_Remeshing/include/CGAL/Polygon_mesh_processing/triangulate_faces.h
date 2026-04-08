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
*   %Default visitor model for the functions `triangulate_face()` and `triangulate_faces()`, model of `PMPTriangulateFaceVisitor`.
*   \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
*   All its functions have an empty body. This class can be used as a
*   base class if only some of the functions of the concept require to be
*   overridden.
*/
template<class PolygonMesh>
struct Default_visitor
  : public Hole_filling::Default_visitor
{
  /// face descriptor type
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  /// vertex descriptor type
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  void before_subface_creations(face_descriptor /*f_old*/) {}
  void after_subface_creations() {}
  void after_subface_created(face_descriptor /*f_new*/) {}
  constexpr
  bool accept_face(face_descriptor /*f*/, vertex_descriptor /*v0*/, vertex_descriptor /*v1*/, vertex_descriptor /*v2*/) const
  {
    return true;
  }
};

template <class PolygonMesh, class Triangulation_visitor>
struct Visitor_wrapper
  : public Triangulation_visitor
{
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  face_descriptor f;
  std::vector<vertex_descriptor> verts;

  bool accept_triangle(int i0, int i1, int i2) const
  {
    if (!static_cast<const Triangulation_visitor*>(this)->accept_triangle(i0, i1, i2))
      return false;
    return static_cast<const Triangulation_visitor*>(this)->accept_face(f, verts[i0], verts[i1], verts[i2]);
  }

  Visitor_wrapper(PolygonMesh& pmesh, face_descriptor f, Triangulation_visitor& tvisitor)
    : Triangulation_visitor(tvisitor), f(f)
  {
    auto h = halfedge(f, pmesh), hstart=h;
    do{
      verts.push_back(target(h, pmesh));
      h = next(h, pmesh);
    }
    while(h!=hstart);
  }
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
  template <typename Point,
            typename Visitor,
            typename NamedParameters>
  bool triangulate_face_with_hole_filling(face_descriptor f,
                                          PolygonMesh& pmesh,
                                          const std::vector<vertex_descriptor>& border_vertices,
                                          const std::vector<Point>& hole_points,
                                          Visitor visitor,
                                          const NamedParameters& np)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // use hole filling
    typedef CGAL::Triple<int, int, int> Face_indices;
    std::vector<Face_indices> patch;

    // visitor and bool are "replaced" if they already exist in `np` because we will
    // return the first one seen while unstacking the NP inheritance stack
    auto new_np = np.visitor(visitor)
                    .use_2d_constrained_delaunay_triangulation(true);

    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch), new_np);

    if(patch.empty())
      return false;

    // triangulate the hole
    std::map<std::pair<int, int>, halfedge_descriptor > halfedge_map;
    int i = static_cast<int>(hole_points.size())-1;
    int j = 0;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      halfedge_map[std::make_pair(i, j)] = h;
      i = j;
      ++j;
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
      {
        f = add_face(pmesh);
        visitor.after_subface_created(f);
      }

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

    using Point = typename boost::property_traits<VPM>::value_type;


    std::vector<Point> hole_points;
    std::vector<vertex_descriptor> border_vertices;
    CGAL_assertion(CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size() > 0);
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      vertex_descriptor v = target(h, pmesh);
      hole_points.push_back(get(vpm, v));
      border_vertices.push_back(v);
    }
    if(border_vertices.size() <= 3)
      return true;

    if(border_vertices.size() == 4)
    {
      std::array<halfedge_descriptor,4> hverts;
      hverts[0] = halfedge(f, pmesh);
      hverts[1] = next(hverts[0], pmesh);
      hverts[2] = next(hverts[1], pmesh);
      hverts[3] = next(hverts[2], pmesh);
      auto verts = make_array(target(hverts[0],pmesh), target(hverts[1],pmesh),
                              target(hverts[2],pmesh), target(hverts[3],pmesh));

      if ( !visitor.accept_face(f,verts[0],verts[1],verts[3]) ||
           !visitor.accept_face(f,verts[1],verts[2],verts[3]) )
      {
        if ( !visitor.accept_face(f,verts[0],verts[1],verts[2]) ||
             !visitor.accept_face(f,verts[0],verts[2],verts[3]) )
          return false;
        visitor.before_subface_creations(f);
        halfedge_descriptor res = CGAL::Euler::split_face(hverts[0], hverts[2], pmesh);
        visitor.after_subface_created(face(res, pmesh));
        visitor.after_subface_created(face(opposite(res, pmesh), pmesh));
        visitor.after_subface_creations();
        return true;
      }
      if ( !visitor.accept_face(f,verts[0],verts[1],verts[2]) ||
           !visitor.accept_face(f,verts[0],verts[2],verts[3]) )
      {
        visitor.before_subface_creations(f);
        halfedge_descriptor res = CGAL::Euler::split_face(hverts[1], hverts[3], pmesh);
        visitor.after_subface_created(face(res, pmesh));
        visitor.after_subface_created(face(opposite(res, pmesh), pmesh));
        visitor.after_subface_creations();
        return true;
      }

      Point_ref p0 = get(vpm, target(hverts[0], pmesh));
      Point_ref p1 = get(vpm, target(hverts[1], pmesh));
      Point_ref p2 = get(vpm, target(hverts[2], pmesh));
      Point_ref p3 = get(vpm, target(hverts[3], pmesh));

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

      typename Traits::Vector_3 p0p1=p1-p0;
      typename Traits::Vector_3 p1p2=p2-p1;
      typename Traits::Vector_3 p2p3=p3-p2;
      typename Traits::Vector_3 p0p3=p3-p0;

      const FT delta1 = cross_product(p1p2, p2p3) * cross_product(-p0p3, p0p1);
      const FT delta2 = cross_product(p0p1, -p1p2) * cross_product(p2p3, p0p3);

      halfedge_descriptor res = boost::graph_traits<PolygonMesh>::null_halfedge();

      if (delta1!=delta2)
        res = (delta2>delta1)
            ?  CGAL::Euler::split_face(hverts[0], hverts[2], pmesh)
            :  CGAL::Euler::split_face(hverts[1], hverts[3], pmesh);
      else
      {
        halfedge_descriptor m =
          *(std::min_element)(hverts.begin(), hverts.end(),
                              [&pmesh,vpm](halfedge_descriptor v0, halfedge_descriptor v1)
                              {
                                return lexicographically_xyz_smaller(get(vpm, target(v0, pmesh)),
                                                                     get(vpm, target(v1, pmesh)));
                              });
        res = (m==hverts[0] || m==hverts[2])
            ? CGAL::Euler::split_face(hverts[0], hverts[2], pmesh)
            : CGAL::Euler::split_face(hverts[1], hverts[3], pmesh);
      }

      visitor.after_subface_created(face(res, pmesh));
      visitor.after_subface_created(face(opposite(res, pmesh), pmesh));

      visitor.after_subface_creations();

      return true;
    }

    Triangulate_faces::Visitor_wrapper<PolygonMesh,Visitor> Visitor_wrapper(pmesh, f, visitor);

    return triangulate_face_with_hole_filling(f, pmesh, border_vertices, hole_points, Visitor_wrapper, np);
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
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`,
*                    with `PMPTriangulateFaceVisitor::face_descriptor` and `PMPTriangulateFaceVisitor::vertex_descriptor` being
*                    the corresponding descriptor types of `booost::graph_traits<PolygonMesh>`}
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
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`,
*                    with `PMPTriangulateFaceVisitor::face_descriptor` and `PMPTriangulateFaceVisitor::vertex_descriptor` being
*                    the corresponding descriptor types of `booost::graph_traits<PolygonMesh>`}
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
bool triangulate_faces(const FaceRange& face_range,
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
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor` and `PMPHolefillingVisitor`,
*                    with `PMPTriangulateFaceVisitor::face_descriptor` and `PMPTriangulateFaceVisitor::vertex_descriptor` being `std::size_t`.}
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
*   %Default visitor for the function `triangulate_polygons()`, model of `PMPTriangulateFaceVisitor`.
*   All its functions have an empty body. This class can be used as a
*   base class if only some of the functions of the concept require to be
*   overridden.
*/
struct Default_visitor
  : public Hole_filling::Default_visitor
{
  /// face id type, referring to the position of the polygon in the input range
  using face_descriptor = std::size_t;
  /// vertex id type
  using vertex_descriptor = std::size_t;

  void before_subface_creations(face_descriptor) {}

  void after_subface_created(face_descriptor) {}

  void after_subface_creations() {}

  constexpr
  bool accept_face(face_descriptor /*f*/, vertex_descriptor /*v0*/, vertex_descriptor /*v1*/, vertex_descriptor /*v2*/) const
  {
    return true;
  }
};

} // namespace Triangulate_polygons

namespace internal {

template <class Triangulation_visitor, class Polygon>
struct Visitor_wrapper
  : public Triangulation_visitor
{
  std::size_t poly_id;
  const Polygon& pid_map;

  bool accept_triangle(int i0, int i1, int i2) const
  {
    if (!static_cast<const Triangulation_visitor*>(this)->accept_triangle(i0, i1, i2))
      return false;
    return static_cast<const Triangulation_visitor*>(this)->accept_face(poly_id, pid_map[i0], pid_map[i1], pid_map[i2]);
  }

  Visitor_wrapper(std::size_t poly_id,  const Polygon& pid_map, Triangulation_visitor& tvisitor)
    : Triangulation_visitor(tvisitor), poly_id(poly_id), pid_map(pid_map)
  {}
};

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
                                             std::size_t poly_id,
                                             PolygonRange& triangulated_polygons, // output
                                             PMap pm,
                                             Visitor& visitor,
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
    Visitor_wrapper new_visitor(poly_id, polygon, visitor);
    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch), np.visitor(new_visitor));

    if(patch.empty())
      return false;

    for(const Face_indices& triangle : patch)
    {
      triangulated_polygons.push_back({hole_points_indices[triangle.first],
                                       hole_points_indices[triangle.second],
                                       hole_points_indices[triangle.third]});
      visitor.after_subface_created(triangulated_polygons.size()-1);
    }
    return true;
  }

public:
  template <typename Polygon,
            typename PointRange,
            typename PolygonRange,
            typename Visitor,
            typename NamedParameters>
  bool operator()(const Polygon& polygon,
                  const PointRange& points,
                  std::size_t poly_id,
                  PolygonRange& triangulated_polygons,
                  Visitor& visitor,
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
    using PID = typename std::iterator_traits<typename Polygon::const_iterator>::value_type;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    PMap pm = choose_parameter<PMap>(get_parameter(np, internal_np::point_map));
    Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

    typename Traits::Construct_cross_product_vector_3 cross_product =
      traits.construct_cross_product_vector_3_object();

    const std::size_t original_size = polygon.size();
    if(original_size == 4)
    {
      Point_ref p0 = get(pm, points[polygon[0]]);
      Point_ref p1 = get(pm, points[polygon[1]]);
      Point_ref p2 = get(pm, points[polygon[2]]);
      Point_ref p3 = get(pm, points[polygon[3]]);

      if (!visitor.accept_face(poly_id, polygon[0], polygon[1], polygon[3]) ||
          !visitor.accept_face(poly_id, polygon[0], polygon[2], polygon[3]))
      {
        if (!visitor.accept_face(poly_id, polygon[0], polygon[1], polygon[2]) ||
            !visitor.accept_face(poly_id, polygon[0], polygon[2], polygon[3]))
        {
          return false;
        }
        triangulated_polygons.push_back({polygon[0], polygon[1], polygon[2]});
        triangulated_polygons.push_back({polygon[0], polygon[2], polygon[3]});
        visitor.after_subface_created(triangulated_polygons.size()-2);
        visitor.after_subface_created(triangulated_polygons.size()-1);
        return true;
      }

      if (!visitor.accept_face(poly_id, polygon[0], polygon[1], polygon[2]) ||
          !visitor.accept_face(poly_id, polygon[0], polygon[2], polygon[3]))
      {
        triangulated_polygons.push_back({polygon[0], polygon[1], polygon[3]});
        triangulated_polygons.push_back({polygon[1], polygon[2], polygon[3]});
        visitor.after_subface_created(triangulated_polygons.size()-2);
        visitor.after_subface_created(triangulated_polygons.size()-1);
        return true;
      }
      /* Chooses the diagonal that will split the quad in two triangles that maximize
       * the scalar product of the un-normalized normals of the two triangles.
       * The lengths of the un-normalized normals (computed using cross-products of two vectors)
       *  are proportional to the area of the triangles.
       * Maximize the scalar product of the two normals will avoid skinny triangles,
       * and will also taken into account the cosine of the angle between the two normals.
       * In particular, if the two triangles are oriented in different directions,
       * the scalar product will be negative.
       */
      typename Traits::Vector_3 p0p1=p1-p0;
      typename Traits::Vector_3 p1p2=p2-p1;
      typename Traits::Vector_3 p2p3=p3-p2;
      typename Traits::Vector_3 p0p3=p3-p0;

      const FT delta1 = cross_product(p1p2, p2p3) * cross_product(-p0p3, p0p1);
      const FT delta2 = cross_product(p0p1, -p1p2) * cross_product(p2p3, p0p3);
      if (delta1!=delta2)
      {
        if(delta2 > delta1)
        {
          triangulated_polygons.push_back({polygon[0], polygon[1], polygon[2]});
          triangulated_polygons.push_back({polygon[0], polygon[2], polygon[3]});
        }
        else
        {
          triangulated_polygons.push_back({polygon[0], polygon[1], polygon[3]});
          triangulated_polygons.push_back({polygon[1], polygon[2], polygon[3]});
        }
      }
      else
      {
        PID mid =
          *(std::min_element)(polygon.begin(), polygon.end(),
                              [&points,pm](PID id1 , PID id2)
                              {
                                return lexicographically_xyz_smaller(get(pm, points[id1]),
                                                                     get(pm, points[id2]));
                              });
        if (mid==0|| mid==2)
        {
          triangulated_polygons.push_back({polygon[0], polygon[1], polygon[2]});
          triangulated_polygons.push_back({polygon[0], polygon[2], polygon[3]});
        }
        else
        {
          triangulated_polygons.push_back({polygon[0], polygon[1], polygon[3]});
          triangulated_polygons.push_back({polygon[1], polygon[2], polygon[3]});
        }
      }
      visitor.after_subface_created(triangulated_polygons.size()-2);
      visitor.after_subface_created(triangulated_polygons.size()-1);

      return true;
    }

    return triangulate_polygon_with_hole_filling(polygon, points, poly_id, triangulated_polygons, pm, visitor, np);
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

  auto visitor = parameters::choose_parameter<Triangulate_polygons::Default_visitor>(
                   parameters::get_parameter(np, internal_np::visitor));

  bool success = true;

  internal::Triangulate_polygon_soup_modifier modifier;
  for(std::size_t poly_id=0; poly_id<polygons.size(); ++poly_id)
  {
    Polygon& polygon = polygons[poly_id];

    visitor.before_subface_creations(poly_id);
    if(polygon.size() <= 3)
    {
      triangulated_polygons.push_back(std::move(polygon));
      visitor.after_subface_created(triangulated_polygons.size()-1);
    }
    else
      if(!modifier(polygon, points, poly_id, triangulated_polygons, visitor, np))
      {
        success = false;
        triangulated_polygons.push_back(std::move(polygon));
        visitor.after_subface_created(triangulated_polygons.size()-1);
      }

    visitor.after_subface_creations();
  }

  std::swap(polygons, triangulated_polygons);

  return success;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
