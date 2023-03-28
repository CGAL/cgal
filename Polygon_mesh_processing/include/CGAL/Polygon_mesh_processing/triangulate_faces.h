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

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/range/size.hpp>

#include <queue>
#include <vector>
#include <utility>
#include <CGAL/array.h>

#ifdef CGAL_TRIANGULATE_FACES_DO_NOT_USE_CDT2
# ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
#   define CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
# endif
#endif

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
{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;

  void before_subface_creations(face_descriptor /*f_old*/) {}
  void after_subface_creations() {}
  void after_subface_created(face_descriptor /*f_new*/) {}
};

} //end namespace Triangulate_faces

namespace internal {

template <class PM
          , typename VertexPointMap
          , typename Kernel
          , typename Visitor>
class Triangulate_polygon_mesh_modifier
{
  typedef Kernel Traits;

  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;
  typedef typename Kernel::Point_3 Point;

  struct Face_info {
    typename boost::graph_traits<PM>::halfedge_descriptor e[3];
    bool is_external;
  };

  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
  VertexPointMap _vpmap;
  Traits _traits;

public:
  Triangulate_polygon_mesh_modifier(VertexPointMap vpmap, const Traits& traits = Traits())
    : _vpmap(vpmap), _traits(traits)
  {
  }

  template <class Face_handle>
  bool is_external(Face_handle fh) const {
    return fh->info().is_external;
  }

  bool triangulate_face(face_descriptor f, PM& pmesh, bool use_cdt, Visitor visitor)
  {
    typedef typename Traits::FT FT;

    typename Traits::Vector_3 normal =
      Polygon_mesh_processing::compute_face_normal(
        f, pmesh, CGAL::parameters::geom_traits(_traits)
                                   .vertex_point_map(_vpmap));

    if(normal == typename Traits::Vector_3(0,0,0))
      return false;

    std::size_t original_size = CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size();
    if(original_size == 4)
    {
      halfedge_descriptor v0, v1, v2, v3;
      v0 = halfedge(f, pmesh);
      Point_ref p0 = get(_vpmap, target(v0, pmesh));
      v1 = next(v0, pmesh);
      Point_ref p1 = get(_vpmap, target(v1, pmesh));
      v2 = next(v1, pmesh);
      Point_ref p2 = get(_vpmap, target(v2, pmesh));
      v3 = next(v2, pmesh);
      Point_ref p3 = get(_vpmap, target(v3, pmesh));

      /* Chooses the diagonal that will split the quad in two triangles that maximize
       * the scalar product of of the un-normalized normals of the two triangles.
       * The lengths of the un-normalized normals (computed using cross-products of two vectors)
       *  are proportional to the area of the triangles.
       * Maximize the scalar product of the two normals will avoid skinny triangles,
       * and will also taken into account the cosine of the angle between the two normals.
       * In particular, if the two triangles are oriented in different directions,
       * the scalar product will be negative.
       */
      FT p1p3 = CGAL::cross_product(p2-p1,p3-p2) * CGAL::cross_product(p0-p3,p1-p0);
      FT p0p2 = CGAL::cross_product(p1-p0,p1-p2) * CGAL::cross_product(p3-p2,p3-p0);
      visitor.before_subface_creations(f);
      halfedge_descriptor res = (p0p2>p1p3)
                              ?  CGAL::Euler::split_face(v0, v2, pmesh)
                              :  CGAL::Euler::split_face(v1, v3, pmesh);

      visitor.after_subface_created(face(res,pmesh));
      visitor.after_subface_created(face(opposite(res,pmesh),pmesh));

      visitor.after_subface_creations();

      return true;
    }

    return triangulate_face_with_hole_filling(f, pmesh, use_cdt, visitor);
  }

  bool triangulate_face_with_hole_filling(face_descriptor f,
                                          PM& pmesh,
                                          const bool use_cdt,
                                          Visitor visitor)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // gather halfedges around the face
    std::vector<Point> hole_points;
    std::vector<vertex_descriptor> border_vertices;
    CGAL_assertion(CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh).size() > 0);
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      vertex_descriptor v = source(h, pmesh);
      hole_points.push_back( get(_vpmap, v) );
      border_vertices.push_back(v);
    }

    // use hole filling
    typedef CGAL::Triple<int, int, int> Face_indices;
    std::vector<Face_indices> patch;
    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch),
                                   parameters::geom_traits(_traits)
                                              .use_2d_constrained_delaunay_triangulation(use_cdt));

    if(patch.empty())
      return false;

    // triangulate the hole
    std::map< std::pair<int, int> , halfedge_descriptor > halfedge_map;
    int i=0;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      int j = std::size_t(i+1) == hole_points.size() ? 0 : i+1;
      halfedge_map[ std::make_pair(i, j) ] = h;
      ++i;
    }

    visitor.before_subface_creations(f);
    bool first = true;
    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(4);
    for(const Face_indices& triangle : patch)
    {
      if (first)
        first=false;
      else
        f=add_face(pmesh);
      visitor.after_subface_created(f);

      std::array<int, 4> indices =
        make_array( triangle.first,
                    triangle.second,
                    triangle.third,
                    triangle.first );
      for (int i=0; i<3; ++i)
      {
        typename std::map< std::pair<int, int> , halfedge_descriptor >::iterator insert_res =
          halfedge_map.insert(
            std::make_pair( std::make_pair(indices[i], indices[i+1]),
                            boost::graph_traits<PM>::null_halfedge() ) ).first;
        if (insert_res->second == boost::graph_traits<PM>::null_halfedge())
        {
          halfedge_descriptor nh = halfedge(add_edge(pmesh), pmesh);
          insert_res->second=nh;
          halfedge_map[std::make_pair(indices[i+1], indices[i])]=opposite(nh, pmesh);
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

  template<typename FaceRange>
  bool operator()(FaceRange face_range, PM& pmesh, bool use_cdt, Visitor visitor)
  {
   bool result = true;
    // One need to store facet handles into a vector, because the list of
    // facets of the polyhedron will be modified during the loop, and
    // that invalidates the range [facets_begin(), facets_end()[.
    std::vector<face_descriptor> facets;
    facets.reserve(std::distance(boost::begin(face_range), boost::end(face_range)));

    //only consider non-triangular faces
    for(face_descriptor fit : face_range)
      if ( next( next( halfedge(fit, pmesh), pmesh), pmesh)
        !=       prev( halfedge(fit, pmesh), pmesh) )
        facets.push_back(fit);

    // Iterates on the vector of face descriptors
    for(face_descriptor f : facets)
    {
      if(!this->triangulate_face(f, pmesh, use_cdt, visitor))
       result = false;
    }
    return result;
  }

  void make_hole(halfedge_descriptor h, PM& pmesh)
  {
    //we are not using Euler::make_hole because it has a precondition
    //that the hole is not made on the boundary of the mesh
    //here we allow making a hole on the boundary, and the pair(s) of
    //halfedges that become border-border are fixed by the connectivity
    //setting made in operator()
    CGAL_assertion(!is_border(h, pmesh));
    face_descriptor fd = face(h, pmesh);

    for(halfedge_descriptor hd : halfedges_around_face(h, pmesh))
    {
      CGAL::internal::set_border(hd, pmesh);
    }
    remove_face(fd, pmesh);
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
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor`}
*     \cgalParamDefault{`Triangulate_faces::Default_visitor<PolygonMesh>`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return `true` if the face has been triangulated.
*
* @see `triangulate_faces()`
*/
template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                      PolygonMesh& pmesh,
                      const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(is_valid_face_descriptor(f, pmesh));

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_property_map(vertex_point, pmesh));

  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;
  Kernel traits = choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits));

  //Option
  bool use_cdt = choose_parameter(get_parameter(np, internal_np::use_2d_constrained_delaunay_triangulation), true);

  typedef typename internal_np::Lookup_named_param_def<
    internal_np::visitor_t,
    NamedParameters,
    Triangulate_faces::Default_visitor<PolygonMesh>//default
  >::type Visitor;
  Visitor visitor = choose_parameter<Visitor>(
                             get_parameter(np, internal_np::visitor),
                             Triangulate_faces::Default_visitor<PolygonMesh>());

  internal::Triangulate_polygon_mesh_modifier<PolygonMesh, VPMap, Kernel, Visitor> modifier(vpmap, traits);
  return modifier.triangulate_face(f, pmesh, use_cdt, visitor);
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
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor`}
*     \cgalParamDefault{`Triangulate_faces::Default_visitor<PolygonMesh>`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*  `\cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return `true` if all the faces have been triangulated.
*
* @see `triangulate_face()`
*/
template <typename FaceRange, typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_faces(FaceRange face_range,
                       PolygonMesh& pmesh,
                       const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;
  Kernel traits = choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits));

  //Option
  bool use_cdt = choose_parameter(get_parameter(np, internal_np::use_2d_constrained_delaunay_triangulation), true);

  typedef typename internal_np::Lookup_named_param_def<
    internal_np::visitor_t,
    NamedParameters,
    Triangulate_faces::Default_visitor<PolygonMesh>//default
  >::type Visitor;
  Visitor visitor = choose_parameter<Visitor>(
                                  get_parameter(np, internal_np::visitor),
                                  Triangulate_faces::Default_visitor<PolygonMesh>());

  internal::Triangulate_polygon_mesh_modifier<PolygonMesh, VPMap, Kernel, Visitor> modifier(vpmap, traits);
  return modifier(face_range, pmesh, use_cdt, visitor);
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
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor`}
*     \cgalParamDefault{`Triangulate_faces::Default_visitor<PolygonMesh>`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return `true` if all the faces have been triangulated.
*
* @see `triangulate_face()`
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

namespace internal {

template <typename Kernel>
class Triangulate_polygon_soup_modifier
{
  using Traits = Kernel;
  using Point = typename Traits::Point_3;
  using Vector = typename Traits::Vector_3;

private:
  Traits _traits;

public:
  Triangulate_polygon_soup_modifier(const Traits& traits = Traits())
    : _traits(traits)
  { }

private:
  template<typename Polygon,
           typename PointRange,
           typename PolygonRange,
           typename PMap,
           typename Visitor>
  bool triangulate_face_with_hole_filling(const Polygon& polygon,
                                          const PointRange& points,
                                          PolygonRange& triangulated_polygons,
                                          PMap pmap,
                                          const bool use_cdt,
                                          Visitor visitor)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    // gather halfedges around the face
    std::vector<Point> hole_points;
    std::vector<std::size_t> hole_points_indices;

    for(std::size_t i : polygon)
    {
      hole_points.push_back(get(pmap, points[i]));
      hole_points_indices.push_back(i);
    }

    // use hole filling
    typedef CGAL::Triple<int, int, int> Face_indices;
    std::vector<Face_indices> patch;
    PMP::triangulate_hole_polyline(hole_points, std::back_inserter(patch),
                                   parameters::geom_traits(_traits)
                                              .use_2d_constrained_delaunay_triangulation(use_cdt));

    if(patch.empty())
    {
      std::cout << "failed hole filling" << std::endl;
      return false;
    }

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

  template <typename Polygon,
            typename PointRange,
            typename PolygonRange,
            typename PMap,
            typename Visitor>
  bool triangulate_face(const Polygon& polygon,
                        const PointRange& points,
                        PolygonRange& triangulated_polygons,
                        PMap pmap,
                        const bool use_cdt,
                        Visitor visitor)
  {
    using FT = typename Traits::FT;
    using Point_ref = typename boost::property_traits<PMap>::reference;

    const std::size_t original_size = polygon.size();
    if(original_size == 4)
    {
      Point_ref p0 = get(pmap, points[polygon[0]]);
      Point_ref p1 = get(pmap, points[polygon[1]]);
      Point_ref p2 = get(pmap, points[polygon[2]]);
      Point_ref p3 = get(pmap, points[polygon[3]]);

      /* Chooses the diagonal that will split the quad in two triangles that maximize
       * the scalar product of of the un-normalized normals of the two triangles.
       * The lengths of the un-normalized normals (computed using cross-products of two vectors)
       *  are proportional to the area of the triangles.
       * Maximize the scalar product of the two normals will avoid skinny triangles,
       * and will also taken into account the cosine of the angle between the two normals.
       * In particular, if the two triangles are oriented in different directions,
       * the scalar product will be negative.
       */
      FT p1p3 = CGAL::cross_product(p2-p1, p3-p2) * CGAL::cross_product(p0-p3, p1-p0);
      FT p0p2 = CGAL::cross_product(p1-p0, p1-p2) * CGAL::cross_product(p3-p2, p3-p0);

      visitor.before_subface_creations(polygon);
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

    return triangulate_face_with_hole_filling(polygon, points, triangulated_polygons, pmap, use_cdt, visitor);
  }

public:
  template<typename PolygonRange,
           typename PointRange,
           typename PMap,
           typename Visitor>
  bool operator()(const PolygonRange& polygons,
                  const PointRange& points,
                  PolygonRange& triangulated_polygons,
                  PMap pmap,
                  const bool use_cdt,
                  Visitor visitor)
  {
    using Polygon = typename boost::range_value<PolygonRange>::type;

    bool result = true;
    triangulated_polygons.reserve(polygons.size());

    for(const Polygon& polygon : polygons)
    {
      if(polygon.size() <= 3)
      {
        triangulated_polygons.push_back(polygon);
        continue;
      }

      if(!triangulate_face(polygon, points, triangulated_polygons, pmap, use_cdt, visitor))
        result = false;
    }

    return result;
  }

}; // class Triangulate_polygon_soup_modifier

} // namespace internal

namespace Triangulate_polygons {
namespace internal {

/** \ingroup PMP_meshing_grp
*   %Default new polygon visitor model of `PMPTriangulateFaceVisitor`.
*   All its functions have an empty body. This class can be used as a
*   base class if only some of the functions of the concept require to be
*   overridden.
*/
struct Default_visitor
{
  template <typename Polygon>
  void before_subface_creations(const Polygon& /*f_old*/) {}

  template <typename Polygon>
  void after_subface_created(const Polygon& /*f_new*/) {}

  void after_subface_creations() {}
};

} // namespace internal
} // namespace Triangulate_polygons

/**
* \ingroup PMP_meshing_grp
*
* triangulates all polygons of a polygon soup. This function depends on the package \ref PkgTriangulation2.
*
* @tparam PointRange a model of `ConstRange`. The value type of its iterator is the point type.
* @tparam PolygonRange a model of the concept `SequenceContainer` and `Swappable`,
*                      whose `value_type` is itself a model of the concepts `SequenceContainer`
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
*     \cgalParamDescription{a visitor that enables to track how polygons are triangulated into triangles}
*     \cgalParamType{a class model of `PMPTriangulateFaceVisitor`}
*     \cgalParamDefault{`Triangulate_polygons::Default_visitor`}
*     \cgalParamExtra{Note that the visitor will be copied, so
*                     it must not have any data member that does not have a reference-like type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return `true` if all the polygons have been triangulated.
*
* @see `triangulate_face()`
*/
template <typename PointRange,
          typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
bool triangulate_polygons(const PointRange& points,
                          PolygonRange& polygons,
                          const NamedParameters& np = parameters::default_values())
{
  using Polygon = typename boost::range_value<PolygonRange>::type;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  //VertexPointMap
  using PMap = typename GetPointMap<PointRange, NamedParameters>::const_type;
  PMap pmap = choose_parameter<PMap>(get_parameter(np, internal_np::point_map));

  //Kernel
  using Point = typename boost::property_traits<PMap>::value_type;
  using Def_Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Kernel = typename internal_np::Lookup_named_param_def<
                            internal_np::geom_traits_t,
                            NamedParameters,
                            Def_Kernel>::type;
  Kernel traits = choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits));

  // Option
  bool use_cdt = choose_parameter(get_parameter(np, internal_np::use_2d_constrained_delaunay_triangulation), true);

  typedef typename internal_np::Lookup_named_param_def<
    internal_np::visitor_t,
    NamedParameters,
    Triangulate_polygons::internal::Default_visitor // default
  >::type Visitor;
  Visitor visitor = choose_parameter<Visitor>(get_parameter(np, internal_np::visitor),
                                              Triangulate_polygons::internal::Default_visitor());

  PolygonRange triangulated_polygons;
  internal::Triangulate_polygon_soup_modifier<Kernel> modifier(traits);
  const bool success = modifier(polygons, points, triangulated_polygons, pmap, use_cdt, visitor);

  std::swap(polygons, triangulated_polygons);

  return success;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIANGULATE_FACES_H
