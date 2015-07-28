// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Jane Tournois


#ifndef CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
#define CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/Origin.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <CGAL/boost/graph/properties.h> 
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL{

namespace Polygon_mesh_processing{

template<typename Point, typename PM, typename VertexPointMap, typename Vector>
void sum_normals(const PM& pmesh,
                 typename boost::graph_traits<PM>::face_descriptor f,
                 VertexPointMap vpmap,
                 Vector& sum)
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor he = halfedge(f, pmesh);
  halfedge_descriptor end = he;
  do
  {
    const Point& prv = get(vpmap, target(prev(he, pmesh), pmesh));
    const Point& curr = get(vpmap, target(he, pmesh));
    const Point& nxt = get(vpmap, target(next(he, pmesh), pmesh));
    Vector n = CGAL::cross_product(nxt - curr, prv - curr);
    sum = sum + n;

    he = next(he, pmesh);
  } while (he != end);
}

/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal to face `f`.
* @tparam PolygonMesh a model of `FaceGraph` that has an internal property map
*         for `CGAL::vertex_point_t`
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param f the face on which the normal is computed
* @param pmesh the polygon mesh containing `f`
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return the computed normal. The return type is a 3D vector type. It is
* either deduced from the `geom_traits` \ref namedparameters if provided,
* or from the geometric traits class deduced from the point property map
* of `pmesh`.
*/
template <typename PolygonMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Vector_3
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3
#endif
compute_face_normal(typename boost::graph_traits<PolygonMesh>::face_descriptor f
                    , const PolygonMesh& pmesh
                    , const NamedParameters& np)
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  using boost::get_param;
  using boost::choose_param;

  Vector normal = CGAL::NULL_VECTOR;
  sum_normals<Point>(pmesh, f
    , choose_param(get_param(np, vertex_point), get(CGAL::vertex_point, pmesh))
    , normal);
 
  return normal / std::sqrt(normal * normal);
}

/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal for all faces of the polygon mesh.
* @tparam PolygonMesh a model of `FaceGraph` that has an internal property map
*         for `CGAL::vertex_point_t`
* @tparam FaceNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
    `Kernel::Vector_3` as value type.
*
* @param pmesh the polygon mesh
* @param fnm the property map in which the normals are written
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
          , typename FaceNormalMap
          , typename NamedParameters>
void
compute_face_normals(const PolygonMesh& pmesh
                   , FaceNormalMap fnm
                   , const NamedParameters& np)
{
  typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type Kernel;

  typename boost::graph_traits<PolygonMesh>::face_descriptor f;
  BOOST_FOREACH(f, faces(pmesh)){
    typename Kernel::Vector_3 vec = compute_face_normal(f, pmesh, np);
    put(fnm, f, vec);
  }
}

/**
* \ingroup PkgPolygonMeshProcessing
* computes the unit normal at vertex `v` as the average of the normals of incident faces.
* @tparam PolygonMesh a model of `FaceGraph` that has an internal property map
*         for `CGAL::vertex_point_t`
*
* @param v the vertex at which the normal is computed
* @param pmesh the polygon mesh containing `v`
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
* @return the computed normal. The return type is a 3D vector type. It is
* either deduced from the `geom_traits` \ref namedparameters if provided,
* or the geometric traits class deduced from the point property map
* of `pmesh`.
*/
template<typename PolygonMesh, typename NamedParameters>
#ifdef DOXYGEN_RUNNING
Vector_3
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Vector_3
#endif
compute_vertex_normal(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                      const PolygonMesh& pmesh,
                      const NamedParameters& np
                      )
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Kernel;

  typedef typename Kernel::Vector_3 Vector;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  Vector normal = CGAL::NULL_VECTOR;
  halfedge_descriptor he = halfedge(v, pmesh);
  halfedge_descriptor end = he;
  do
    {
    if (!is_border(he, pmesh))
    {
      Vector n = compute_face_normal(face(he, pmesh), pmesh, np);
      normal = normal + (n / std::sqrt(n*n));
    }
    he = opposite(next(he, pmesh), pmesh);
  } while (he != end);

  return normal / std::sqrt(normal * normal);
}

/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal for all vertices of the polygon mesh.
* @tparam PolygonMesh a model of `FaceListGraph` that has an internal property map
*         for `CGAL::vertex_point_t`
* @tparam VertexNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
    the return type of `compute_vertex_normal()` as value type.
*
* @param pmesh the polygon mesh
* @param vnm the property map in which the normals are written
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
          , typename VertexNormalMap
          , typename NamedParameters
          >
void
compute_vertex_normals(const PolygonMesh& pmesh
                      , VertexNormalMap vnm
                      , const NamedParameters& np
                      )
{
  typedef typename GetGeomTraits<PolygonMesh,NamedParameters>::type Kernel;

  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v;
  BOOST_FOREACH(v, vertices(pmesh)){
    typename Kernel::Vector_3 vec = compute_vertex_normal(v, pmesh, np);
    put(vnm, v, vec);
  }
}

/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal for all vertices and faces of the polygon mesh.
* @tparam PolygonMesh a model of `FaceListGraph` that has an internal property map
*         for `CGAL::vertex_point_t`
* @tparam VertexNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
    `Kernel::Vector_3` as value type.
* @tparam FaceNormalMap a model of `WritablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
    `Kernel::Vector_3` as value type.
*
* @param pmesh the polygon mesh
* @param vnm the property map in which the vertex normals are written
* @param fnm the property map in which the face normals are written
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
*    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
          , typename VertexNormalMap
          , typename FaceNormalMap
          , typename NamedParameters
          >
void
compute_normals(const PolygonMesh& pmesh
                , VertexNormalMap vnm
                , FaceNormalMap fnm
                , const NamedParameters& np
                )
{
  compute_vertex_normals(pmesh, vnm, np);
  compute_face_normals(pmesh, fnm, np);
}

///\cond SKIP_IN_MANUAL
// compute_vertex_normal overloads
template <typename PolygonMesh>
typename CGAL::Kernel_traits< typename property_map_value<PolygonMesh, CGAL::vertex_point_t>::type>::Kernel::Vector_3
compute_vertex_normal(
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
  const PolygonMesh& pmesh)
{
  return compute_vertex_normal(v, pmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

// compute_vertex_normals overloads
template <typename PolygonMesh, typename VertexNormalMap>
void
compute_vertex_normals(const PolygonMesh& pmesh,
                      VertexNormalMap vnm)
{
  compute_vertex_normals(pmesh, vnm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

// compute_face_normal overload
template <typename PolygonMesh>
typename CGAL::Kernel_traits < typename property_map_value<PolygonMesh, CGAL::vertex_point_t>::type>::Kernel::Vector_3
compute_face_normal(
  typename boost::graph_traits<PolygonMesh>::face_descriptor f,
  const PolygonMesh& pmesh)
{
  return compute_face_normal(f, pmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

// compute_face_normals overload
template <typename PolygonMesh, typename FaceNormalMap>
void
compute_face_normals(const PolygonMesh& pmesh, FaceNormalMap fnm)
{
  compute_face_normals(pmesh, fnm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

// compute_normals overload
template <typename PolygonMesh, typename VertexNormalMap, typename FaceNormalMap>
void
compute_normals(const PolygonMesh& pmesh,
                VertexNormalMap vnm,
                FaceNormalMap fnm)
{
  compute_normals(pmesh, vnm, fnm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

}

} // end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
