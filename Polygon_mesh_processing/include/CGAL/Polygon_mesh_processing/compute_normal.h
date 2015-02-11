// Copyright (c) 2013 INRIA Sophia-Anitpolis (France).
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
// Author(s)     : Pierre Alliez


#ifndef CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
#define CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H

#include <CGAL/boost/graph/helpers.h>

namespace CGAL{

namespace Polygon_mesh_processing{
///\cond SKIP_IN_MANUAL

  template <typename PolygonMesh, typename VertexNormalPropertyMap, typename FaceNormalPropertyMap>
void
compute_normals(const PolygonMesh& pmesh,
                VertexNormalPropertyMap vnpm,
                FaceNormalPropertyMap fnpm)
{  
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  compute_normals(pmesh, vnpm, fnpm, get(vertex_point, pmesh), Kernel());
}


template <typename PolygonMesh, typename VertexNormalPropertyMap, typename FaceNormalPropertyMap, typename PointPropertyMap>
void
compute_normals(const PolygonMesh& pmesh,
                VertexNormalPropertyMap vnpm,
                FaceNormalPropertyMap fnpm,
                PointPropertyMap ppmap)
{  
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  compute_normals(pmesh, vnpm, fnpm, ppmap, Kernel());
}

///\endcond
/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal for all vertices and faces of the polygon mesh.
* @tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the point type of `PolygonMesh`.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam VertexNormalPropertyMap the property map in which the vertex normals are written.
* @tparam FaceNormalPropertyMap the property map in which the face normals are written.
* @tparam PointPropertyMap the property map with the points associated to the vertices.
*
* @param pmesh the polygon mesh
*/
template <typename PolygonMesh
          , typename VertexNormalPropertyMap
          , typename FaceNormalPropertyMap
          , typename PointPropertyMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, vertex_point_t>::type
#endif
          , typename Kernel
#ifdef DOXYGEN_RUNNING
          = typename Kernel_traits<
              typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
          >
void
compute_normals(const PolygonMesh& pmesh
                , VertexNormalPropertyMap vnpm
                , FaceNormalPropertyMap fnpm
                , PointPropertyMap ppmap
#ifdef DOXYGEN_RUNNING
                = get(vertex_point, pmesh)
#endif
                ,const Kernel& k
#ifdef DOXYGEN_RUNNING
                = Kernel()
#endif
                )
{
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v;
  BOOST_FOREACH(v, vertices(pmesh)){
    typename Kernel::Vector_3 vec = compute_vertex_normal(v,pmesh, ppmap, k);
    put(vnpm, v, vec);
  }
  typename boost::graph_traits<PolygonMesh>::face_descriptor f;
  BOOST_FOREACH(f, faces(pmesh)){
    typename Kernel::Vector_3 vec = compute_face_normal(f,pmesh, ppmap, k);
    put(fnpm, f, vec);
  }
}



///\cond SKIP_IN_MANUAL

template <typename PolygonMesh, typename VertexNormalPropertyMap>
void
compute_vertex_normals(const PolygonMesh& pmesh,
                      VertexNormalPropertyMap npm)
{  
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  compute_vertex_normals(pmesh, npm, get(vertex_point, pmesh), Kernel());
}


template <typename PolygonMesh, typename VertexNormalPropertyMap, typename PointPropertyMap>
void
compute_vertex_normals(const PolygonMesh& pmesh,
                      VertexNormalPropertyMap npm,
                      PointPropertyMap ppmap)
{  
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  compute_vertex_normals(pmesh, npm, ppmap, Kernel());
}

///\endcond
/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal for all vertices of the polygon mesh.
* @tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the point type of `PolygonMesh`.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam VertexNormalPropertyMap the property map in which the normals are written.
* @tparam PointPropertyMap the property map with the points associated to the vertices.
*
* @param f the face on which the normal is computed
* @param pmesh the polygon mesh
*/
template <typename PolygonMesh
          , typename VertexNormalPropertyMap
          , typename PointPropertyMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, vertex_point_t>::type
#endif
          , typename Kernel
#ifdef DOXYGEN_RUNNING
          = typename Kernel_traits<
              typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
          >
void
compute_vertex_normals(const PolygonMesh& pmesh
                      , VertexNormalPropertyMap vnpm
                      , PointPropertyMap ppmap
#ifdef DOXYGEN_RUNNING
                      = get(vertex_point, pmesh)
#endif
                      ,const Kernel& k
#ifdef DOXYGEN_RUNNING
                      = Kernel()
#endif
                      )
{
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v;
  BOOST_FOREACH(v, vertices(pmesh)){
    typename Kernel::Vector_3 vec = compute_vertex_normal(v,pmesh, ppmap, k);
    put(vnpm, v, vec);
  }
}

///\cond SKIP_IN_MANUAL

template <typename PolygonMesh, typename FaceNormalPropertyMap>
void
compute_face_normals(const PolygonMesh& pmesh,
                      FaceNormalPropertyMap npm)
{  
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  compute_face_normals(pmesh, npm, get(vertex_point, pmesh), Kernel());
}


template <typename PolygonMesh, typename FaceNormalPropertyMap, typename PointPropertyMap>
void
compute_face_normals(const PolygonMesh& pmesh,
                      FaceNormalPropertyMap npm,
                      PointPropertyMap ppmap)
{  
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  compute_face_normals(pmesh, npm, ppmap, Kernel());
}

///\endcond

/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal for all faces of the polygon mesh.
* @tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the point type of `PolygonMesh`.
* @tparam PolygonMesh a model of `FaceGraph`
* @tparam FaceNormalPropertyMap the property map in which the normals are written.
* @tparam PointPropertyMap the property map with the points associated to the vertices.
*

* @param pmesh the polygon mesh
*/
template <typename PolygonMesh
          , typename FaceNormalPropertyMap
          , typename PointPropertyMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, vertex_point_t>::type
#endif
          , typename Kernel
#ifdef DOXYGEN_RUNNING
          = typename Kernel_traits<
              typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
          >
void
compute_face_normals(const PolygonMesh& pmesh
                      , FaceNormalPropertyMap npm
                      , PointPropertyMap ppmap
#ifdef DOXYGEN_RUNNING
                      = get(vertex_point, pmesh)
#endif
                      ,const Kernel& k
#ifdef DOXYGEN_RUNNING
                      = Kernel()
#endif
                      )
{
  typename boost::graph_traits<PolygonMesh>::face_descriptor f;
  BOOST_FOREACH(f, faces(pmesh)){
    typename Kernel::Vector_3 vec = compute_face_normal(f,pmesh, ppmap, k);
    put(npm, f, vec);
  }
}


///\cond SKIP_IN_MANUAL

template <typename PolygonMesh>
typename Kernel_traits<
  typename boost::property_traits<
    typename boost::property_map<
      PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel::Vector_3
compute_face_normal(
  typename boost::graph_traits<PolygonMesh>::face_descriptor f,
  const PolygonMesh& pmesh)
{
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  return compute_face_normal(f, pmesh, get(vertex_point, pmesh), Kernel());
}
 


template <typename PolygonMesh, typename PointPropertyMap>
typename Kernel_traits<
  typename boost::property_traits<PointPropertyMap>::value_type>::Kernel::Vector_3
compute_face_normal(
  typename boost::graph_traits<PolygonMesh>::face_descriptor f,
  const PolygonMesh& pmesh,
  PointPropertyMap ppmap)
{
  typedef typename Kernel_traits<
    typename boost::property_traits<PointPropertyMap>::value_type>::Kernel Kernel;
  return compute_face_normal(f, pmesh, ppmap, Kernel());
}
 
/// \endcond 


/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal to face `f`.
* @tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the point type of `PolygonMesh`.
* @tparam PolygonMesh a model of `FaceGraph`
*
* @param f the face on which the normal is computed
* @param pmesh the polygon mesh to which `f` belongs
*/
template <typename PolygonMesh
          , typename PointPropertyMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
#endif
          , typename Kernel
#ifdef DOXYGEN_RUNNING
          = typename Kernel_traits<
              typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
          >
typename Kernel::Vector_3
compute_face_normal(typename boost::graph_traits<PolygonMesh>::face_descriptor f,
                     const PolygonMesh& pmesh
                     , PointPropertyMap ppmap 
#ifdef DOXYGEN_RUNNING
                     = get(vertex_point_t, pmesh)
#endif
                     , const Kernel& k
#ifdef DOXYGEN_RUNNING 
                     = Kernel()
#endif
                     )
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  Vector normal = CGAL::NULL_VECTOR;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor he = halfedge(f, pmesh);
  halfedge_descriptor end = he;
  do
  {
    const Point& prv = ppmap[target(prev(he, pmesh), pmesh)];
    const Point& curr = ppmap[target(he, pmesh)];
    const Point& nxt = ppmap[target(next(he, pmesh), pmesh)];
    Vector n = CGAL::cross_product(nxt - curr, prv - curr);
    normal = normal + n;

    he = next(he, pmesh);
  } while (he != end);

  return normal / std::sqrt(normal * normal);
}

///\cond SKIP_IN_MANUAL  

template <typename PolygonMesh>
typename Kernel_traits<
  typename boost::property_traits<
    typename boost::property_map<
      PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel::Vector_3
compute_vertex_normal(
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
  const PolygonMesh& pmesh)
{
  typedef typename Kernel_traits<
    typename boost::property_traits<
      typename boost::property_map<
        PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  return compute_vertex_normal(v, pmesh, get(vertex_point, pmesh), Kernel());
}
 


template <typename PolygonMesh, typename PointPropertyMap>
typename Kernel_traits<
  typename boost::property_traits<
    typename boost::property_map<
      PolygonMesh, CGAL::vertex_point_t>::type>::value_type>::Kernel::Vector_3
compute_vertex_normal(
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
  const PolygonMesh& pmesh,
  PointPropertyMap ppmap)
{ typedef typename Kernel_traits<
    typename boost::property_traits<PointPropertyMap>::value_type>::Kernel Kernel;
  return compute_vertex_normal(v, pmesh, ppmap, Kernel);
}
  
/// \endcond 



/**
* \ingroup PkgPolygonMeshProcessing
* computes the unit normal at vertex `v` as the average of the normals of incident faces.
* @tparam Kernel a %CGAL `Kernel` with `FT` a model of `FieldWithSqrt`
* @tparam PolygonMesh a model of `FaceGraph`
*
* @param v the vertex around which the normal is computed
* @param pmesh the polygon mesh to which `v` belongs
*/
template<typename PolygonMesh
          , typename PointPropertyMap
#ifdef DOXYGEN_RUNNING
          = typename boost::property_map<PolygonMesh, vertex_point_t>::type
#endif
         , typename Kernel
#ifdef DOXYGEN_RUNNING
         = typename Kernel_traits<
             typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
         >
typename Kernel::Vector_3
compute_vertex_normal(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                      const PolygonMesh& pmesh,
                      PointPropertyMap ppmap 
#ifdef DOXYGEN_RUNNING
                      = get(vertex_point_t, pmesh)
#endif
                      , const Kernel& k
#ifdef DOXYGEN_RUNNING 
                      = Kernel()
#endif
                      )
{
  typedef typename Kernel::Vector_3 Vector;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  Vector normal = CGAL::NULL_VECTOR;
  halfedge_descriptor he = halfedge(v, pmesh);
  halfedge_descriptor end = he;
  do
    {
    if (!is_border(he, pmesh))
    {
      Vector n = compute_face_normal(face(he, pmesh), pmesh, ppmap, k);
      normal = normal + (n / std::sqrt(n*n));
    }
    he = opposite(next(he, pmesh), pmesh);
  } while (he != end);

  return normal / std::sqrt(normal * normal);
}

} 

} // end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
