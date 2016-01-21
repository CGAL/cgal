// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_PARAMETERIZE_H
#define CGAL_PARAMETERIZE_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>

/// \file parameterize.h

namespace CGAL {


/// \ingroup  PkgSurfaceParameterizationMainFunction
///
/// Compute a one-to-one mapping from a 3D triangle surface `mesh` to a
/// 2D circle, using Floater Mean Value Coordinates algorithm.
/// A one-to-one mapping is guaranteed.
///
/// The mapping is piecewise linear on the input mesh triangles.
/// The result is a (u,v) pair of parameter coordinates for each vertex of the input mesh.
///
/// \pre `mesh` must be a surface with one connected component.
/// \pre `mesh` must be a triangular mesh.
///
template <class TriangleMesh>
typename Parameterizer_traits_3<TriangleMesh>::Error_code
parameterize(TriangleMesh& mesh)  ///< 3D mesh, model of TriangleMesh concept
{
    Mean_value_coordinates_parameterizer_3<TriangleMesh> parameterizer;
    return parameterizer.parameterize(mesh);
}


  namespace Parameterization {
    
    template <typename Mesh, typename Map>
    struct Vertices {

      Vertices(const Mesh& mesh, Map& map)
        : mesh(mesh), map(&map), index(0)
      {}
    
      void operator()(const typename boost::graph_traits<Mesh>::face_descriptor& fd)
      {
        BOOST_FOREACH(typename boost::graph_traits<Mesh>::vertex_descriptor vd, vertices_around_face(halfedge(fd,mesh),mesh)){
          if(map->find(vd) == map->end()){
            (*map)[vd] = index++;
          }
        }
      }

      const Mesh& mesh;
      mutable Map* map;
      int index;
    };
  }


/// \ingroup  PkgSurfaceParameterizationMainFunction
///
/// Compute a one-to-one mapping from a 3D triangle surface `mesh` to a
/// simple 2D domain.
/// The mapping is piecewise linear on the triangle mesh.
/// The result is a pair (u,v) of parameter coordinates for each vertex of the input mesh.
///
/// One-to-one mapping may be guaranteed or
/// not, depending on the chosen Parameterizer algorithm.
///
/// \pre `mesh` must be a surface with one connected component.
/// \pre `mesh` must be a triangular mesh.
/// \pre The mesh border must be mapped onto a convex polygon
///   (for fixed border parameterizations).

template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
typename Parameterizer_traits_3<TriangleMesh>::Error_code
parameterize(TriangleMesh& mesh,
             Parameterizer parameterizer,
             HD bhd,
             VertexUVmap uvm,
             VertexIndexMap vimap,
             VertexParameterizedMap vpm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef std::map<vertex_descriptor,int> Indices;
  Indices indices;
  CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd,mesh),mesh),
                                                     mesh,
                                                     boost::make_function_output_iterator(Parameterization::Vertices<TriangleMesh,Indices>(mesh,indices)));
  return parameterizer.parameterize(mesh, bhd, uvm, boost::make_assoc_property_map(indices), vpm);
}

  
template <class TriangleMesh, class HD, class VertexUVmap, typename VertexParameterizedMap>
typename Parameterizer_traits_3<TriangleMesh>::Error_code
parameterize(TriangleMesh& mesh,
             HD bhd,
             VertexUVmap uvm,
             VertexParameterizedMap vpm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef std::map<vertex_descriptor,int> Indices;
  Indices indices;
  BOOST_FOREACH(vertex_descriptor v, vertices(mesh)){
    indices[v] = -1;
  }
  CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd,mesh),mesh),
                                                     mesh,
                                                     boost::make_function_output_iterator(Parameterization::Vertices<TriangleMesh,Indices>(mesh,indices)));
  Mean_value_coordinates_parameterizer_3<TriangleMesh> parameterizer;
  return parameterizer.parameterize(mesh, bhd, uvm, boost::make_assoc_property_map(indices), vpm);
}

template <class TM>
class Seam_mesh;


template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap, typename VertexParameterizedMap>
typename Parameterizer_traits_3<Seam_mesh<TriangleMesh> >::Error_code
parameterize(Seam_mesh<TriangleMesh>& mesh,
             Parameterizer parameterizer,
             HD bhd,
             VertexUVmap uvm,
             VertexParameterizedMap vpm)
{ 
  typedef typename boost::graph_traits<Seam_mesh<TriangleMesh> >::vertex_descriptor vertex_descriptor;
  typedef std::map<vertex_descriptor,int> Vertex_index_map;
  Vertex_index_map vim;
  boost::associative_property_map<Vertex_index_map> vipm(vim);
  mesh.initialize_vertex_index_map(bhd,vipm);
  Seam_mesh_uv_map<TriangleMesh,VertexUVmap>  putter(mesh,uvm);
  return parameterizer.parameterize(mesh, bhd, putter, vipm, vpm);
}


} //namespace CGAL

#endif //CGAL_PARAMETERIZE_H
