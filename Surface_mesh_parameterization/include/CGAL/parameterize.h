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
template <class ParameterizationMesh_3>
typename Parameterizer_traits_3<ParameterizationMesh_3>::Error_code
parameterize(ParameterizationMesh_3& mesh)  ///< 3D mesh, model of ParameterizationMesh_3 concept
{
    Mean_value_coordinates_parameterizer_3<ParameterizationMesh_3> parameterizer;
    return parameterizer.parameterize(mesh);
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
///
template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
typename Parameterizer_traits_3<TriangleMesh>::Error_code
parameterize(TriangleMesh& mesh,
             Parameterizer parameterizer,
             HD bhd,
             VertexUVmap uvm,
             VertexIndexMap vimap,
             VertexParameterizedMap vpm)
{
  return parameterizer.parameterize(mesh, bhd, uvm, vimap, vpm);
}


template <class TM>
class Seam_mesh;


template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
typename Parameterizer_traits_3<Seam_mesh<TriangleMesh> >::Error_code
parameterize(Seam_mesh<TriangleMesh>& mesh,
             Parameterizer parameterizer,
             HD bhd,
             VertexUVmap uvm,
             VertexIndexMap vimap,
             VertexParameterizedMap vpm)
{
  std::cerr << "treat Seam_mesh"<< std::endl;
  Seam_mesh_uv_map<TriangleMesh,VertexUVmap>  putter(mesh,uvm);
  return parameterizer.parameterize(mesh, bhd, putter, vimap, vpm);
}


} //namespace CGAL

#endif //CGAL_PARAMETERIZE_H
