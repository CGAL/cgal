// This file has been adopted from CGAL, to integrate the below mentioned paper
//
// Paper         : Learning to Reconstruct Symmetric Shapes using Planar Parameterization of 3D Surface
// Author(s)     : Hardik Jain, Manuel Wöllhaf, Olaf Hellwich
// Conference    : IEEE International Conference on Computer Vision Workshops (ICCVW) 2019
//

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_PARAMETERIZE_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_PARAMETERIZE_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/Bool_property_map.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

/// \file Iterative_parameterize.h
namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMainFunction
///
/// Compute a one-to-one mapping from a 3D triangle surface `mesh` to a
/// simple 2D domain.
/// The mapping is piecewise linear on the triangle mesh.
/// The result is a pair `(u,v)` of parameter coordinates for each vertex of the input mesh.
///
/// A one-to-one mapping may be guaranteed or not, depending on
/// the chosen Parameterizer algorithm.
///
/// \tparam TriangleMesh must be a model of `FaceGraph`.
/// \tparam Parameterizer must be a model of `Parameterizer_3`.
/// \tparam HD must be the halfedge_descriptor type corresponding to the graph
///         traits of TriangleMesh.
/// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
///         as value type.
///
/// \param mesh a triangulated surface.
/// \param parameterizer a parameterizer.
/// \param bhd a halfedge descriptor on the boundary of `mesh`.
/// \param uvm an instanciation of the class `VertexUVmap`.
/// \param iterations an integer number of iterations to run the parameterization.
/// \param error return error value of the iterative process.
///
/// \pre `mesh` must be a triangular mesh.
/// \pre The mesh border must be mapped onto a convex polygon
///   (for fixed border parameterizations).
/// \pre The vertices must be indexed (vimap must be initialized).
///
template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap>
Error_code parameterize(TriangleMesh& mesh,
    Parameterizer parameterizer,
    HD bhd,
    VertexUVmap uvm,
    int& iterations,
    double& error)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  CGAL::Polygon_mesh_processing::connected_component(
      face(opposite(bhd, mesh), mesh),
      mesh,
      boost::make_function_output_iterator(
          internal::Index_map_filler<TriangleMesh, Indices>(mesh, indices)));
  boost::associative_property_map<Indices> vipm(indices);

  boost::unordered_set<vertex_descriptor> vs;
  internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);

  return parameterizer.parameterize(mesh, bhd, uvm, vipm, vpm, iterations, error);
}

template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap>
Error_code parameterize(TriangleMesh& mesh,
    Parameterizer parameterizer,
    HD bhd,
    VertexUVmap uvm,
    int& iterations)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  CGAL::Polygon_mesh_processing::connected_component(
      face(opposite(bhd, mesh), mesh),
      mesh,
      boost::make_function_output_iterator(
          internal::Index_map_filler<TriangleMesh, Indices>(mesh, indices)));
  boost::associative_property_map<Indices> vipm(indices);

  boost::unordered_set<vertex_descriptor> vs;
  internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);
  double error;
  return parameterizer.parameterize(mesh, bhd, uvm, vipm, vpm, iterations, error);
}

template <class TriangleMesh, class Parameterizer, class HD, class VertexUVmap>
Error_code parameterize(TriangleMesh& mesh,
    Parameterizer parameterizer,
    HD bhd,
    VertexUVmap uvm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  CGAL::Polygon_mesh_processing::connected_component(
      face(opposite(bhd, mesh), mesh),
      mesh,
      boost::make_function_output_iterator(
          internal::Index_map_filler<TriangleMesh, Indices>(mesh, indices)));
  boost::associative_property_map<Indices> vipm(indices);

  boost::unordered_set<vertex_descriptor> vs;
  internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpm(vs);
  int iterations = 0;
  double error;
  return parameterizer.parameterize(mesh, bhd, uvm, vipm, vpm, iterations, error);
}

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ITERATIVE_PARAMETERIZE_H
