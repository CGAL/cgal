// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_CONVERT_IMAGE_TO_SURFACE_MESH_H
#define CGAL_FEATURE_GRAPH_INTERNAL_CONVERT_IMAGE_TO_SURFACE_MESH_H

#include <CGAL/license/Feature_graph.h>

#include <CGAL/Surface_mesh.h>

namespace CGAL
{

namespace Feature_graph
{

namespace internal
{

template <typename Image>
float image_value_with_check(const Image& image, const int& i, const int&  j, const int&  k, const float& value_outside = 0)
{
  if ( i >= 0 && static_cast<std::size_t>(i) < image.xdim() &&
        j >= 0 && static_cast<std::size_t>(j) < image.ydim() &&
        k >= 0 && static_cast<std::size_t>(k) < image.zdim() )
    return image.value(i, j, k);
  else
    return value_outside;
}

template <typename Point_3>
struct Get_points_property_from_surface_mesh
{
  typedef Surface_mesh<Point_3> Mesh;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typename Mesh::Property_map<vertex_descriptor, Point_3>& operator()(const std::size_t mesh_index, Mesh& mesh) const
  {
    return mesh.points();
  }
};

struct Get_surface_mesh_domains_and_orientations
{
  template <typename DomainOutputIterator, typename OrientationOutputIterator>
  void operator()(const int& axis, const int& i, const int& j, const int& k, const float& base_value, const float& neighbor_value,
    DomainOutputIterator domain_out, OrientationOutputIterator orientation_out) const
  {
    const bool orientation = (base_value < neighbor_value);
    std::array<float,2> surface_domain = {base_value, neighbor_value};
    if (!orientation)
      std::swap(surface_domain[0], surface_domain[1]);

    *domain_out++ = surface_domain;
    *orientation_out++ = orientation;
  }
};

struct Get_volume_mesh_domains_and_orientations
{
  template <typename DomainOutputIterator, typename OrientationOutputIterator>
  void operator()(const int& axis, const int& i, const int& j, const int& k, const float& base_value, const float& neighbor_value,
    DomainOutputIterator domain_out, OrientationOutputIterator orientation_out) const
  {
    *domain_out++ = base_value;
    *orientation_out++ = false;
    *domain_out++ = neighbor_value;
    *orientation_out++ = true;
  }
};

// constructs one or more surface mesh from the labelled image.
//   The parameters enable the user to create a mesh per surface, a mesh per domain, or a more complex scenario
// PolygonMesh is model of MutableFaceGraph
// GetPointsProperty returns a property (by adress) that holds the points of a mesh
// PointsPropertyBuilder constructs the PointsPropertyMap from a mesh
// PointelToMeshVertexMap will hold the mesh vertex_descriptor from a pointel. Note: this is a map, it should define `find(...)`, `end()`, `emplace(...)`.
// MeshDomain is an identifier representing the domain of the mesh. it is constructed from a surfel (see GetMeshDomainsAndOrientationsFromSurfel)
// GetMeshDomainsAndOrientationsFromSurfel constructs one or more mesh domains from a surfel. The surfel will be added to all the meshes with the same domain.
template <typename K, typename Image, typename PolygonMesh,
      typename GetPointsProperty, typename PointelToMeshVertexMap, typename MeshDomain, typename GetMeshDomainsAndOrientationsFromSurfel>
void convert_labelled_image_to_face_graphs(const Image& image,
  std::vector<PolygonMesh>& meshes,
  std::vector<MeshDomain>& mesh_domains,
  std::vector<PointelToMeshVertexMap>& image_point_to_mesh_vertex,
  GetPointsProperty get_points_property = GetPointsProperty(),
  GetMeshDomainsAndOrientationsFromSurfel get_mesh_domains_and_orientation_from_surfel = GetMeshDomainsAndOrientationsFromSurfel(),
  const float& value_outside = 0)
{
  typedef typename K::Point_3 Point_3;
  typedef PolygonMesh Mesh;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;

  std::map<MeshDomain, std::size_t> mesh_domain_to_mesh_index;

  int i,j,k;

  const double& vx = image.vx();
  const double& vy = image.vy();
  const double& vz = image.vz();
  const double& tx = image.tx();
  const double& ty = image.ty();
  const double& tz = image.tz();

  const int xdim = static_cast<int>(image.xdim());
  const int ydim = static_cast<int>(image.ydim());
  const int zdim = static_cast<int>(image.zdim());
  // on each dimension, add one to the image point index to represent the border
  const int xdim_p = xdim+1;
  const int ydim_p = ydim+1;
  const int xydim_p = ydim_p*xdim_p;
  std::vector<MeshDomain> current_surfel_mesh_domains;
  std::vector<bool> current_surfel_mesh_orientations;
  std::vector<std::size_t> current_surfel_mesh_indices;
  std::vector<std::array<vertex_descriptor, 4>> current_surfel_vertices;

  auto add_surfel_to_mesh = [&](const int& axis, const int& i, const int& j, const int& k, const float& base_value, const float& neighbor_value)
  {
    // get or create mesh
    current_surfel_mesh_domains.clear();
    current_surfel_mesh_orientations.clear();
    current_surfel_mesh_indices.clear();
    get_mesh_domains_and_orientation_from_surfel(axis, i, j, k, base_value, neighbor_value, std::back_inserter(current_surfel_mesh_domains), std::back_inserter(current_surfel_mesh_orientations));
    for (const MeshDomain& mesh_domain : current_surfel_mesh_domains)
    {
      const auto& [mesh_domain_and_mesh_iterator, inserted_surface] = mesh_domain_to_mesh_index.try_emplace(mesh_domain, meshes.size());
      if (inserted_surface)
      {
        const std::size_t new_size = meshes.size()+1;
        meshes.resize(new_size);
        mesh_domains.push_back(mesh_domain);
        image_point_to_mesh_vertex.resize(new_size);
      }
      const std::size_t& mesh_index = mesh_domain_and_mesh_iterator->second;
      current_surfel_mesh_indices.push_back(mesh_index);
    }
    // get or create surfel vertices
    std::array<int, 3> base_image_point{0,0,0};
    std::array<int, 3> u_axis{0,0,0};
    std::array<int, 3> v_axis{0,0,0};
    base_image_point[axis] = 1;
    u_axis[(axis+1)%3] = 1;
    v_axis[(axis+2)%3] = 1;
    current_surfel_vertices.clear();
    current_surfel_vertices.resize(current_surfel_mesh_indices.size());
    const std::size_t i_base = i + base_image_point[0];
    const std::size_t j_base = j + base_image_point[1];
    const std::size_t k_base = k + base_image_point[2];
    for (int v = 0; v < 2; v++)
      for (int u = 0; u < 2; u++)
      {
        const std::size_t surfel_vertices_index = v*2+u;
        const std::size_t i_vertex = i_base + u * u_axis[0] + v * v_axis[0];
        const std::size_t j_vertex = j_base + u * u_axis[1] + v * v_axis[1];
        const std::size_t k_vertex = k_base + u * u_axis[2] + v * v_axis[2];
        const std::size_t vertex_index = k_vertex*xydim_p + j_vertex *xdim_p + i_vertex;
        const Point_3 world_pointel = Point_3(i_vertex * vx + tx, j_vertex * vy + ty, k_vertex * vz + tz);
        for (std::size_t current_surfel_mesh_index = 0, number_of_meshes = current_surfel_mesh_indices.size(); current_surfel_mesh_index < number_of_meshes; current_surfel_mesh_index++)
        {
          const std::size_t& mesh_index = current_surfel_mesh_indices[current_surfel_mesh_index];
          Mesh& mesh = meshes[mesh_index];
          PointelToMeshVertexMap& image_point_to_current_mesh_vertex = image_point_to_mesh_vertex[mesh_index];
          const auto& point_and_vertex_iterator = image_point_to_current_mesh_vertex.find(vertex_index);
          vertex_descriptor surfel_vertex;
          if (point_and_vertex_iterator == image_point_to_current_mesh_vertex.end())
          {
            surfel_vertex = add_vertex(mesh);
            put(get_points_property(mesh_index, mesh), surfel_vertex, world_pointel);
            image_point_to_current_mesh_vertex.emplace(vertex_index, surfel_vertex);
          }
          else
            surfel_vertex = point_and_vertex_iterator->second;
          current_surfel_vertices[current_surfel_mesh_index][surfel_vertices_index] = surfel_vertex;
        }
      }
    // add surfel to each mesh
    for (std::size_t current_surfel_mesh_index = 0, number_of_meshes = current_surfel_mesh_indices.size(); current_surfel_mesh_index < number_of_meshes; current_surfel_mesh_index++)
    {
      const std::size_t& mesh_index = current_surfel_mesh_indices[current_surfel_mesh_index];
      const bool& orientation = current_surfel_mesh_orientations[current_surfel_mesh_index];
      std::array<vertex_descriptor, 4>& surfel_vertices = current_surfel_vertices[current_surfel_mesh_index];
      if (orientation)
        std::swap(surfel_vertices[0], surfel_vertices[1]); // 0 1 2 3 --> 1 0 2 3
      else
        std::swap(surfel_vertices[2], surfel_vertices[3]); // 0 1 2 3 --> 0 1 3 2
      face_descriptor face = CGAL::Euler::add_face(surfel_vertices, meshes[mesh_index]);
    }
  };

  for ( k = -1 ; k < zdim ; k+=1 )
  {
    for ( j = -1 ; j < ydim ; j+=1 )
    {
      for ( i = -1 ; i < xdim ; i+=1 )
      {
        const float & value = image_value_with_check(image, i, j, k, value_outside);

        const float & value_x = image_value_with_check(image, i+1, j, k, value_outside);
        const float & value_y = image_value_with_check(image, i, j+1, k, value_outside);
        const float & value_z = image_value_with_check(image, i, j, k+1, value_outside);

        // x direction
        if (value != value_x)
          add_surfel_to_mesh(0, i, j, k, value, value_x);
        // y direction
        if (value != value_y)
          add_surfel_to_mesh(1, i, j, k, value, value_y);
        // z direction
        if (value != value_z)
          add_surfel_to_mesh(2, i, j, k, value, value_z);
      }
    }
  }
}

} /* namespace internal */

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_INTERNAL_CONVERT_IMAGE_TO_SURFACE_MESH_H