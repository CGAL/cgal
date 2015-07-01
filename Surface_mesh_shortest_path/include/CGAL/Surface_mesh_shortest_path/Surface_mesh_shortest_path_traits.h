// Copyright (c) 2014 GeometryFactory
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
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_TRAITS_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_TRAITS_H

#include <CGAL/Surface_mesh_shortest_path/barycentric.h>
#include <CGAL/Surface_mesh_shortest_path/function_objects.h>

#include <ostream>
#include <cstddef>

#include <boost/array.hpp>

namespace CGAL {

/*!
\ingroup PkgSurfaceMeshShortestPathTraitsClasses

\brief A model of the concept `SurfaceMeshShortestPathTraits`
as required by the `Surface_mesh_shortest_path` class.

\tparam K A \cgal Kernel

\tparam G A model of `FaceListGraph`

\cgalModels `SurfaceMeshShortestPathTraits`
*/
template <
  class K,
  class G>
class Surface_mesh_shortest_path_traits : public K
{
public:

  /// Kernel type
  typedef K Kernel;

  /// Triangle mesh type
  typedef G Triangle_mesh;

  typedef typename Kernel::FT FT;

  /// Barycentric coordinate type
  typedef typename CGAL::cpp11::array<FT,3> Barycentric_coordinate;

  // Predicates
public:
  typedef typename Surface_mesh_shortest_paths_3::Compare_relative_intersection_along_segment_2<Kernel> Compare_relative_intersection_along_segment_2;
  typedef typename Surface_mesh_shortest_paths_3::Is_saddle_vertex<Kernel, Triangle_mesh> Is_saddle_vertex;

  // Constructions
public:
  class Construct_barycentric_coordinate
  {
  public:
    typedef Barycentric_coordinate result_type;

    result_type operator() (const FT& a, const FT& b, const FT& c) const
    {
      Barycentric_coordinate output;
      output[0] = a;
      output[1] = b;
      output[2] = c;
      return output;
    }
  };

  class Construct_barycentric_coordinate_weight
  {
  public:
    typedef FT result_type;

    result_type operator() (const Barycentric_coordinate b, std::size_t i) const
    {
      return b[i % 3];
    }
  };

  typedef typename Surface_mesh_shortest_paths_3::Construct_triangle_3_to_triangle_2_projection<K> Construct_triangle_3_to_triangle_2_projection;
  typedef typename Surface_mesh_shortest_paths_3::Construct_triangle_3_along_segment_2_flattening<K> Construct_triangle_3_along_segment_2_flattening;
  typedef typename Surface_mesh_shortest_paths_3::Compute_parametric_distance_along_segment_2<K> Compute_parametric_distance_along_segment_2;
  typedef typename Surface_mesh_shortest_paths_3::Construct_barycentric_coordinate_in_triangle_2<K, Barycentric_coordinate, Construct_barycentric_coordinate> Construct_barycentric_coordinate_in_triangle_2;
  typedef typename Surface_mesh_shortest_paths_3::Construct_barycentric_coordinate_in_triangle_3<K, Barycentric_coordinate, Construct_barycentric_coordinate> Construct_barycentric_coordinate_in_triangle_3;
  typedef typename Surface_mesh_shortest_paths_3::Classify_barycentric_coordinate<Barycentric_coordinate, Construct_barycentric_coordinate_weight> Classify_barycentric_coordinate;

private:
  Kernel m_kernel;
  Construct_barycentric_coordinate m_construct_barycentric_coordinate_object;
  Construct_barycentric_coordinate_weight m_construct_barycentric_coordinate_weight_object;
  Classify_barycentric_coordinate m_classify_barycentric_coordinate_object;
  Construct_triangle_3_to_triangle_2_projection m_construct_triangle_3_to_triangle_2_projection_object;
  Construct_triangle_3_along_segment_2_flattening m_construct_triangle_3_along_segment_2_flattening_object;
  Construct_barycentric_coordinate_in_triangle_2 m_construct_barycentric_coordinate_in_triangle_2_object;
  Construct_barycentric_coordinate_in_triangle_3 m_construct_barycentric_coordinate_in_triangle_3_object;
  Compare_relative_intersection_along_segment_2 m_compare_relative_intersection_along_segment_2_object;
  Is_saddle_vertex m_is_saddle_vertex_object;
  Compute_parametric_distance_along_segment_2 m_compute_parametric_distance_along_segment_2_object;

public:

  Surface_mesh_shortest_path_traits()
  {
  }

  Surface_mesh_shortest_path_traits(const Kernel& kernel)
    : m_kernel(kernel)
    , m_classify_barycentric_coordinate_object(m_construct_barycentric_coordinate_weight_object)
    , m_construct_triangle_3_to_triangle_2_projection_object(m_kernel)
    , m_construct_triangle_3_along_segment_2_flattening_object(m_kernel)
    , m_is_saddle_vertex_object(m_kernel, m_construct_triangle_3_to_triangle_2_projection_object, m_construct_triangle_3_along_segment_2_flattening_object)
    , m_compare_relative_intersection_along_segment_2_object(m_kernel)
    , m_construct_barycentric_coordinate_in_triangle_2_object(m_construct_barycentric_coordinate_object, m_kernel.construct_vector_2_object(), m_kernel.compute_scalar_product_2_object())
    , m_construct_barycentric_coordinate_in_triangle_3_object(m_construct_barycentric_coordinate_object, m_kernel.construct_vector_3_object(), m_kernel.compute_scalar_product_3_object())
  {
  }

  Construct_barycentric_coordinate construct_barycentric_coordinate_object() const { return m_construct_barycentric_coordinate_object; }
  Construct_barycentric_coordinate_weight construct_barycentric_coordinate_weight_object() const { return m_construct_barycentric_coordinate_weight_object; }
  Classify_barycentric_coordinate classify_barycentric_coordinate_object() const { return m_classify_barycentric_coordinate_object; }

  Is_saddle_vertex is_saddle_vertex_object() const { return m_is_saddle_vertex_object; }

  Compare_relative_intersection_along_segment_2 compare_relative_intersection_along_segment_2_object() const { return m_compare_relative_intersection_along_segment_2_object; }
  Construct_triangle_3_to_triangle_2_projection construct_triangle_3_to_triangle_2_projection_object() const { return m_construct_triangle_3_to_triangle_2_projection_object; }
  Construct_triangle_3_along_segment_2_flattening construct_triangle_3_along_segment_2_flattening_object() const { return m_construct_triangle_3_along_segment_2_flattening_object; }
  Construct_barycentric_coordinate_in_triangle_2 construct_barycentric_coordinate_in_triangle_2_object() const { return m_construct_barycentric_coordinate_in_triangle_2_object; }
  Construct_barycentric_coordinate_in_triangle_3 construct_barycentric_coordinate_in_triangle_3_object() const { return m_construct_barycentric_coordinate_in_triangle_3_object; }
  Compute_parametric_distance_along_segment_2 compute_parametric_distance_along_segment_2_object() const { return m_compute_parametric_distance_along_segment_2_object; }
};

template <class Kernel, class Triangle_mesh>
std::ostream& operator<<(std::ostream& os, typename Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh>::Barycentric_coordinate b)
{
  return os << b[0] << " " << b[1] << " " << b[2];
}

#ifndef DOXYGEN_RUNNING // needed due to a bug in doxygen
/*!
\ingroup PkgSurfaceMeshShortestPathTraitsClasses

\internal

\brief Provides an implementation of the SurfaceMeshShortestPathTraits
model which uses an exact Kernel during the unfolding operations to achieve better overall precision

\tparam K Kernel Type

\tparam G triangle mesh type

\cgalModels `SurfaceMeshShortestPathTraits`
*/
template <
  class K,
  class G>
class Surface_mesh_shortest_path_traits_with_robust_unfolding : public Surface_mesh_shortest_path_traits<K,G>
{
public:
  typedef K Kernel;
  typedef typename Surface_mesh_shortest_paths_3::Robust_project_triangle_3_to_triangle_2<K> Construct_triangle_3_to_triangle_2_projection;
  typedef typename Surface_mesh_shortest_paths_3::Robust_flatten_triangle_3_along_segment_2<K> Construct_triangle_3_along_segment_2_flattening;

private:

  Construct_triangle_3_to_triangle_2_projection m_robust_construct_triangle_3_to_triangle_2_projection_object;
  Construct_triangle_3_along_segment_2_flattening m_robust_flatten_triangle_3_along_segment_2;

public:

  Surface_mesh_shortest_path_traits_with_robust_unfolding()
  {
  }

  Surface_mesh_shortest_path_traits_with_robust_unfolding(const Kernel& kernel)
    : Surface_mesh_shortest_path_traits<K,G>(kernel)
    , m_robust_construct_triangle_3_to_triangle_2_projection_object(kernel)
    , m_robust_flatten_triangle_3_along_segment_2(kernel)
  {
  }

  Construct_triangle_3_to_triangle_2_projection construct_triangle_3_to_triangle_2_projection_object() const { return m_robust_construct_triangle_3_to_triangle_2_projection_object; }
  Construct_triangle_3_along_segment_2_flattening construct_triangle_3_along_segment_2_flattening_object() const { return m_robust_flatten_triangle_3_along_segment_2; }
};
#endif

} // namespace CGAL


#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_TRAITS_H
