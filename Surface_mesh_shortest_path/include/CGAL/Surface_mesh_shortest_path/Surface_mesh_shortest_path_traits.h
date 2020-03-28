// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_TRAITS_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_TRAITS_H

#include <CGAL/license/Surface_mesh_shortest_path.h>

#include <CGAL/disable_warnings.h>

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

\tparam TriangleMesh A model of `FaceListGraph`

\cgalModels `SurfaceMeshShortestPathTraits`
*/
template <
  class K,
  class TriangleMesh>
class Surface_mesh_shortest_path_traits : public K
{
public:

  /// Kernel type
  typedef K Kernel;

  /// Triangle mesh type
  typedef TriangleMesh Triangle_mesh;

  typedef typename Kernel::FT FT;

  /// Barycentric coordinates type
  typedef typename std::array<FT,3> Barycentric_coordinates;

  // Predicates
  typedef typename Surface_mesh_shortest_paths_3::Compare_relative_intersection_along_segment_2<Kernel> Compare_relative_intersection_along_segment_2;
  typedef typename Surface_mesh_shortest_paths_3::Is_saddle_vertex<Kernel, Triangle_mesh> Is_saddle_vertex;

  // Constructions
public:
  class Construct_barycentric_coordinates
  {
  public:
    typedef Barycentric_coordinates result_type;

    result_type operator() (const FT& a, const FT& b, const FT& c) const
    {
      Barycentric_coordinates output;
      output[0] = a;
      output[1] = b;
      output[2] = c;
      return output;
    }
  };

  class Construct_barycentric_coordinates_weight
  {
  public:
    typedef FT result_type;

    result_type operator() (const Barycentric_coordinates b, std::size_t i) const
    {
      return b[i % 3];
    }
  };

  typedef typename Surface_mesh_shortest_paths_3::Construct_triangle_3_to_triangle_2_projection<K> Construct_triangle_3_to_triangle_2_projection;
  typedef typename Surface_mesh_shortest_paths_3::Construct_triangle_3_along_segment_2_flattening<K> Construct_triangle_3_along_segment_2_flattening;
  typedef typename Surface_mesh_shortest_paths_3::Compute_parametric_distance_along_segment_2<K> Compute_parametric_distance_along_segment_2;
  typedef typename Surface_mesh_shortest_paths_3::Construct_barycentric_coordinates_in_triangle_2<K, Barycentric_coordinates, Construct_barycentric_coordinates> Construct_barycentric_coordinates_in_triangle_2;
  typedef typename Surface_mesh_shortest_paths_3::Construct_barycentric_coordinates_in_triangle_3<K, Barycentric_coordinates, Construct_barycentric_coordinates> Construct_barycentric_coordinates_in_triangle_3;
  typedef typename Surface_mesh_shortest_paths_3::Classify_barycentric_coordinates<Barycentric_coordinates, Construct_barycentric_coordinates_weight> Classify_barycentric_coordinates;

private:
  Kernel m_kernel;
  Construct_barycentric_coordinates m_construct_barycentric_coordinates_object;
  Construct_barycentric_coordinates_weight m_construct_barycentric_coordinates_weight_object;
  Classify_barycentric_coordinates m_classify_barycentric_coordinates_object;
  Construct_triangle_3_to_triangle_2_projection m_construct_triangle_3_to_triangle_2_projection_object;
  Construct_triangle_3_along_segment_2_flattening m_construct_triangle_3_along_segment_2_flattening_object;
  Construct_barycentric_coordinates_in_triangle_2 m_construct_barycentric_coordinates_in_triangle_2_object;
  Construct_barycentric_coordinates_in_triangle_3 m_construct_barycentric_coordinates_in_triangle_3_object;
  Compare_relative_intersection_along_segment_2 m_compare_relative_intersection_along_segment_2_object;
  Is_saddle_vertex m_is_saddle_vertex_object;
  Compute_parametric_distance_along_segment_2 m_compute_parametric_distance_along_segment_2_object;

public:

  Surface_mesh_shortest_path_traits()
  {
  }

  Surface_mesh_shortest_path_traits(const Kernel& kernel)
    : m_kernel(kernel)
    , m_classify_barycentric_coordinates_object(m_construct_barycentric_coordinates_weight_object)
    , m_construct_triangle_3_to_triangle_2_projection_object(m_kernel)
    , m_construct_triangle_3_along_segment_2_flattening_object(m_kernel)
    , m_is_saddle_vertex_object(m_kernel, m_construct_triangle_3_to_triangle_2_projection_object, m_construct_triangle_3_along_segment_2_flattening_object)
    , m_compare_relative_intersection_along_segment_2_object(m_kernel)
    , m_construct_barycentric_coordinates_in_triangle_2_object(m_construct_barycentric_coordinates_object, m_kernel.construct_vector_2_object(), m_kernel.compute_scalar_product_2_object())
    , m_construct_barycentric_coordinates_in_triangle_3_object(m_construct_barycentric_coordinates_object, m_kernel.construct_vector_3_object(), m_kernel.compute_scalar_product_3_object())
  {
  }

  Construct_barycentric_coordinates construct_barycentric_coordinates_object() const { return m_construct_barycentric_coordinates_object; }
  Construct_barycentric_coordinates_weight construct_barycentric_coordinates_weight_object() const { return m_construct_barycentric_coordinates_weight_object; }
  Classify_barycentric_coordinates classify_barycentric_coordinates_object() const { return m_classify_barycentric_coordinates_object; }

  Is_saddle_vertex is_saddle_vertex_object() const { return m_is_saddle_vertex_object; }

  Compare_relative_intersection_along_segment_2 compare_relative_intersection_along_segment_2_object() const { return m_compare_relative_intersection_along_segment_2_object; }
  Construct_triangle_3_to_triangle_2_projection construct_triangle_3_to_triangle_2_projection_object() const { return m_construct_triangle_3_to_triangle_2_projection_object; }
  Construct_triangle_3_along_segment_2_flattening construct_triangle_3_along_segment_2_flattening_object() const { return m_construct_triangle_3_along_segment_2_flattening_object; }
  Construct_barycentric_coordinates_in_triangle_2 construct_barycentric_coordinates_in_triangle_2_object() const { return m_construct_barycentric_coordinates_in_triangle_2_object; }
  Construct_barycentric_coordinates_in_triangle_3 construct_barycentric_coordinates_in_triangle_3_object() const { return m_construct_barycentric_coordinates_in_triangle_3_object; }
  Compute_parametric_distance_along_segment_2 compute_parametric_distance_along_segment_2_object() const { return m_compute_parametric_distance_along_segment_2_object; }

  #ifndef CGAL_NO_DEPRECATED_CODE
  //deprecated in CGAL 4.10
  typedef Barycentric_coordinates Barycentric_coordinate;
  typedef Construct_barycentric_coordinates_weight Construct_barycentric_coordinate_weight;
  typedef Construct_barycentric_coordinates Construct_barycentric_coordinate;
  typedef Construct_barycentric_coordinates_in_triangle_2 Construct_barycentric_coordinate_in_triangle_2;
  typedef Construct_barycentric_coordinates_in_triangle_3 Construct_barycentric_coordinate_in_triangle_3;
  typedef Classify_barycentric_coordinates Classify_barycentric_coordinate;

  Construct_barycentric_coordinates construct_barycentric_coordinate_object() const { return m_construct_barycentric_coordinates_object; }
  Construct_barycentric_coordinates_weight construct_barycentric_coordinate_weight_object() const { return m_construct_barycentric_coordinates_weight_object; }
  Classify_barycentric_coordinates classify_barycentric_coordinate_object() const { return m_classify_barycentric_coordinates_object; }
  Construct_barycentric_coordinates_in_triangle_2 construct_barycentric_coordinate_in_triangle_2_object() const { return m_construct_barycentric_coordinates_in_triangle_2_object; }
  Construct_barycentric_coordinates_in_triangle_3 construct_barycentric_coordinate_in_triangle_3_object() const { return m_construct_barycentric_coordinates_in_triangle_3_object; }
  #endif

};

template <class Kernel, class Triangle_mesh>
std::ostream& operator<<(std::ostream& os, typename Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh>::Barycentric_coordinates b)
{
  return os << b[0] << " " << b[1] << " " << b[2];
}

#ifdef CGAL_SMSP_USE_ROBUST_TRAITS_CODE
#if defined(CGAL_USE_LEDA) || defined(CGAL_USE_CORE)
#ifndef DOXYGEN_RUNNING // needed due to a bug in doxygen
/*!
\ingroup PkgSurfaceMeshShortestPathTraitsClasses

\internal

\brief Provides an implementation of the SurfaceMeshShortestPathTraits
model which uses an exact Kernel during the unfolding operations to achieve better overall precision

\tparam K Kernel Type

\tparam TriangleMesh triangle mesh type

\cgalModels `SurfaceMeshShortestPathTraits`
*/
template <
  class K,
  class TriangleMesh>
class Surface_mesh_shortest_path_traits_with_robust_unfolding : public Surface_mesh_shortest_path_traits<K,TriangleMesh>
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
    : Surface_mesh_shortest_path_traits<K,TriangleMesh>(kernel)
    , m_robust_construct_triangle_3_to_triangle_2_projection_object(kernel)
    , m_robust_flatten_triangle_3_along_segment_2(kernel)
  {
  }

  Construct_triangle_3_to_triangle_2_projection construct_triangle_3_to_triangle_2_projection_object() const { return m_robust_construct_triangle_3_to_triangle_2_projection_object; }
  Construct_triangle_3_along_segment_2_flattening construct_triangle_3_along_segment_2_flattening_object() const { return m_robust_flatten_triangle_3_along_segment_2; }
};
#endif
#endif
#endif

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_TRAITS_H
