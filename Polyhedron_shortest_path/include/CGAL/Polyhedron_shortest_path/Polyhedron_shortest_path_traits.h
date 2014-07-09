// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H

#include <CGAL/Polyhedron_shortest_path/Internal/Barycentric.h>
#include <CGAL/Polyhedron_shortest_path/Internal/function_objects.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

namespace CGAL {

/*!
\ingroup PkgPolyhedronShortestPathTraitsClasses

\brief Provides an implementation of the PolyhedronShortestPathTraits 
model as required by the Polyhedron_shortest_path algorithm

\tparam K The kernel type whose geometric primitives to use

\tparam P The polyhedron type the algorithm is to act on

\cgalModels `PolyhedronShortestPathTraits`
*/
template <
  class K, 
  class P>
class Polyhedron_shortest_path_default_traits
{
public:

  typedef K Kernel;
  typedef P Polyhedron;

  typedef typename Kernel::FT FT;
  
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Triangle_2 Triangle_2;
  
  typedef typename Kernel::Vector_3 Barycentric_coordinate;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  
  // Predicates
public:
  typedef typename Kernel::Orientation_2 Orientation_2;
  typedef typename internal::Compare_relative_intersection_along_segment_2<K> Compare_relative_intersection_along_segment_2;
  typedef typename internal::Is_saddle_vertex<K,P> Is_saddle_vertex;
  
  // Constructions
public:
  typedef typename Kernel::Intersect_2 Intersect_2;
  typedef typename Kernel::Construct_projected_point_2 Construct_projected_point_2;
  typedef typename Kernel::Compute_squared_distance_2 Compute_squared_distance_2;
  typedef typename Kernel::Compute_squared_distance_3 Compute_squared_distance_3;
  typedef typename internal::Project_triangle_3_to_triangle_2<K> Project_triangle_3_to_triangle_2;
  typedef typename internal::Flatten_triangle_3_along_segment_2<K> Flatten_triangle_3_along_segment_2;
  typedef typename internal::Parametric_distance_along_segment_2<K> Parametric_distance_along_segment_2;

  class Construct_barycentric_coordinate_2
  {
  public:
    Barycentric_coordinate operator()(const Triangle_2& t, const Point_2& p) const
    {
      return internal::Construct_barycentric_coordinate_any<FT, Barycentric_coordinate, Point_2, Vector_2, Triangle_2>(t, p);
    }
  };
  
  class Construct_triangle_location_2
  {
  public:
    Point_2 operator()(const Triangle_2& t, const Barycentric_coordinate& a) const
    {
      return internal::Construct_triangle_location_any<FT, Barycentric_coordinate, Point_2, Triangle_2>(t, a);
    }
  };
  
  class Construct_barycentric_coordinate_3
  {
  public:
    Barycentric_coordinate operator()(const Triangle_3& t, const Point_3& p) const
    {
      return internal::Construct_barycentric_coordinate_any<FT, Barycentric_coordinate, Point_3, Vector_3, Triangle_3>(t, p);
    }
  };
  
  class Construct_triangle_location_3
  {
  public:
    Point_3 operator()(const Triangle_3& t, const Barycentric_coordinate& a) const
    {
      return internal::Construct_triangle_location_any<FT, Barycentric_coordinate, Point_3, Triangle_3>(t, a);
    }
  };
  
private:
  Kernel m_kernel;
  Project_triangle_3_to_triangle_2 m_project_triangle_3_to_triangle_2_object;
  Flatten_triangle_3_along_segment_2 m_flatten_triangle_3_along_segment_2_object;
  Construct_barycentric_coordinate_2 m_construct_barycentric_coordinate_2_object;
  Construct_triangle_location_2 m_construct_triangle_location_2_object;
  Construct_barycentric_coordinate_3 m_construct_barycentric_coordinate_3_object;
  Construct_triangle_location_3 m_construct_triangle_location_3_object;
  Compare_relative_intersection_along_segment_2 m_compare_relative_intersection_along_segment_2_object;
  Is_saddle_vertex m_is_saddle_vertex_object;
  Parametric_distance_along_segment_2 m_parametric_distance_along_segment_2_object;
  
public:

  Polyhedron_shortest_path_default_traits()
  {
  }
  
  Polyhedron_shortest_path_default_traits(const Kernel& kernel)
    : m_kernel(kernel)
    , m_project_triangle_3_to_triangle_2_object(m_kernel.compute_squared_distance_3_object())
    , m_flatten_triangle_3_along_segment_2_object(m_kernel.compute_squared_distance_3_object())
    , m_is_saddle_vertex_object(m_project_triangle_3_to_triangle_2_object, m_flatten_triangle_3_along_segment_2_object, m_kernel.orientation_2_object())
    , m_compare_relative_intersection_along_segment_2_object(m_kernel.compute_squared_distance_2_object(), m_kernel.intersect_2_object())
  {
  }
  
  Orientation_2 orientation_2_object() const { return m_kernel.orientation_2_object(); }
  Intersect_2 intersect_2_object() const { return m_kernel.intersect_2_object(); }
  Construct_projected_point_2 construct_projected_point_2_object() const { return m_kernel.construct_projected_point_2_object(); }
  Compute_squared_distance_2 compute_squared_distance_2_object() const { return m_kernel.compute_squared_distance_2_object(); }
  Compute_squared_distance_3 compute_squared_distance_3_object() const { return m_kernel.compute_squared_distance_3_object(); }
  Compare_relative_intersection_along_segment_2 compare_relative_intersection_along_segment_2_object() const { return m_compare_relative_intersection_along_segment_2_object; }
  Is_saddle_vertex is_saddle_vertex_object() const { return m_is_saddle_vertex_object; }
  
  Project_triangle_3_to_triangle_2 project_triangle_3_to_triangle_2_object() const { return m_project_triangle_3_to_triangle_2_object; }
  Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2_object() const { return m_flatten_triangle_3_along_segment_2_object; }
  Construct_barycentric_coordinate_2 construct_barycentric_coordinate_2_object() const { return m_construct_barycentric_coordinate_2_object; }
  Construct_triangle_location_2 construct_triangle_location_2_object() const { return m_construct_triangle_location_2_object; }
  Construct_barycentric_coordinate_3 construct_barycentric_coordinate_3_object() const { return m_construct_barycentric_coordinate_3_object; }
  Construct_triangle_location_3 construct_triangle_location_3_object() const { return m_construct_triangle_location_3_object; }
  Parametric_distance_along_segment_2 parametric_distance_along_segment_2_object() const { return m_parametric_distance_along_segment_2_object; }
};

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H
