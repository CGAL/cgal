// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H

#include <CGAL/Polyhedron_shortest_path/internal/Barycentric.h>
#include <CGAL/Polyhedron_shortest_path/internal/function_objects.h>

#include <CGAL/boost/graph/properties.h>
//#include <CGAL/boost/graph/properties_Polyhedron_3.h>
//#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <ostream>
#include <boost/array.hpp>

namespace CGAL {

/*!
\ingroup PkgPolyhedronShortestPathTraitsClasses

\brief Provides an implementation of the FaceGraphShortestPathTraits 
model as required by the Polyhedron_shortest_path algorithm

\tparam K The kernel type whose geometric primitives to use

\tparam F The faceGraph type the algorithm is to act on

\cgalModels `FaceGraphShortestPathTraits`
*/
template <
  class K, 
  class F>
class Polyhedron_shortest_path_default_traits
{
public:

  typedef K Kernel;
  typedef F FaceGraph;

  typedef typename Kernel::FT FT;
  
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Triangle_2 Triangle_2;

  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  
  class Barycentric_coordinate
  {
  private:
    boost::array<FT,3> m_coords;
  public:
    Barycentric_coordinate()
    {
    }
    
    Barycentric_coordinate(const Barycentric_coordinate& other)
      : m_coords(other.m_coords)
    {
    }
  
    Barycentric_coordinate(const FT& a, const FT& b, const FT& c)
    {
      m_coords[0] = a;
      m_coords[1] = b;
      m_coords[2] = c;
    }
    
    const FT& operator [] (size_t i) const 
    {
      return m_coords[i % 3];
    }
  };
  
public:
  #define CGAL_kernel_functor(X, Y) \
    typedef typename Kernel::X X;     \
    X Y() const { return m_kernel.Y(); }

  // Predicates
  CGAL_kernel_functor(Compare_distance_2, compare_distance_2_object)
  CGAL_kernel_functor(Orientation_2, orientation_2_object);
    
  // Constructions
  CGAL_kernel_functor(Construct_point_2, construct_point_2_object)
  CGAL_kernel_functor(Construct_vector_2, construct_vector_2_object)
  CGAL_kernel_functor(Construct_ray_2, construct_ray_2_object)
  CGAL_kernel_functor(Construct_line_2, construct_line_2_object)
  CGAL_kernel_functor(Construct_segment_2, construct_segment_2_object)
  CGAL_kernel_functor(Construct_triangle_2, construct_triangle_2_object)

  CGAL_kernel_functor(Construct_vertex_2, construct_vertex_2_object)
  CGAL_kernel_functor(Construct_source_2, construct_source_2_object)
  CGAL_kernel_functor(Construct_target_2, construct_target_2_object)
  
  CGAL_kernel_functor(Construct_barycenter_2, construct_barycenter_2_object)
  CGAL_kernel_functor(Compute_squared_distance_2, compute_squared_distance_2_object)
  CGAL_kernel_functor(Intersect_2, intersect_2_object)
  CGAL_kernel_functor(Construct_point_on_2, construct_point_on_2_object)

  CGAL_kernel_functor(Construct_vector_3, construct_vector_3_object)
  CGAL_kernel_functor(Construct_triangle_3, construct_triangle_3_object)
  
  CGAL_kernel_functor(Construct_vertex_3, construct_vertex_3_object)
  CGAL_kernel_functor(Construct_source_3, construct_source_3_object)
  CGAL_kernel_functor(Construct_target_3, construct_target_3_object)
  
  CGAL_kernel_functor(Construct_barycenter_3, construct_barycenter_3_object)
  CGAL_kernel_functor(Compute_squared_distance_3, compute_squared_distance_3_object)
  
  #undef CGAL_kernel_functor
  
  // Predicates
public:
  typedef typename internal::Compare_relative_intersection_along_segment_2<Kernel> Compare_relative_intersection_along_segment_2;
  typedef typename internal::Is_saddle_vertex<Kernel, FaceGraph> Is_saddle_vertex;
  
  // Constructions
public:
  class Construct_barycentric_coordinate
  {
  public:
    typedef Barycentric_coordinate result_type;
    
    result_type operator() (const FT& a, const FT& b, const FT& c) const
    {
      return Barycentric_coordinate(a, b, c);
    }
  };
  
  class Construct_barycentric_coordinate_weight
  {
  public:
    typedef FT result_type;
    
    result_type operator() (const Barycentric_coordinate b, int i) const
    {
      return b[i % 3];
    }
  };

  typedef typename internal::Project_triangle_3_to_triangle_2<K> Project_triangle_3_to_triangle_2;
  typedef typename internal::Flatten_triangle_3_along_segment_2<K> Flatten_triangle_3_along_segment_2;
  typedef typename internal::Parametric_distance_along_segment_2<K> Parametric_distance_along_segment_2;
  typedef typename internal::Construct_barycentric_coordinate_in_triangle_2<K, Barycentric_coordinate, Construct_barycentric_coordinate> Construct_barycentric_coordinate_in_triangle_2;
  typedef typename internal::Construct_barycentric_coordinate_in_triangle_3<K, Barycentric_coordinate, Construct_barycentric_coordinate> Construct_barycentric_coordinate_in_triangle_3;
  typedef typename internal::Classify_barycentric_coordinate<Barycentric_coordinate, Construct_barycentric_coordinate_weight> Classify_barycentric_coordinate;
  
private:
  Kernel m_kernel;
  Construct_barycentric_coordinate m_construct_barycentric_coordinate_object;
  Construct_barycentric_coordinate_weight m_construct_barycentric_coordinate_weight_object;
  Classify_barycentric_coordinate m_classify_barycentric_coordinate_object;
  Project_triangle_3_to_triangle_2 m_project_triangle_3_to_triangle_2_object;
  Flatten_triangle_3_along_segment_2 m_flatten_triangle_3_along_segment_2_object;
  Construct_barycentric_coordinate_in_triangle_2 m_construct_barycentric_coordinate_in_triangle_2_object;
  Construct_barycentric_coordinate_in_triangle_3 m_construct_barycentric_coordinate_in_triangle_3_object;
  Compare_relative_intersection_along_segment_2 m_compare_relative_intersection_along_segment_2_object;
  Is_saddle_vertex m_is_saddle_vertex_object;
  Parametric_distance_along_segment_2 m_parametric_distance_along_segment_2_object;
  
public:

  Polyhedron_shortest_path_default_traits()
  {
  }
  
  Polyhedron_shortest_path_default_traits(const Kernel& kernel)
    : m_kernel(kernel)
    , m_classify_barycentric_coordinate_object(m_construct_barycentric_coordinate_weight_object)
    , m_project_triangle_3_to_triangle_2_object(m_kernel.compute_squared_distance_3_object())
    , m_flatten_triangle_3_along_segment_2_object(m_kernel.compute_squared_distance_3_object())
    , m_is_saddle_vertex_object(m_project_triangle_3_to_triangle_2_object, m_flatten_triangle_3_along_segment_2_object, m_kernel.orientation_2_object())
    , m_compare_relative_intersection_along_segment_2_object(m_kernel.compute_squared_distance_2_object(), m_kernel.intersect_2_object())
    , m_construct_barycentric_coordinate_in_triangle_2_object(m_construct_barycentric_coordinate_object, m_kernel.construct_difference_of_vectors_2_object(), m_kernel.compute_scalar_product_2_object())
    , m_construct_barycentric_coordinate_in_triangle_3_object(m_construct_barycentric_coordinate_object, m_kernel.construct_difference_of_vectors_3_object(), m_kernel.compute_scalar_product_3_object())
  {
  }

  Construct_barycentric_coordinate construct_barycentric_coordinate_object() const { return m_construct_barycentric_coordinate_object; }
  Construct_barycentric_coordinate_weight construct_barycentric_coordinate_weight_object() const { return m_construct_barycentric_coordinate_weight_object; }
  Classify_barycentric_coordinate classify_barycentric_coordinate_object() const { return m_classify_barycentric_coordinate_object; }
  
  Is_saddle_vertex is_saddle_vertex_object() const { return m_is_saddle_vertex_object; }
  
  Compare_relative_intersection_along_segment_2 compare_relative_intersection_along_segment_2_object() const { return m_compare_relative_intersection_along_segment_2_object; }
  Project_triangle_3_to_triangle_2 project_triangle_3_to_triangle_2_object() const { return m_project_triangle_3_to_triangle_2_object; }
  Flatten_triangle_3_along_segment_2 flatten_triangle_3_along_segment_2_object() const { return m_flatten_triangle_3_along_segment_2_object; }
  Construct_barycentric_coordinate_in_triangle_2 construct_barycentric_coordinate_in_triangle_2_object() const { return m_construct_barycentric_coordinate_in_triangle_2_object; }
  Construct_barycentric_coordinate_in_triangle_3 construct_barycentric_coordinate_in_triangle_3_object() const { return m_construct_barycentric_coordinate_in_triangle_3_object; }
  Parametric_distance_along_segment_2 parametric_distance_along_segment_2_object() const { return m_parametric_distance_along_segment_2_object; }
};

template <class Kernel, class FaceGraph>
std::ostream& operator<<(std::ostream& os, typename Polyhedron_shortest_path_default_traits<Kernel, FaceGraph>::Barycentric_coordinate b)
{
  return os << b[0] << " " << b[1] << " " << b[2];
}

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H
