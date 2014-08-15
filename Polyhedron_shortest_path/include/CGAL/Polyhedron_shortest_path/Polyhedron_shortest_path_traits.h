// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_TRAITS_H

#include <CGAL/Polyhedron_shortest_path/barycentric.h>
#include <CGAL/Polyhedron_shortest_path/function_objects.h>

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
class Polyhedron_shortest_path_default_traits : public K
{
public:

  typedef K Kernel;
  typedef F FaceGraph;

  typedef typename Kernel::FT FT;
  
  class Barycentric_coordinate
  {
  private:
    CGAL::cpp11::array<FT,3> m_coords;
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
  
  // Predicates
public:
  typedef typename PolyhedronShortestPath::Compare_relative_intersection_along_segment_2<Kernel> Compare_relative_intersection_along_segment_2;
  typedef typename PolyhedronShortestPath::Is_saddle_vertex<Kernel, FaceGraph> Is_saddle_vertex;
  
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

  typedef typename PolyhedronShortestPath::Project_triangle_3_to_triangle_2<K> Project_triangle_3_to_triangle_2;
  typedef typename PolyhedronShortestPath::Flatten_triangle_3_along_segment_2<K> Flatten_triangle_3_along_segment_2;
  typedef typename PolyhedronShortestPath::Parametric_distance_along_segment_2<K> Parametric_distance_along_segment_2;
  typedef typename PolyhedronShortestPath::Construct_barycentric_coordinate_in_triangle_2<K, Barycentric_coordinate, Construct_barycentric_coordinate> Construct_barycentric_coordinate_in_triangle_2;
  typedef typename PolyhedronShortestPath::Construct_barycentric_coordinate_in_triangle_3<K, Barycentric_coordinate, Construct_barycentric_coordinate> Construct_barycentric_coordinate_in_triangle_3;
  typedef typename PolyhedronShortestPath::Classify_barycentric_coordinate<Barycentric_coordinate, Construct_barycentric_coordinate_weight> Classify_barycentric_coordinate;
  
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
    , m_project_triangle_3_to_triangle_2_object(m_kernel)
    , m_flatten_triangle_3_along_segment_2_object(m_kernel)
    , m_is_saddle_vertex_object(m_kernel, m_project_triangle_3_to_triangle_2_object, m_flatten_triangle_3_along_segment_2_object)
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
