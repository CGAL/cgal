#include <algorithm>
#include <CGAL/circulator.h>

namespace CGAL {
namespace internal {

template<unsigned int axis>
struct Axis_compare {
  template<class Vertex>
  bool operator()(const Vertex& v0, const Vertex& v1) const
  { return v0.point()[axis] < v1.point()[axis]; }
};

// Taken from compute_normal.h inside Polyhedron demo //
template <class Facet, class Kernel>
typename Kernel::Vector_3 compute_facet_normal(const Facet& f)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Facet::Halfedge_around_facet_const_circulator HF_circulator;
  Vector normal = CGAL::NULL_VECTOR;
  HF_circulator he = f.facet_begin();
  HF_circulator end = he;
  CGAL_For_all(he,end)
  {
    const Point& prev = he->prev()->vertex()->point();
    const Point& curr = he->vertex()->point();
    const Point& next = he->next()->vertex()->point();
    Vector n = CGAL::cross_product(next-curr,prev-curr);
    normal = normal + n;
  }
  return normal / std::sqrt(normal * normal);
}

template <class Kernel, class Vertex>
typename Kernel::Vector_3 compute_vertex_normal(const Vertex& v)
{
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Vertex::Halfedge_around_vertex_const_circulator HV_circulator;
  typedef typename Vertex::Facet Facet;
  Vector normal = CGAL::NULL_VECTOR;
  HV_circulator he = v.vertex_begin();
  HV_circulator end = he;
  CGAL_For_all(he,end)
  {
    if(!he->is_border())
    {
      Vector n = compute_facet_normal<Facet, Kernel>(*he->facet());
      normal = normal + (n / std::sqrt(n*n));
    }
  }
  return normal; // No need to normalize here
}

} // namespace internal

/** 
 * Test whether a polyhedron has correct orientation
 * @pre @a polyhedron.is_closed()
 *
 * @tparam Polyhedron a %CGAL polyhedron
 *
 * @param polyhedron a closed polyhedron to be tested
 *
 * @return true if orientation is OK
 * @code
 * // use inside out to fix orientation
 * if(!is_oriented(polyhedron)) {
 *   polyhedron.inside_out();
 * }
 * @endcode
 */
template<class Polyhedron>
bool is_oriented(const Polyhedron& polyhedron) {
  CGAL_assertion(polyhedron.is_closed());
  const unsigned int axis = 0;

  typename Polyhedron::Vertex_const_iterator v_min 
    = std::min_element(polyhedron.vertices_begin(), polyhedron.vertices_end(), internal::Axis_compare<axis>());
  
  typedef typename Polyhedron::Traits K;
  const typename K::Vector_3& normal_v_min = internal::compute_vertex_normal<K>(*v_min);

  CGAL_warning(normal_v_min[axis] != 0);
  return normal_v_min[axis] < 0;
} 
} // namespace CGAL