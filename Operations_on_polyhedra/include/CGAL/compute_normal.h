#ifndef CGAL_COMPUTE_NORMAL
#define CGAL_COMPUTE_NORMAL
#include <CGAL/circulator.h>

namespace CGAL {

namespace internal {

template<bool normalize>
struct Compute_normal { };

// Specialize for false
template<>
struct Compute_normal<false> {

  template <class Kernel, class Facet>
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
    return normal;
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
        Vector n = compute_facet_normal<Kernel>(*he->facet());
        normal = normal + (n / std::sqrt(n*n));
      }
    }
    return normal;
  }
};

// Specialize for true
template<>
struct Compute_normal<true> {

  template <class Kernel, class Facet>
  typename Kernel::Vector_3 compute_facet_normal(const Facet& f)
  {
    const typename Kernel::Vector_3& n = Compute_normal<false>().compute_facet_normal<Kernel>(f);
    return n / std::sqrt(n * n);
  }

  template <class Kernel, class Vertex>
  typename Kernel::Vector_3 compute_vertex_normal(const Vertex& v)
  {
    const typename Kernel::Vector_3& n = Compute_normal<false>().compute_vertex_normal<Kernel>(v);
    return n / std::sqrt(n * n);
  }
};

} // namespace internal

/** 
 * Computes facet normal.
 * @tparam Kernel kernel
 * @tparam normalize tag specifying whether the returned vector is normalized. It is default to false and can be omitted.
 * @tparam Facet type of facet
 *
 * @param f facet
 */
template <class Kernel, bool normalize, class Facet>
typename Kernel::Vector_3 compute_facet_normal(const Facet& f)
{ 
  return internal::Compute_normal<normalize>().compute_facet_normal<Kernel>(f);
}

// normalize = false
template <class Kernel, class Facet>
typename Kernel::Vector_3 compute_facet_normal(const Facet& f)
{ 
  return internal::Compute_normal<false>().compute_facet_normal<Kernel>(f);
}

/** 
 * Computes vertex normal by averaging normalized normals of neighbor facets.
 * @tparam Kernel kernel
 * @tparam normalize tag specifying whether the returned vector is normalized. It is default to false and can be omitted.
 *         Regardless of this tag normals of neighbor facets are normalized while averaging.
 * @tparam Vertex type of vertex
 *
 * @param v vertex
 */
template <class Kernel, bool normalize, class Vertex>
typename Kernel::Vector_3 compute_vertex_normal(const Vertex& v)
{
  return internal::Compute_normal<normalize>().compute_vertex_normal<Kernel>(v);
}

// normalize = false
template <class Kernel, class Vertex>
typename Kernel::Vector_3 compute_vertex_normal(const Vertex& v)
{
  return internal::Compute_normal<false>().compute_vertex_normal<Kernel>(v);
}

} // namespace CGAL
#endif // _COMPUTE_NORMAL_
