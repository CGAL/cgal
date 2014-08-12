#ifndef _GET_NORMAL_
#define _GET_NORMAL_

namespace CGAL {

namespace internal {

template <class Facet, class Kernel>
typename Kernel::Vector_3 get_facet_normal(const Facet& f)
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

template <class Vertex, class Kernel>
typename Kernel::Vector_3 get_vertex_normal(const Vertex& v)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Vertex::Halfedge_around_vertex_const_circulator HV_circulator;
  Vector normal = CGAL::NULL_VECTOR;
  HV_circulator he = v.vertex_begin();
  HV_circulator end = he;
  CGAL_For_all(he,end)
  {
    if(!he->is_border())
    {
      const Point& prev = he->prev()->vertex()->point();
      const Point& curr = he->vertex()->point();
      const Point& next = he->next()->vertex()->point();

      Vector p1 = next - curr;
      p1 = p1 / std::sqrt(p1 * p1);
      Vector p2 = prev - curr;
      p2 = p2 / std::sqrt(p2 * p2);

      double cosine = p1 * p2;
      if      (cosine < -1.0) cosine = -1.0;
      else if (cosine >  1.0) cosine =  1.0;
      double angle = acos(cosine);

      Vector n = CGAL::cross_product(next-curr,prev-curr);
      n = n / std::sqrt(n * n);
      n = n * angle;

//      Vector n = get_facet_normal<Facet,Kernel>(*he->facet());
      normal = normal + n;
    }
  }
  return normal / std::sqrt(normal * normal);
}

} //namespace internal

} //namespace CGAL

#endif // _GET_NORMAL_
