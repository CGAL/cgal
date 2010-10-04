#ifndef _COMPUTE_NORMAL_
#define _COMPUTE_NORMAL_

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

template <class Vertex, class Kernel>
typename Kernel::Vector_3 compute_vertex_normal(const Vertex& v)
{
  typedef typename Kernel::Point_3 Point;
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
      Vector n = compute_facet_normal<Facet,Kernel>(*he->facet());
      normal = normal + (n / std::sqrt(n*n));
    }
  }
  return normal / std::sqrt(normal * normal);
}

#endif // _COMPUTE_NORMAL_
