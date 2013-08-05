namespace CGAL {
namespace internal {

template<class Polyhedron, class edge_descriptor, class Point>
edge_descriptor mesh_split(Polyhedron& polyhedron, edge_descriptor ei, Point pn)
{
  edge_descriptor en = polyhedron.split_edge(ei);
  en->vertex()->point() = pn;
  polyhedron.split_facet(en, ei->next());

  en->id() = -1;
  en->opposite()->id() = -1;
  ei->id() = -1;
  ei->opposite()->id() = -1;
  en->next()->id() = -1;
  en->next()->opposite()->id() = -1;
  en->next()->next()->id() = -1;
  ei->next()->id() = -1;
  edge_descriptor ej = en->opposite();
  if (!(ej->is_border()))
  {
    polyhedron.split_facet(ei->opposite(), ej->next());
    ej->next()->id() = -1;
    edge_descriptor ei_op_next = ei->opposite()->next();
    ei_op_next->id() = -1;
    ei_op_next->opposite()->id() = -1;
    ei_op_next->next()->id() = -1;
  }

  return en;
}

template<class Polyhedron>
typename Polyhedron::Traits::FT volume(const Polyhedron& p)
{
  typedef typename Polyhedron::Facet_const_iterator Facet_iterator;
  typedef typename Polyhedron::Traits::Point_3 Point;
  typename Kernel_traits<Point>::Kernel::Compute_volume_3 volume;
  typename Polyhedron::Traits::FT res=0;
  Point origin(0,0,0);
  for (Facet_iterator f = p.facets_begin(); f != p.facets_end(); ++f)
  {
    typename Polyhedron::Halfedge_const_handle he = f->facet_begin();
    const Point& v1 = he->vertex()->point();
    const Point& v2 = he->next()->vertex()->point();
    const Point& v3 = he->next()->next()->vertex()->point();
    res+=volume(origin,v1,v2,v3);
  }
  return res;
}

} //namespace internal
} //namespace CGAL
