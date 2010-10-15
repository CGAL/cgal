#ifndef CGAL_POLY_OP_H_
#define CGAL_POLY_OP_H_

//the facet stores the normal
struct Facet_unit_normal {
  template < class Facet >
  void operator() (Facet & f)
    {
      typename Facet::Halfedge_handle h = f.halfedge();
      typename Facet::Vector_3 normal =
	CGAL::cross_product(h->vertex()->point() -
			    h->opposite()->vertex()->point(),
			    h->next()->vertex()->point() -
			    h->opposite()->vertex()->point());
      f.setNormal( normal / CGAL::sqrt(normal * normal));
    }
};

#endif
