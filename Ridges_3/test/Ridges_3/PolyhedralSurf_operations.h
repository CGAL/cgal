#ifndef _POLY_OP_H_
#define _POLY_OP_H_

/* struct Edge_length { */
/*   template < class HalfEdge >  */
/*   void operator() (HalfEdge & h)  */
/*     { */
/*       double d = */
/* 	CGAL::squared_distance(h.prev()->vertex()->point(), */
/* 			 h.vertex()->point()); */
/*       h.setLength(CGAL::sqrt(d)); */
/*   } */
/* }; */

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
