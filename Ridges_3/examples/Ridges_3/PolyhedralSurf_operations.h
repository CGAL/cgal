#ifndef CGAL_POLY_OP_H_
#define CGAL_POLY_OP_H_

//the facet stores the normal
template <typename TriangleMesh>
struct Facet_unit_normal {
  TriangleMesh * P;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor Halfedge_handle;
  Facet_unit_normal(TriangleMesh* P)
    : P(P)
  {}

  template < class Facet >
  void operator() (Facet & f)
    {
      typename Halfedge_handle h = halfedge(f,*P);
      typename Facet::Vector_3 normal =
	CGAL::cross_product(target(h,*P)->point() -
			    target(opposite(h,*P),*P)->point(),
			    target(next(h,*P),*P)->point() -
			    target(opposite(h,*P),*P)->point());
      f.setNormal( normal / CGAL::sqrt(normal * normal));
    }
};

#endif
