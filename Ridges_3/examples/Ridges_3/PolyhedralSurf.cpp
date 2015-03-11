#include "PolyhedralSurf.h"

void PolyhedralSurf::compute_facets_normals()
{
  BOOST_FOREACH(face_descriptor f, faces(*this)){
    halfedge_descriptor h = halfedge(f,*this);
      Vector_3 normal =
	CGAL::cross_product(target(h,*this)->point() -
			    target(opposite(h,*this),*this)->point(),
			    target(next(h,*this),*this)->point() -
			    target(opposite(h,*this),*this)->point());
      f->setNormal( normal / CGAL::sqrt(normal * normal));
  }
}

const Vector_3 PolyhedralSurf::computeFacetsAverageUnitNormal(vertex_descriptor v)
{
  halfedge_descriptor h;
  face_descriptor f;
  Vector_3 sum(0., 0., 0.), n;

  CGAL::Halfedge_around_target_circulator<PolyhedralSurf> hedgeb(halfedge(v,*this),*this), hedgee = hedgeb;

  do
    {
      h = *hedgeb;
      if (is_border_edge(h,*this))
	{
	  hedgeb++;
	  continue;
	}

      f =  face(h,*this);
      n = f->getUnitNormal();
      sum = (sum + n);
      hedgeb++;
    }
  while (hedgeb != hedgee);
  sum = sum / std::sqrt(sum * sum);
  return sum;
}
