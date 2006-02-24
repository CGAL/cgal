#include "PolyhedralSurf.h"



Vector_3 PolyhedralSurf::getHalfedge_vector(Halfedge * h)
{
  Vector_3 v = h->opposite()->vertex()->point() - h->vertex()->point();
  return v;
}

///???????myterious behavior????????
void PolyhedralSurf::compute_facets_normals()
{
  std::for_each(this->facets_begin(), this->facets_end(),
		Facet_unit_normal()); 
}



Vector_3 PolyhedralSurf::computeFacetsAverageUnitNormal(Vertex * v)
{
  Halfedge *h;
  Facet *f;
  Vector_3 sum(0., 0., 0.), n;

  Halfedge_around_vertex_circulator
    hedgeb = v->vertex_begin(), hedgee = hedgeb;

  do
    {
      h = &(*hedgeb);
      if (h->is_border_edge())
	{
	  hedgeb++;
	  continue;
	}

      f = &(*h->facet());
      n = f->getUnitNormal();
      sum = (sum + n);
      hedgeb++;
    }
  while (hedgeb != hedgee);
  sum = sum / std::sqrt(sum * sum);
  return sum;
}


