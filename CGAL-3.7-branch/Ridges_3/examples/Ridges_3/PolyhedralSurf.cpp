#include "PolyhedralSurf.h"

void PolyhedralSurf::compute_facets_normals()
{
  std::for_each(this->facets_begin(), this->facets_end(),
		Facet_unit_normal());
}

const Vector_3 PolyhedralSurf::computeFacetsAverageUnitNormal(const Vertex_const_handle v)
{
  Halfedge_const_handle h;
  Facet_const_handle f;
  Vector_3 sum(0., 0., 0.), n;

  Halfedge_around_vertex_const_circulator
    hedgeb = v->vertex_begin(), hedgee = hedgeb;

  do
    {
      h = hedgeb;
      if (h->is_border_edge())
	{
	  hedgeb++;
	  continue;
	}

      f =  h->facet();
      n = f->getUnitNormal();
      sum = (sum + n);
      hedgeb++;
    }
  while (hedgeb != hedgee);
  sum = sum / std::sqrt(sum * sum);
  return sum;
}
