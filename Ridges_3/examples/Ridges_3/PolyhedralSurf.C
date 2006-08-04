#include "PolyhedralSurf.h"

// Vector_3 PolyhedralSurf::getHalfedge_vector(Halfedge * h)
// {
//   Vector_3 v = h->opposite()->vertex()->point() - h->vertex()->point();
//   return v;
// }

// double PolyhedralSurf::
// compute_mean_edges_length_around_vertex(Vertex* v)
// {
//   Halfedge_around_vertex_circulator
//     hedgeb = v->vertex_begin(), hedgee = hedgeb;
//   int  count_he = 0;
//   double sum = 0.;
//   CGAL_For_all(hedgeb, hedgee)
//     {
//       sum += hedgeb->getLength();
//       count_he++;
//     }
//   return sum/count_he;
// }

// void PolyhedralSurf::compute_edges_length()
// {
//   std::for_each(this->halfedges_begin(), this->halfedges_end(),
// 		Edge_length());
// }


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

