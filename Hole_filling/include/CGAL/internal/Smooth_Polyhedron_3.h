#ifndef CGAL_HOLE_FILLING_SMOOTH_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_SMOOTH_POLYHEDRON_3_H
// going to be moved some place more relevant
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <CGAL/internal/Weights.h>

namespace CGAL {

// simple Laplacian smoothing
template<class WeightCalculator, class Polyhedron, class InputIterator>
void smooth(  Polyhedron& polyhedron,
  InputIterator vertex_begin, 
  InputIterator vertex_end, 
  WeightCalculator weight_calculator = WeightCalculator()
  ) 
{
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator in_edge_iterator;
  typedef typename Polyhedron::Traits::Point_3   Point;

  for(; vertex_begin !=vertex_end; ++vertex_begin) {
    vertex_descriptor v = *vertex_begin;
    
    double total_w = 0.0;
    Point smoothed(0.0, 0.0, 0.0);
    in_edge_iterator e, e_end;
    for (boost::tie(e,e_end) = boost::in_edges(v, polyhedron); e != e_end; e++)
    {
      double w_ij = weight_calculator(*e, polyhedron);
      total_w += w_ij;
      smoothed = smoothed + ((boost::source(*e, polyhedron)->point() - CGAL::ORIGIN) * w_ij);
    }
    v->point() =  CGAL::ORIGIN + ((smoothed - CGAL::ORIGIN) / total_w);
  }
}

template<class Polyhedron, class InputIterator>
void smooth(Polyhedron& polyhedron,
  InputIterator vertex_begin, 
  InputIterator vertex_end
  ) 
{
  CGAL::smooth(polyhedron, vertex_begin, vertex_end,
    CGAL::internal::Cotangent_weight<Polyhedron, CGAL::internal::Cotangent_value_Meyer<Polyhedron> >());
}

}//namespace CGAL

#endif