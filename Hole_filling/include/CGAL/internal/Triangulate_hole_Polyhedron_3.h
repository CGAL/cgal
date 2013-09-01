#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <CGAL/internal/Triangulate_hole_polyline.h>
#include <CGAL/Timer.h>

#include <vector>

namespace CGAL {
namespace internal {


template<class Polyhedron, class OutputIterator>
struct Tracer_polyhedron 
{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle    Facet_handle;

  Tracer_polyhedron(OutputIterator& out,
                    Polyhedron& polyhedron,
                    std::vector<Halfedge_handle>& P)
    : out(out), polyhedron(polyhedron), P(P)
  { }

  template <class LookupTable>
  Halfedge_handle 
  operator()(const LookupTable& lambda, 
             int i, int k,
             bool last = true)
  {
    if(i + 1 == k) { return P[i+1]; }

    Halfedge_handle h, g;
    if(i+2 == k){
      if(last)
      { h = polyhedron.fill_hole(P[i+1]); }
      else 
      { h = polyhedron.add_facet_to_border(P[i+1]->prev(), P[i+2/*k*/]); }
      
      CGAL_assertion(h->facet() != Facet_handle());
      *out++ = h->facet();
      return h->opposite();
    } 
    else 
    {
      int la = lambda.get(i, k);
      h = operator()(lambda, i, la, false);
      g = operator()(lambda, la, k, false);

      if(last)
      { h = polyhedron.fill_hole(g); }
      else 
      { h = polyhedron.add_facet_to_border(h->prev(), g); }

      CGAL_assertion(h->facet() != Facet_handle());
      *out++ = h->facet();
      return h->opposite();
    }
  }

  OutputIterator& out;
  Polyhedron& polyhedron;
  std::vector<Halfedge_handle>& P;
};

template<class Polyhedron, class OutputIterator>
std::pair<OutputIterator, Weight_min_max_dihedral_and_area> 
triangulate_hole_Polyhedron(Polyhedron& polyhedron, 
                            typename Polyhedron::Halfedge_handle border_halfedge, 
                            OutputIterator out,
                            bool use_delaunay_triangulation = false)
{
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator   Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle                    Halfedge_handle;
  typedef typename Polyhedron::Traits::Point_3                    Point_3;

  typedef typename std::map<Vertex_handle, int>::iterator         Vertex_set_it;

  std::vector<Point_3>         P, Q;
  std::vector<Halfedge_handle> P_edges;
  std::map<Vertex_handle, int> vertex_set;
  
  int n = 0;
  Halfedge_around_facet_circulator circ(border_halfedge), done(circ);
  do{
    P.push_back(circ->vertex()->point());
    Q.push_back(circ->next()->opposite()->next()->vertex()->point());
    P_edges.push_back(circ);
    if(!vertex_set.insert(std::make_pair(circ->vertex(), n++)).second) {
      CGAL_warning(!"Returning no output. Non-manifold vertex is found on boundary!");
      return std::make_pair(out, Weight_min_max_dihedral_and_area::NOT_VALID());
    }
  } while (++circ != done);
  
  // existing_edges contains neighborhood information between boundary vertices
  // more precisely if v_i is neighbor to any other vertex than v_(i-1) and v_(i+1),
  // this edge is put into existing_edges
  Edge_set existing_edges;
  for(Vertex_set_it v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
    int v_it_id = v_it->second;
    int v_it_prev = v_it_id == 0   ? n-1 : v_it_id-1;
    int v_it_next = v_it_id == n-1 ? 0   : v_it_id+1;

    Halfedge_around_vertex_circulator circ_vertex(v_it->first->vertex_begin()), done_vertex(circ_vertex);
    do {
      Vertex_set_it v_it_neigh_it = vertex_set.find(circ_vertex->opposite()->vertex());

      if(v_it_neigh_it != vertex_set.end()) {
        int v_it_neigh_id = v_it_neigh_it->second;
        if( v_it_neigh_id != v_it_prev && v_it_neigh_id != v_it_next ) {
          existing_edges.insert(std::make_pair(v_it_id, v_it_neigh_id));
        }
      }
    } while(++circ_vertex != done_vertex);
  }

  CGAL::Timer timer; timer.start();

  typedef Weight_calculator<Weight_min_max_dihedral_and_area, Is_valid_existing_edges_and_degenerate_triangle> WC;
  Is_valid_existing_edges_and_degenerate_triangle iv(existing_edges);
  
  // fill hole using polyline function, with custom tracer for Polyhedron
  Tracer_polyhedron<Polyhedron, OutputIterator> tracer(out, polyhedron, P_edges);
  Weight_min_max_dihedral_and_area weight = 
    internal::triangulate_hole_polyline(P.begin(), P.end(), 
                                        Q.begin(), Q.end(), 
                                        tracer, WC(iv), use_delaunay_triangulation);

  CGAL_TRACE_STREAM << "Hole filling: " << timer.time() << " sc." << std::endl; timer.reset();
  return std::make_pair(tracer.out, weight);
}

}// namespace internal

/*!
 \ingroup PkgHoleFilling
 Function triangulating a hole in surface mesh.
 The hole should contain no non-manifold vertex. Generated patch is guaranteed to not to break edge manifoldness and contain no degenerate triangle.
 If no possible patch is found, @a polyhedron is not altered in any way, and no facet handle is put into @a out.

 @tparam Polyhedron a %CGAL polyhedron
 @tparam OutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets.

 @param polyhedron surface mesh containing the hole
 @param border_halfedge a border halfedge incident to the hole
 @param output iterator over patch facets.
 @param use_delaunay_triangulation

 @return @a out
*/
template<class Polyhedron, class OutputIterator>
OutputIterator 
triangulate_hole(Polyhedron& polyhedron, 
                 typename Polyhedron::Halfedge_handle border_halfedge, 
                 OutputIterator out,
                 bool use_delaunay_triangulation = false)
{
  CGAL_precondition(border_halfedge->is_border());
  return internal::triangulate_hole_Polyhedron
    (polyhedron, border_halfedge, out, use_delaunay_triangulation).first;
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H