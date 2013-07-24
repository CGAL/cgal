#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <vector>
#include <CGAL/assertions.h>
#include <CGAL/internal/Triangulate_hole_polyline.h>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Timer.h>
namespace CGAL {

/**
 * @brief Function triangulating a hole in surface mesh.
 * 
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam OutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 *
 * @param[in, out] polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] output iterator over patch facets.
 * 
 * @return `output`
 */
template<class Polyhedron, class OutputIterator>
OutputIterator 
triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  OutputIterator output
  )
{
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator   Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Vertex_handle                      Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle                    Halfedge_handle;
  typedef typename Polyhedron::Traits::Point_3                    Point_3;

  std::vector<Point_3>         P, Q;
  std::vector<Halfedge_handle> P_edges;
  std::set<Vertex_handle>      vertex_set;

  Halfedge_around_facet_circulator circ(border_halfedge), done(circ);
  do{
    P.push_back(circ->vertex()->point());
    Q.push_back(circ->next()->opposite()->next()->vertex()->point());
    P_edges.push_back(circ);
    vertex_set.insert(circ->vertex());
  } while (++circ != done);

  CGAL::Timer timer; timer.start();

  // fill hole using polyline function
  std::vector<boost::tuple<int, int, int> > tris;
  triangulate_hole_polyline(P.begin(), P.end(), Q.begin(), Q.end(), back_inserter(tris));

  CGAL_TRACE_STREAM << "Hole filling: " << timer.time() << " sc." << std::endl; timer.reset();

  // construct patch
  // TODO: performance might be improved,
  //       although test on RedCircleBox.off shows that cost of this part is really insignificant compared to hole filling
  polyhedron.fill_hole(border_halfedge);
  *output++ = border_halfedge->facet();
  for(std::vector<boost::tuple<int, int, int> >::iterator tris_it = tris.begin(); tris_it != tris.end(); ++tris_it) {

    for(int i = 0; i < 3; ++i) { // iterate over each edge in triangle
      int v0_id, v1_id;
      if(i == 0)      { v0_id = tris_it->get<0>(); v1_id = tris_it->get<1>(); }
      else if(i == 1) { v0_id = tris_it->get<1>(); v1_id = tris_it->get<2>(); }
      else            { v0_id = tris_it->get<2>(); v1_id = tris_it->get<0>(); }

      bool border_edge = (v0_id + 1 == v1_id) || (v0_id == 0 && v1_id == P.size() -1); 
      bool smaller = v0_id < v1_id;

      if( smaller &&  /* to process one edge just one time */  
         !border_edge /* border edges don't require split  */ ) 
      { // this part is going to be executed exactly N-3 times

        Vertex_handle v0 = P_edges[v0_id]->vertex();
        Vertex_handle v1 = P_edges[v1_id]->vertex();

        Halfedge_around_vertex_circulator circ_v0(v0->vertex_begin()), circ_v1(v1->vertex_begin());
        bool found = false;
        // find halfedges sharing the common facet for split
        do { 
          circ_v1 = v1->vertex_begin();
          do {
            if(circ_v0->facet() == circ_v1->facet()) 
            {
              if(vertex_set.find(circ_v0->opposite()->vertex()) != vertex_set.end() &&
                 vertex_set.find(circ_v1->opposite()->vertex()) != vertex_set.end() )
              { found = true; }
            }
          } while(!found && ++circ_v1 != v1->vertex_begin());
        } while(!found && ++circ_v0 != v0->vertex_begin());

        CGAL_assertion(found);
        polyhedron.split_facet(circ_v0, circ_v1);
        *output++ = circ_v1->facet();
      }
    }
  }
  CGAL_TRACE_STREAM << "Patch construction: " << timer.time() << " sc." << std::endl;
  return output;
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H