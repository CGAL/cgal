#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <CGAL/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>
#include <CGAL/boost/graph/iterator.h>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {


template<class PolygonMesh, class OutputIterator>
struct Tracer_polyhedron 
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_handle;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor    Facet_handle;

  Tracer_polyhedron(OutputIterator out,
                    PolygonMesh& pmesh,
                    std::vector<Halfedge_handle>& P)
    : out(out), pmesh(pmesh), P(P)
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
        {
          h = P[i+1];
          Euler::fill_hole(h,pmesh); }
      else 
        { h = Euler::add_face_to_border(prev(P[i+1],pmesh), P[i+2/*k*/], pmesh); }
      
      CGAL_assertion(face(h,pmesh) != boost::graph_traits<PolygonMesh>::null_face());
      *out++ = face(h,pmesh);
      return opposite(h,pmesh);
    } 
    else 
    {
      int la = lambda.get(i, k);
      h = operator()(lambda, i, la, false);
      g = operator()(lambda, la, k, false);

      if(last)
        {
          h = g;
          Euler::fill_hole(g,pmesh);
        }
      else 
        { h = Euler::add_face_to_border(prev(h,pmesh), g, pmesh); }

      CGAL_assertion(face(h,pmesh) != boost::graph_traits<PolygonMesh>::null_face());
      *out++ = face(h,pmesh);
      return opposite(h,pmesh);
    }
  }

  OutputIterator out;
  PolygonMesh& pmesh;
  std::vector<Halfedge_handle>& P;
};

// This function is used in test cases (since it returns not just OutputIterator but also Weight)
template<class PolygonMesh, class OutputIterator>
std::pair<OutputIterator, CGAL::internal::Weight_min_max_dihedral_and_area> 
triangulate_hole_Polyhedron(PolygonMesh& pmesh, 
                            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, 
                            OutputIterator out,
                            bool use_delaunay_triangulation)
{
  typedef Halfedge_around_face_circulator<PolygonMesh>   Hedge_around_face_circulator;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor Vertex_handle;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_handle;
  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point_3;
  
  typedef typename std::map<Vertex_handle, int>::iterator         Vertex_set_it;

  std::vector<Point_3>         P, Q;
  std::vector<Halfedge_handle> P_edges;
  std::map<Vertex_handle, int> vertex_set;
  Point_property_map ppmap = get(vertex_point,pmesh);

  int n = 0;
  Hedge_around_face_circulator circ(border_halfedge,pmesh), done(circ);
  do{
    P.push_back(ppmap[target(*circ, pmesh)]);
    Q.push_back(ppmap[target(next(opposite(next(*circ,pmesh),pmesh),pmesh),pmesh)]);
    P_edges.push_back(*circ);
    if(!vertex_set.insert(std::make_pair(target(*circ,pmesh), n++)).second) {
      CGAL_warning(!"Returning no output. Non-manifold vertex is found on boundary!");
      return std::make_pair(out,
                            CGAL::internal::Weight_min_max_dihedral_and_area::NOT_VALID());
    }
  } while (++circ != done);
  
  // existing_edges contains neighborhood information between boundary vertices
  // more precisely if v_i is neighbor to any other vertex than v_(i-1) and v_(i+1),
  // this edge is put into existing_edges
  std::vector<std::pair<int, int> > existing_edges;
  for(Vertex_set_it v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
    int v_it_id = v_it->second;
    int v_it_prev = v_it_id == 0   ? n-1 : v_it_id-1;
    int v_it_next = v_it_id == n-1 ? 0   : v_it_id+1;

    Halfedge_around_target_circulator<PolygonMesh> circ_vertex(halfedge(v_it->first,pmesh),pmesh), done_vertex(circ_vertex);
    do {
      Vertex_set_it v_it_neigh_it = vertex_set.find(source(*circ_vertex,pmesh));

      if(v_it_neigh_it != vertex_set.end()) {
        int v_it_neigh_id = v_it_neigh_it->second;
        if( v_it_neigh_id != v_it_prev && v_it_neigh_id != v_it_next ) {
          if(v_it_id < v_it_neigh_id) // to include one edge just one time
          { existing_edges.push_back(std::make_pair(v_it_id, v_it_neigh_id)); }
        }
      }
    } while(++circ_vertex != done_vertex);
  }

  CGAL::Timer timer; timer.start();

  //#define CGAL_USE_WEIGHT_INCOMPLETE
  #ifdef CGAL_USE_WEIGHT_INCOMPLETE
  typedef CGAL::internal::Weight_calculator<Weight_incomplete<CGAL::internal::Weight_min_max_dihedral_and_area>,
        CGAL::internal::Is_valid_existing_edges_and_degenerate_triangle> WC;
  #else
  typedef CGAL::internal::Weight_calculator<CGAL::internal::Weight_min_max_dihedral_and_area,
        CGAL::internal::Is_valid_existing_edges_and_degenerate_triangle> WC;
  #endif

  CGAL::internal::Is_valid_existing_edges_and_degenerate_triangle iv(existing_edges);
  
  // fill hole using polyline function, with custom tracer for PolygonMesh
  Tracer_polyhedron<PolygonMesh, OutputIterator>
    tracer(out, pmesh, P_edges);
  CGAL::internal::Weight_min_max_dihedral_and_area weight = 
        triangulate_hole_polyline(P, Q, 
                                        tracer, WC(iv), use_delaunay_triangulation)
                                      #ifdef CGAL_USE_WEIGHT_INCOMPLETE
                                        .weight; // get actual weight in Weight_incomplete
                                      #else
                                        ;
                                      #endif

  CGAL_TRACE_STREAM << "Hole filling: " << timer.time() << " sc." << std::endl; timer.reset();
  return std::make_pair(tracer.out, weight);
}

}// namespace internal
}// namespace Polygon_mesh_processing
}// namespace CGAL
#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
