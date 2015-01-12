#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <CGAL/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/iterator.h>
#include <vector>

namespace CGAL {
namespace internal {


template<class Polyhedron, class OutputIterator>
struct Tracer_polyhedron 
{
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor Halfedge_handle;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor    Facet_handle;

  Tracer_polyhedron(OutputIterator out,
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
        { h = polyhedron.add_facet_to_border(prev(P[i+1],polyhedron), P[i+2/*k*/]); }
      
      CGAL_assertion(face(h,polyhedron) != boost::graph_traits<Polyhedron>::null_face());
      *out++ = face(h,polyhedron);
      return opposite(h,polyhedron);
    } 
    else 
    {
      int la = lambda.get(i, k);
      h = operator()(lambda, i, la, false);
      g = operator()(lambda, la, k, false);

      if(last)
      { h = polyhedron.fill_hole(g); }
      else 
        { h = polyhedron.add_facet_to_border(prev(h,polyhedron), g); }

      CGAL_assertion(face(h,polyhedron) != boost::graph_traits<Polyhedron>::null_face());
      *out++ = face(h,polyhedron);
      return opposite(h,polyhedron);
    }
  }

  OutputIterator out;
  Polyhedron& polyhedron;
  std::vector<Halfedge_handle>& P;
};

// This function is used in test cases (since it returns not just OutputIterator but also Weight)
template<class Polyhedron, class OutputIterator>
std::pair<OutputIterator, Weight_min_max_dihedral_and_area> 
triangulate_hole_Polyhedron(Polyhedron& polyhedron, 
                            typename boost::graph_traits<Polyhedron>::halfedge_descriptor border_halfedge, 
                            OutputIterator out,
                            bool use_delaunay_triangulation = false)
{
  typedef Halfedge_around_target_circulator<Polyhedron>  Halfedge_around_vertex_circulator;
  typedef Halfedge_around_face_circulator<Polyhedron>   Halfedge_around_face_circulator;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor Vertex_handle;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor Halfedge_handle;
  typedef typename boost::property_map<Polyhedron,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point_3;
  
  typedef typename std::map<Vertex_handle, int>::iterator         Vertex_set_it;

  std::vector<Point_3>         P, Q;
  std::vector<Halfedge_handle> P_edges;
  std::map<Vertex_handle, int> vertex_set;
  Point_property_map ppmap = get(vertex_point,polyhedron);

  int n = 0;
  Halfedge_around_face_circulator circ(border_halfedge,polyhedron), done(circ);
  do{
    P.push_back(ppmap[target(*circ, polyhedron)]);
    Q.push_back(ppmap[target(next(opposite(next(*circ,polyhedron),polyhedron),polyhedron),polyhedron)]);
    P_edges.push_back(*circ);
    if(!vertex_set.insert(std::make_pair(target(*circ,polyhedron), n++)).second) {
      CGAL_warning(!"Returning no output. Non-manifold vertex is found on boundary!");
      return std::make_pair(out, Weight_min_max_dihedral_and_area::NOT_VALID());
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

    Halfedge_around_target_circulator<Polyhedron> circ_vertex(v_it->first->vertex_begin(),polyhedron), done_vertex(circ_vertex);
    do {
      Vertex_set_it v_it_neigh_it = vertex_set.find(source(*circ_vertex,polyhedron));

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
  typedef Weight_calculator<Weight_incomplete<Weight_min_max_dihedral_and_area>, Is_valid_existing_edges_and_degenerate_triangle> WC;
  #else
  typedef Weight_calculator<Weight_min_max_dihedral_and_area, Is_valid_existing_edges_and_degenerate_triangle> WC;
  #endif

  Is_valid_existing_edges_and_degenerate_triangle iv(existing_edges);
  
  // fill hole using polyline function, with custom tracer for Polyhedron
  Tracer_polyhedron<Polyhedron, OutputIterator>
    tracer(out, polyhedron, P_edges);
  Weight_min_max_dihedral_and_area weight = 
    internal::triangulate_hole_polyline(P, Q, 
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

/*!
 \ingroup PkgPolygonMeshProcessing
 Function triangulating a hole in surface mesh.
 The hole should contain no non-manifold vertex. Generated patch is guaranteed to not to break edge manifoldness and contain no degenerate triangle.
 If no possible patch is found, @a polyhedron is not altered in any way, and no face descriptor is put into @a out.

 @tparam Polyhedron a %CGAL polyhedron
 @tparam OutputIterator iterator holding `boost::graph_traits<Polyhedron>::face_descriptor` for patch faces.

 @param polyhedron surface mesh containing the hole
 @param border_halfedge a border halfedge incident to the hole
 @param out over patch faces.
 @param use_delaunay_triangulation

 @return @a out

  \todo move to a non-internal header file
  \todo handle islands
  \todo `Polyhedron` should be a model of `MutableFaceGraph`
*/
template<class Polyhedron, class OutputIterator>
OutputIterator 
triangulate_hole(Polyhedron& polyhedron, 
                 typename boost::graph_traits<Polyhedron>::halfedge_descriptor border_halfedge, 
                 OutputIterator out,
                 bool use_delaunay_triangulation = false)
{
  CGAL_precondition(face(border_halfedge,polyhedron) == boost::graph_traits<Polyhedron>::null_face());
  return internal::triangulate_hole_Polyhedron
    (polyhedron, border_halfedge, out, use_delaunay_triangulation).first;
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
