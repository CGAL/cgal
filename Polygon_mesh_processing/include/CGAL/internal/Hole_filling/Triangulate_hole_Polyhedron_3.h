#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <CGAL/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/iterator.h>
#include <vector>

namespace CGAL {
namespace internal {


template<class PolygonMesh, class OutputIterator>
struct Tracer_polyhedron 
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor Halfedge_handle;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor    Facet_handle;

  Tracer_polyhedron(OutputIterator out,
                    PolygonMesh& poly,
                    std::vector<Halfedge_handle>& P)
    : out(out), poly(poly), P(P)
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
          Euler::fill_hole(h,poly); }
      else 
        { h = Euler::add_face_to_border(prev(P[i+1],poly), P[i+2/*k*/], poly); }
      
      CGAL_assertion(face(h,poly) != boost::graph_traits<PolygonMesh>::null_face());
      *out++ = face(h,poly);
      return opposite(h,poly);
    } 
    else 
    {
      int la = lambda.get(i, k);
      h = operator()(lambda, i, la, false);
      g = operator()(lambda, la, k, false);

      if(last)
        {
          h = g;
          Euler::fill_hole(g,poly);
        }
      else 
        { h = Euler::add_face_to_border(prev(h,poly), g, poly); }

      CGAL_assertion(face(h,poly) != boost::graph_traits<PolygonMesh>::null_face());
      *out++ = face(h,poly);
      return opposite(h,poly);
    }
  }

  OutputIterator out;
  PolygonMesh& poly;
  std::vector<Halfedge_handle>& P;
};

// This function is used in test cases (since it returns not just OutputIterator but also Weight)
template<class PolygonMesh, class OutputIterator>
std::pair<OutputIterator, Weight_min_max_dihedral_and_area> 
triangulate_hole_Polyhedron(PolygonMesh& poly, 
                            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, 
                            OutputIterator out,
                            bool use_delaunay_triangulation = false)
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
  Point_property_map ppmap = get(vertex_point,poly);

  int n = 0;
  Hedge_around_face_circulator circ(border_halfedge,poly), done(circ);
  do{
    P.push_back(ppmap[target(*circ, poly)]);
    Q.push_back(ppmap[target(next(opposite(next(*circ,poly),poly),poly),poly)]);
    P_edges.push_back(*circ);
    if(!vertex_set.insert(std::make_pair(target(*circ,poly), n++)).second) {
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

    Halfedge_around_target_circulator<PolygonMesh> circ_vertex(halfedge(v_it->first,poly),poly), done_vertex(circ_vertex);
    do {
      Vertex_set_it v_it_neigh_it = vertex_set.find(source(*circ_vertex,poly));

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
  
  // fill hole using polyline function, with custom tracer for PolygonMesh
  Tracer_polyhedron<PolygonMesh, OutputIterator>
    tracer(out, poly, P_edges);
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
 If no possible patch is found, @a poly is not altered in any way, and no face descriptor is put into @a out.

 @tparam PolygonMesh must be a model of `MutableFaceGraph`
 @tparam OutputIterator iterator holding `boost::graph_traits<PolygonMesh>::face_descriptor` for patch faces.

 @param poly surface mesh containing the hole
 @param border_halfedge a border halfedge incident to the hole
 @param out over patch faces.
 @param use_delaunay_triangulation

 @return @a out

  \todo move to a non-internal header file
  \todo handle islands
*/
template<class PolygonMesh, class OutputIterator>
OutputIterator 
triangulate_hole(PolygonMesh& poly, 
                 typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge, 
                 OutputIterator out,
                 bool use_delaunay_triangulation = false)
{
  CGAL_precondition(face(border_halfedge,poly) == boost::graph_traits<PolygonMesh>::null_face());
  return internal::triangulate_hole_Polyhedron
    (poly, border_halfedge, out, use_delaunay_triangulation).first;
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
