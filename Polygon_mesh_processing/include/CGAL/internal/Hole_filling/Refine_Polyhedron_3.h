#ifndef CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
#include <cmath>
#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>

namespace CGAL {
namespace internal {

template<class Polyhedron>
class Refine_Polyhedron_3 {
// typedefs
  typedef typename boost::property_map<Polyhedron,vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point_3;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor Vertex_handle;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor Halfedge_handle;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor Facet_handle;

  typedef Halfedge_around_face_circulator<Polyhedron>   Halfedge_around_facet_circulator;
  typedef Halfedge_around_target_circulator<Polyhedron>  Halfedge_around_vertex_circulator;

private:
  Polyhedron& poly;
  Point_property_map ppmap;
  
  bool flippable(Halfedge_handle h) {
    // this check is added so that edge flip does not break manifoldness
    // it might happen when there is an edge where flip_edge(h) will be placed (i.e. two edges collide after flip)
    Vertex_handle v_tip_0 = target(next(h,poly),poly);
    Vertex_handle v_tip_1 = target(next(opposite(h,poly),poly),poly);
    Halfedge_around_vertex_circulator v_cir(next(h,poly), poly), v_end(v_cir);
    do {
      if(target(opposite(*v_cir, poly),poly) == v_tip_1) { return false; }
    } while(++v_cir != v_end);
    
    // also eliminate collinear triangle generation
    if( CGAL::collinear(ppmap[v_tip_0], ppmap[v_tip_1], ppmap[target(h, poly)]) ||
        CGAL::collinear(ppmap[v_tip_0], ppmap[v_tip_1], ppmap[target(opposite(h, poly),poly)]) ) {
      return false;
    }

    return true;
  }

  bool relax(Halfedge_handle h)
  {
    const Point_3& p = ppmap[target(h,poly)];
    const Point_3& q = ppmap[target(opposite(h,poly),poly)];
    const Point_3& r = ppmap[target(next(h,poly),poly)];
    const Point_3& s = ppmap[target(next(opposite(h,poly),poly),poly)];
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
        (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) )
    {
      if(flippable(h)) {
        Euler::flip_edge(h,poly);
        return true;
      }
    }
    return false;
  }

  template<class VertexOutputIterator, class FacetOutputIterator>
  bool subdivide(std::vector<Facet_handle>& facets, 
                 const std::set<Halfedge_handle>& border_edges,
                 std::map<Vertex_handle, double>& scale_attribute, 
                 VertexOutputIterator& vertex_out,
                 FacetOutputIterator& facet_out,
                 double alpha)
  {
    std::size_t facet_size = facets.size();
    for(std::size_t i = 0; i < facet_size; ++i){
      CGAL_assertion(facets[i]  != Facet_handle());

      Halfedge_handle hh =  halfedge(facets[i],poly);
      Vertex_handle vi = target(halfedge(facets[i],poly),poly);
      Vertex_handle vj = target(next(halfedge(facets[i],poly),poly),poly);
      Vertex_handle vk = target(prev(halfedge(facets[i],poly),poly),poly);
      Point_3 c = CGAL::centroid(ppmap[vi], ppmap[vj], ppmap[vk]);
      double sac  = (scale_attribute[vi] + scale_attribute[vj] + scale_attribute[vk])/3.0;
      double dist_c_vi = std::sqrt(CGAL::squared_distance(c,ppmap[vi]));
      double dist_c_vj = std::sqrt(CGAL::squared_distance(c,ppmap[vj]));
      double dist_c_vk = std::sqrt(CGAL::squared_distance(c,ppmap[vk]));
      if((alpha * dist_c_vi > sac) &&
         (alpha * dist_c_vj > sac) &&
         (alpha * dist_c_vk > sac) &&
         (alpha * dist_c_vi > scale_attribute[vi]) &&
         (alpha * dist_c_vj > scale_attribute[vj]) &&
         (alpha * dist_c_vk > scale_attribute[vk])){
        Halfedge_handle h = Euler::add_center_vertex(halfedge(facets[i],poly),poly);
        put(ppmap, target(h,poly), c);
          scale_attribute[target(h,poly)] = sac;
          *vertex_out++ = target(h,poly);

          // collect 2 new facets for next round 
          Facet_handle h1 = face(opposite(next(h,poly),poly),poly);
          Facet_handle h2 = face(opposite(h,poly),poly);
          facets.push_back(h1); facets.push_back(h2);
          *facet_out++ = h1;    *facet_out++ = h2;
          // relax edges of the  patching mesh 
          Halfedge_handle e_ij = prev(h,poly);
          Halfedge_handle e_ik = next(opposite(h,poly),poly);
          Halfedge_handle e_jk = prev(opposite(next(h,poly),poly),poly);

          if(border_edges.find(e_ij) == border_edges.end()){
            relax(e_ij);
          }
          if(border_edges.find(e_ik) == border_edges.end()){
            relax(e_ik);
          }
          if(border_edges.find(e_jk) == border_edges.end()){
            relax(e_jk);
          }
      }
    }
    return facets.size() != facet_size;
  }

  bool relax(const std::vector<Facet_handle>& facets, 
             const std::set<Halfedge_handle>& border_edges)
  {
    int flips = 0;
    std::list<Halfedge_handle> interior_edges;
    std::set<Halfedge_handle> included_map; 

    for(typename std::vector<Facet_handle>::const_iterator it = facets.begin(); it!= facets.end(); ++it) {
      Halfedge_around_face_circulator<Polyhedron>  circ(halfedge(*it,poly),poly), done(circ);
      do {
        Halfedge_handle h = *circ;        
        if(border_edges.find(h) == border_edges.end()){
          // do not remove included_map and use if(&*h < &*oh) { interior_edges.push_back(h) } 
          // which will change the order of edges from run to run
          Halfedge_handle oh = opposite(h,poly);
          Halfedge_handle h_rep = (h < oh) ? h : oh; // AF: was &*h < &*oh
          if(included_map.insert(h_rep).second) {
            interior_edges.push_back(h_rep);
          }
        }
      } while(++circ != done);
    }

    CGAL_TRACE_STREAM << "Test " << interior_edges.size() << " edges " << std::endl;
    //do not just use std::set (included_map) for iteration, the order effects the output (we like to make it deterministic)
    for(typename std::list<Halfedge_handle>::iterator it = interior_edges.begin(); it != interior_edges.end();++it) {
      if(relax(*it)) {
        ++flips;
      }
    }

    CGAL_TRACE_STREAM << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

  double average_length(Vertex_handle vh,
                        const std::set<Facet_handle>& interior_map, 
                        bool accept_internal_facets)
  {
    const Point_3& vp = ppmap[vh]; 
    Halfedge_around_target_circulator<Polyhedron> circ(halfedge(vh,poly),poly), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      Facet_handle f(face(*circ,poly)), f_op(face(opposite(*circ,poly),poly));

      if(!accept_internal_facets) {
        if(interior_map.find(f) != interior_map.end() && interior_map.find(f_op) != interior_map.end())
        { continue; } // which means current edge is an interior edge and should not be included in scale attribute calculation
      }

      const Point_3& vq = ppmap[target(opposite(*circ,poly),poly)];
      sum += std::sqrt(CGAL::squared_distance(vp, vq));
      ++deg;
    } while(++circ != done);

    CGAL_assertion(deg != 0); // this might happen when accept_internal_facets = false but there is
    return sum/deg;
  }

  void calculate_scale_attribute(const std::vector<Facet_handle>& facets, 
                                 const std::set<Facet_handle>& interior_map,
                                 std::map<Vertex_handle, double>& scale_attribute,
                                 bool accept_internal_facets) 
  {
    for(typename std::vector<Facet_handle>::const_iterator f_it = facets.begin(); f_it != facets.end(); ++f_it) {
      Halfedge_around_face_circulator<Polyhedron> circ(halfedge(*f_it,poly),poly), done(circ);
      do {
        Vertex_handle v = target(*circ,poly);
        std::pair<typename std::map<Vertex_handle, double>::iterator, bool> v_insert 
          = scale_attribute.insert(std::make_pair(v, 0));
        if(!v_insert.second) { continue; } // already calculated
        v_insert.first->second = average_length(v, interior_map, accept_internal_facets);
      } while(++circ != done);
    }
  }

  bool contain_internal_facets(const std::vector<Facet_handle>& facets,
                               const std::set<Facet_handle>& interior_map) const
  {
    for(typename std::vector<Facet_handle>::const_iterator f_it = facets.begin(); f_it != facets.end(); ++f_it) {
      Halfedge_around_face_circulator<Polyhedron> circ(halfedge(*f_it,poly),poly), done(circ);
      do {
        Vertex_handle v = target(*circ,poly);
        Halfedge_around_target_circulator<Polyhedron> circ_v(*circ,poly), done_v(circ_v);
        bool internal_v = true;
        do {
          Facet_handle f(face(*circ,poly)), f_op(face(opposite(*circ_v,poly),poly));

          if(interior_map.find(f) == interior_map.end() || interior_map.find(f_op) == interior_map.end()) {
            internal_v = false;
            break;
          } 
        } while(++circ_v != done_v);

        if(internal_v) { return true; }
      } while(++circ != done);
    }
    return false;
  }

public:
  Refine_Polyhedron_3(Polyhedron& poly)
    : poly(poly), ppmap(get(vertex_point, poly))
  {}

  template<class InputIterator, class FacetOutputIterator, class VertexOutputIterator>
  void refine(InputIterator facet_begin, 
              InputIterator facet_end, 
              FacetOutputIterator& facet_out,
              VertexOutputIterator& vertex_out,
              double alpha)
  {
    std::vector<Facet_handle> facets(facet_begin, facet_end); // do not use just std::set, the order effects the output (for the same input we want to get same output)
    std::set<Facet_handle> interior_map(facet_begin, facet_end);

    // store boundary edges - to be used in relax 
    std::set<Halfedge_handle> border_edges;
    for(typename std::vector<Facet_handle>::const_iterator it = facets.begin(); it!= facets.end(); ++it){
      Halfedge_around_face_circulator<Polyhedron>  circ(halfedge(*it,poly),poly), done(circ);
      do {
        if(interior_map.find(face(opposite(*circ,poly),poly)) == interior_map.end()) {
          border_edges.insert(*circ);
        }
      } while(++circ != done);
    }

    // check whether there is any need to accept internal facets
    bool accept_internal_facets = contain_internal_facets(facets, interior_map);
    std::map<Vertex_handle, double> scale_attribute;
    calculate_scale_attribute(facets, interior_map, scale_attribute, accept_internal_facets);

    CGAL::Timer total_timer; total_timer.start();
    for(int i = 0; i < 10; ++i) {
      CGAL::Timer timer; timer.start();
      bool is_subdivided = subdivide(facets, border_edges, scale_attribute, vertex_out, facet_out, alpha);
      CGAL_TRACE_STREAM << "**Timer** subdivide() :" << timer.time() << std::endl; timer.reset();
      if(!is_subdivided) { break; }

      bool is_relaxed = relax(facets, border_edges);
      CGAL_TRACE_STREAM << "**Timer** relax() :" << timer.time() << std::endl;
      if(!is_relaxed) { break; }
    }

    CGAL_TRACE_STREAM << "**Timer** TOTAL: " << total_timer.time() << std::endl;
  }
};

}//namespace internal

/*!
\ingroup PkgPolygonMeshProcessing
@brief Function refining a region on surface mesh

@tparam Polyhedron a %CGAL polyhedron
@tparam InputIterator iterator over input facets
@tparam FacetOutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets
@tparam VertexOutputIterator iterator holding `Polyhedron::Vertex_handle` for patch vertices

@param poly surface mesh to be refined
@param facet_begin first iterator of the range of facets
@param facet_end past-the-end iterator of the range of facets
@param facet_out iterator over newly created facets
@param vertex_out iterator over newly created vertices
@param density_control_factor factor for density where larger values cause denser refinements

@return pair of @a facet_out and @a vertex_out

@todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
\todo move to a non-internal header file
\todo `Polyhedron` should be a model of `MutableFaceGraph`
 */
template<
  class Polyhedron,
  class InputIterator,
  class FacetOutputIterator,
  class VertexOutputIterator
>
std::pair<FacetOutputIterator, VertexOutputIterator>
refine(Polyhedron& poly,
       InputIterator facet_begin, 
       InputIterator facet_end,
       FacetOutputIterator facet_out,
       VertexOutputIterator vertex_out,
       double density_control_factor = std::sqrt(2.0))
{
  internal::Refine_Polyhedron_3<Polyhedron> refine_functor(poly);
  refine_functor.refine
    (facet_begin, facet_end, facet_out, vertex_out, density_control_factor);
  return std::make_pair(facet_out, vertex_out);
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
