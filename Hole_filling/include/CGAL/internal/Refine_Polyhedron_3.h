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

namespace CGAL {
namespace internal {

template<class Polyhedron>
class Refine_Polyhedron_3 {
// typedefs
  typedef typename Polyhedron::Traits::Point_3   Point_3;
  typedef typename Polyhedron::Vertex_handle     Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle   Halfedge_handle;
  typedef typename Polyhedron::Facet_handle      Facet_handle;

  typedef typename Polyhedron::Halfedge_around_facet_circulator   Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

private:

  bool flippable(Halfedge_handle h) {
    // this check is added so that edge flip does not break manifoldness
    // it might happen when there is an edge where flip_edge(h) will be placed (i.e. two edges collide after flip)
    Vertex_handle v_tip_0 = h->next()->vertex();
    Vertex_handle v_tip_1 = h->opposite()->next()->vertex();
    Halfedge_around_vertex_circulator v_cir(v_tip_0->vertex_begin()), v_end(v_cir);
    do {
      if(v_cir->opposite()->vertex() == v_tip_1) { return false; }
    } while(++v_cir != v_end);
    
    // also eliminate collinear triangle generation
    if( CGAL::collinear(v_tip_0->point(), v_tip_1->point(), h->vertex()->point()) ||
        CGAL::collinear(v_tip_0->point(), v_tip_1->point(), h->opposite()->vertex()->point()) ) {
      return false;
    }

    return true;
  }

  bool relax(Polyhedron& poly, Halfedge_handle h)
  {
    const Point_3& p = h->vertex()->point();
    const Point_3& q = h->opposite()->vertex()->point();
    const Point_3& r = h->next()->vertex()->point();
    const Point_3& s = h->opposite()->next()->vertex()->point();
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
      (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) ){
      
      if(flippable(h)) {
        poly.flip_edge(h);
        return true;
      }
    }
    return false;
  }

  template<class VertexOutputIterator, class FacetOutputIterator>
  bool subdivide(Polyhedron& poly, 
                 std::vector<Facet_handle>& facets, 
                 std::set<Facet_handle>& interior_map,
                 std::map<Vertex_handle, double>& scale_attribute, 
                 VertexOutputIterator& vertex_out,
                 FacetOutputIterator& facet_out,
                 double alpha)
  {
    std::list<Facet_handle> new_facets;
    for(typename std::vector<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      CGAL_assertion(*it  != Facet_handle());

      Halfedge_handle hh =  (*it)->halfedge();
      Vertex_handle vi = (*it)->halfedge()->vertex();
      Vertex_handle vj = (*it)->halfedge()->next()->vertex();
      Vertex_handle vk = (*it)->halfedge()->prev()->vertex();
      Point_3 c = CGAL::centroid(vi->point(), vj->point(), vk->point());
      double sac  = (scale_attribute[vi] + scale_attribute[vj] + scale_attribute[vk])/3.0;
      double dist_c_vi = std::sqrt(CGAL::squared_distance(c,vi->point()));
      double dist_c_vj = std::sqrt(CGAL::squared_distance(c,vj->point()));
      double dist_c_vk = std::sqrt(CGAL::squared_distance(c,vk->point()));
      if((alpha * dist_c_vi > sac) &&
        (alpha * dist_c_vj > sac) &&
        (alpha * dist_c_vk > sac) &&
        (alpha * dist_c_vi > scale_attribute[vi]) &&
        (alpha * dist_c_vj > scale_attribute[vj]) &&
        (alpha * dist_c_vk > scale_attribute[vk])){
          Halfedge_handle h = poly.create_center_vertex((*it)->halfedge());
          h->vertex()->point() = c;
          scale_attribute[h->vertex()] = sac;
          *vertex_out++ = h->vertex();

          // collect 2 new facets for next round 
          Facet_handle h1 = h->next()->opposite()->face();
          Facet_handle h2 = h->opposite()->face();
          new_facets.push_back(h1); interior_map.insert(h1);
          new_facets.push_back(h2); interior_map.insert(h2);
          *facet_out++ = h1; *facet_out++ = h2;
          // relax edges of the  patching mesh 
          Halfedge_handle e_ij = h->prev();
          Halfedge_handle e_ik = h->opposite()->next();
          Halfedge_handle e_jk = h->next()->opposite()->prev();

          if(interior_map.find(e_ij->opposite()->face()) != interior_map.end()){
            relax(poly, e_ij);
          }
          if(interior_map.find(e_ik->opposite()->face()) != interior_map.end()){
            relax(poly, e_ik);
          }
          if(interior_map.find(e_jk->opposite()->face()) != interior_map.end()){
            relax(poly, e_jk);
          }
      }
    }
    facets.insert(facets.end(), new_facets.begin(), new_facets.end());
    return ! new_facets.empty();
  }

  bool relax(Polyhedron& poly, const std::vector<Facet_handle>& facets, const std::set<Facet_handle>& interior_map)
  {
    int flips = 0;
    std::list<Halfedge_handle> interior_edges;
    std::set<Halfedge_handle> included_map; 

    for(typename std::vector<Facet_handle>::const_iterator it = facets.begin(); it!= facets.end(); ++it) {
      Halfedge_around_facet_circulator  circ = (*it)->facet_begin(), done(circ);
      do {
        Halfedge_handle h = circ;
        Halfedge_handle oh = h->opposite();
        if(interior_map.find(oh->face()) != interior_map.end()){
          // do not remove included_map and use if(h < oh) { interior_edges.push_back(h) } 
          // which will change the order of edges from run to run
          Halfedge_handle h_rep = (h < oh) ? h : oh;
          if(included_map.insert(h_rep).second) {
            interior_edges.push_back(h_rep);
          }
        }
      } while(++circ != done);
    }

    CGAL_TRACE_STREAM << "Test " << interior_edges.size() << " edges " << std::endl;
    //do not just use std::set (included_map) for iteration, the order effects the output (we like to make it deterministic)
    for(typename std::list<Halfedge_handle>::iterator it = interior_edges.begin(); it != interior_edges.end();++it) {
      if(relax(poly,*it)) {
        ++flips;
      }
    }

    CGAL_TRACE_STREAM << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

  double average_length(Vertex_handle vh, const std::set<Facet_handle>& interior_map, bool accept_internal_facets)
  {
    const Point_3& vp = vh->point(); 
    Halfedge_around_vertex_circulator circ(vh->vertex_begin()), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      Facet_handle f(circ->facet()), f_op(circ->opposite()->facet());

      if(!accept_internal_facets) {
        if(interior_map.find(f) != interior_map.end() && interior_map.find(f_op) != interior_map.end())
        { continue; } // which means current edge is an interior edge and should not be included in scale attribute calculation
      }

      const Point_3& vq = circ->opposite()->vertex()->point();
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
      Halfedge_around_facet_circulator circ((*f_it)->facet_begin()), done(circ);
      do {
        Vertex_handle v = circ->vertex();
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
      Halfedge_around_facet_circulator circ((*f_it)->facet_begin()), done(circ);
      do {
        Vertex_handle v = circ->vertex();
        Halfedge_around_vertex_circulator circ_v(v->vertex_begin()), done_v(circ_v);
        bool internal_v = true;
        do {
          Facet_handle f(circ_v->facet()), f_op(circ_v->opposite()->facet());

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
  template<class InputIterator, class FacetOutputIterator, class VertexOutputIterator>
  void refine(Polyhedron& poly, 
              InputIterator facet_begin, 
              InputIterator facet_end, 
              FacetOutputIterator& facet_out,
              VertexOutputIterator& vertex_out,
              double alpha)
  {
    std::vector<Facet_handle> facets(facet_begin, facet_end); // do not use just std::set, the order effects the output (for the same input we want to get same output)
    std::set<Facet_handle> interior_map(facet_begin, facet_end);

    std::map<Vertex_handle, double> scale_attribute;
    bool accept_internal_facets = contain_internal_facets(facets, interior_map);
    calculate_scale_attribute(facets, interior_map, scale_attribute, accept_internal_facets);

    CGAL::Timer total_timer; total_timer.start();
    for(int i = 0; i < 10; ++i) {
      CGAL::Timer timer; timer.start();
      bool is_subdivided = subdivide(poly, facets, interior_map, scale_attribute, vertex_out, facet_out, alpha);
      CGAL_TRACE_STREAM << "**Timer** subdivide() :" << timer.time() << std::endl; timer.reset();
      if(!is_subdivided) { break; }

      bool is_relaxed = relax(poly,facets, interior_map);
      CGAL_TRACE_STREAM << "**Timer** relax() :" << timer.time() << std::endl;
      if(!is_relaxed) { break; }
    }

    CGAL_TRACE_STREAM << "**Timer** TOTAL: " << total_timer.time() << std::endl;
  }
};

}//namespace internal

/*!
\ingroup PkgHoleFilling
@brief Function refining a region on surface mesh

@tparam Polyhedron a %CGAL polyhedron
@tparam InputIterator iterator over input facets
@tparam FacetOutputIterator iterator holding `Polyhedron::Facet_handle` for patch facets
@tparam VertexOutputIterator iterator holding `Polyhedron::Vertex_handle` for patch vertices

@param polyhedron surface mesh to be refined
@param facet_begin first iterator of the range of facets
@param facet_end past-the-end iterator of the range of facets
@param facet_out iterator over newly created facets
@param vertex_out iterator over newly created vertices
@param density_control_factor factor for density where larger values cause denser refinements

@return pair of @a facet_out and @a vertex_out

@todo current algorithm iterates 10 times at most, since (I guess) there is no termination proof.
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
  internal::Refine_Polyhedron_3<Polyhedron> refine_functor;
  refine_functor.refine
    (poly, facet_begin, facet_end, facet_out, vertex_out, density_control_factor);
  return std::make_pair(facet_out, vertex_out);
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
