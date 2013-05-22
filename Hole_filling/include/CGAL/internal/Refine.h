#ifndef CGAL_HOLE_FILLING_REFINE_H
#define CGAL_HOLE_FILLING_REFINE_H
#include <cmath>
#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL {
namespace internal {

template<class Polyhedron>
class Refine_Polyhedron_3 {
// typedefs
  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
// members
  double alpha;

public:
  Refine_Polyhedron_3(double alpha) 
    : alpha(alpha) { }

private:
  bool relax(Polyhedron& poly, Halfedge_handle h)
  {
    const Point_3& p = h->vertex()->point();
    const Point_3& q = h->opposite()->vertex()->point();
    const Point_3& r = h->next()->vertex()->point();
    const Point_3& s = h->opposite()->next()->vertex()->point();
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
      (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) ){
        poly.flip_edge(h);
        return true;
    }
    return false;
  }

  bool subdivide(Polyhedron& poly, std::set<Facet_handle>& facets, std::map<Vertex_handle, double>& scale_attribute)
  {
    std::list<Facet_handle> new_facets;
    for(typename std::set<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      if(*it  == Facet_handle()){
      } else {
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

            // collect 2 new facets for next round 
            Facet_handle h1 = h->next()->opposite()->face();
            Facet_handle h2 = h->opposite()->face();
            new_facets.push_back(h1);
            new_facets.push_back(h2);
            // relax edges of the  patching mesh 
            Halfedge_handle e_ij = h->prev();
            Halfedge_handle e_ik = h->opposite()->next();
            Halfedge_handle e_jk = h->next()->opposite()->prev();

            if(facets.find(e_ij->opposite()->face()) != facets.end()){
              relax(poly, e_ij);
            }
            if(facets.find(e_ik->opposite()->face()) != facets.end()){
              relax(poly, e_ik);
            }
            if(facets.find(e_jk->opposite()->face()) != facets.end()){
              relax(poly, e_jk);
            }
        }
      }
    }
    for(typename std::list<Facet_handle>::iterator it = new_facets.begin();
      it!= new_facets.end();
      ++it){
        facets.insert(*it);
    }
    return ! new_facets.empty();
  }

  bool relax(Polyhedron& poly, std::set<Facet_handle>& facets)
  {
    int flips = 0;
    std::list<Halfedge_handle> interior_edges;
    for(typename std::set<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      Halfedge_around_facet_circulator  circ = (*it)->facet_begin(), done(circ);
      do {
        Halfedge_handle h = circ;
        Halfedge_handle oh = h->opposite();
        if(facets.find(oh->face()) != facets.end()){
          // it's an interior edge
          interior_edges.push_back((h < oh) ? h : oh );
        }
        ++circ;
      } while(circ != done);
    }
    std::cerr << "Test " << interior_edges.size() << " edges " << std::endl;
    for(typename std::list<Halfedge_handle>::iterator it = interior_edges.begin();
      it != interior_edges.end();
      ++it){
        if(relax(poly,*it)){
          ++flips;
        }
    }
    std::cerr << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

public:
  void operator()(std::map<Vertex_handle, double>& scale_attribute, std::set<Facet_handle>& facets, Polyhedron& poly)
  {
    int i = 0;
    do {
      if(i == 10){
        break;
      }
      i++;
    } while( subdivide(poly, facets, scale_attribute) && relax(poly,facets) );

    // according to paper it should be like below (?) IOY
    //while(true) {
    //  bool subdiv = subdivide(poly, facets);
    //  if(!subdiv) { break; }
    //  while(relax(poly,facets)) {}
    //}
  }
};

}//namespace internal
}//namespace CGAL
#endif //CGAL_HOLE_FILLING_REFINE_H