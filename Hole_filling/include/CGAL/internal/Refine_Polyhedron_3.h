#ifndef CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
#include <cmath>
#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/trace.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL {
namespace internal {

//// use custom comp to prevent random ordering in set
//template<class Polyhedron>
//struct Facet_comparator {
//  typedef typename Polyhedron::Facet_handle Facet_handle;
//  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
//  typedef typename Polyhedron::Vertex_handle Vertex_handle;
//  typedef typename Polyhedron::Traits::Point_3 Point_3;
//
//  bool operator()(Facet_handle f1, Facet_handle f2) const {
//    return &*f1 < &*f2;
//    //const Point_3& p11 = f1->halfedge()->vertex()->point();
//    //const Point_3& p21 = f1->halfedge()->next()->vertex()->point();
//    //const Point_3& p31 = f1->halfedge()->next()->next()->vertex()->point();
//
//    //const Point_3& p12 = f2->halfedge()->vertex()->point();
//    //const Point_3& p22 = f2->halfedge()->next()->vertex()->point();
//    //const Point_3& p32 = f2->halfedge()->next()->next()->vertex()->point();
//
//    //if(p11 != p12) { return point_comp(p11, p12); }
//    //if(p21 != p22) { return point_comp(p21, p22); }
//    //return point_comp(p31, p32);
//
//    //Halfedge_handle hf1 = f1->halfedge(), hf2 = f2->halfedge();
//    //do {
//    //  if(hf1->vertex()->point() != hf2->vertex()->point()) {
//    //    return point_comp(hf1->vertex()->point(), hf2->vertex()->point());
//    //  }
//    //  ++hf1; ++hf2;
//    //} while(hf1 != f1->halfedge());
//    //return false;
//  }
//
//  bool point_comp(const Point_3& p1, const Point_3& p2) const {
//    if(p1.x() != p2.x()) return p1.x() < p2.x();
//    if(p1.y() != p2.y()) return p1.y() < p2.y();
//    return p1.z() < p2.z();
//  }
//};

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

  //typedef std::set<Facet_handle, Facet_comparator<Polyhedron> > Facet_set;
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

        //if(!poly.is_valid()) { 
        //  std::cout << "before flip not valid" << std::endl; 
        //}

        poly.flip_edge(h);

        //if(!poly.is_valid()) {
        //  std::cout << "after flip not valid" << std::endl; 
        //}

        return true;
    }
    return false;
  }

  bool subdivide(Polyhedron& poly, std::vector<Facet_handle>& facets, 
    std::set<Facet_handle>& interior_map, std::map<Vertex_handle, double>& scale_attribute)
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

          // collect 2 new facets for next round 
          Facet_handle h1 = h->next()->opposite()->face();
          Facet_handle h2 = h->opposite()->face();
          new_facets.push_back(h1); interior_map.insert(h1);
          new_facets.push_back(h2); interior_map.insert(h2);
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

  bool relax(Polyhedron& poly, std::vector<Facet_handle>& facets, std::set<Facet_handle>& interior_map)
  {
    int flips = 0;
    std::list<Halfedge_handle> interior_edges;
    std::set<Halfedge_handle> included_map; // do not use just std::set, the order effects the output (for the same input we want to get same output)

    for(typename std::vector<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      Halfedge_around_facet_circulator  circ = (*it)->facet_begin(), done(circ);
      do {
        Halfedge_handle h = circ;
        Halfedge_handle oh = h->opposite();
        if(interior_map.find(oh->face()) != interior_map.end()){
          // it's an interior edge
          Halfedge_handle h_rep = (h < oh) ? h : oh;
          if(included_map.insert(h_rep).second) {
            interior_edges.push_back(h_rep);
          }
        }
        ++circ;
      } while(circ != done);
    }

    CGAL_TRACE_STREAM << "Test " << interior_edges.size() << " edges " << std::endl;
    for(typename std::list<Halfedge_handle>::iterator it = interior_edges.begin();
      it != interior_edges.end();
      ++it){
        if(relax(poly,*it)){
          ++flips;
        }
    }

    CGAL_TRACE_STREAM << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

public:
  template<class InputIterator, class OutputIterator>
  void operator()(Polyhedron& poly, std::map<Vertex_handle, double>& scale_attribute, 
    InputIterator facet_begin, InputIterator facet_end, OutputIterator out)
  {
    std::vector<Facet_handle> facets(facet_begin, facet_end); // do not use just std::set, the order effects the output (for the same input we want to get same output)
    std::set<Facet_handle> interior_map(facet_begin, facet_end);

    int i = 0;
    do {
      if(i == 10){
        break;
      }
      i++;
    } while( subdivide(poly, facets, interior_map, scale_attribute) && relax(poly,facets, interior_map) );

    std::copy(facets.begin(), facets.end(), out);
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
#endif //CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H