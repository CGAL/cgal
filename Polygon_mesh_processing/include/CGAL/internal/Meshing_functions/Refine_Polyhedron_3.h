#ifndef CGAL_POLYGON_MESH_PROCESSING_REFINE_POLYHEDRON_3_H
#define CGAL_POLYGON_MESH_PROCESSING_REFINE_POLYHEDRON_3_H

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
#include <CGAL/boost/graph/properties.h>
#include <boost/property_map/property_map.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template<class PolygonMesh>
class Refine_Polyhedron_3 {
//// typedefs
  typedef typename boost::property_map<PolygonMesh,
                                       boost::vertex_point_t>::type Point_property_map;
  typedef typename boost::property_traits<Point_property_map>::value_type Point_3;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  typedef Halfedge_around_face_circulator<PolygonMesh>   Halfedge_around_facet_circulator;
  typedef Halfedge_around_target_circulator<PolygonMesh>  Halfedge_around_vertex_circulator;

private:
  PolygonMesh& pmesh;
  Point_property_map ppmap;
  
  bool flippable(halfedge_descriptor h) {
    // this check is added so that edge flip does not break manifoldness
    // it might happen when there is an edge where flip_edge(h) will be placed (i.e. two edges collide after flip)
    vertex_descriptor v_tip_0 = target(next(h,pmesh),pmesh);
    vertex_descriptor v_tip_1 = target(next(opposite(h,pmesh),pmesh),pmesh);
    Halfedge_around_vertex_circulator v_cir(next(h,pmesh), pmesh), v_end(v_cir);
    do {
      if(target(opposite(*v_cir, pmesh),pmesh) == v_tip_1) { return false; }
    } while(++v_cir != v_end);
    
    // also eliminate collinear triangle generation
    if( CGAL::collinear(ppmap[v_tip_0], ppmap[v_tip_1], ppmap[target(h, pmesh)]) ||
        CGAL::collinear(ppmap[v_tip_0], ppmap[v_tip_1], ppmap[target(opposite(h, pmesh),pmesh)]) ) {
      return false;
    }

    return true;
  }

  bool relax(halfedge_descriptor h)
  {
    const Point_3& p = ppmap[target(h,pmesh)];
    const Point_3& q = ppmap[target(opposite(h,pmesh),pmesh)];
    const Point_3& r = ppmap[target(next(h,pmesh),pmesh)];
    const Point_3& s = ppmap[target(next(opposite(h,pmesh),pmesh),pmesh)];
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
        (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) )
    {
      if(flippable(h)) {
        Euler::flip_edge(h,pmesh);
        return true;
      }
    }
    return false;
  }

  template<class VertexOutputIterator, class FacetOutputIterator>
  bool subdivide(std::vector<face_descriptor>& facets, 
                 const std::set<halfedge_descriptor>& border_edges,
                 std::map<vertex_descriptor, double>& scale_attribute, 
                 VertexOutputIterator& vertex_out,
                 FacetOutputIterator& facet_out,
                 double alpha)
  {
    std::size_t facet_size = facets.size();
    for(std::size_t i = 0; i < facet_size; ++i){
      CGAL_assertion(facets[i]  != boost::graph_traits<PolygonMesh>::null_face());

      halfedge_descriptor hh =  halfedge(facets[i],pmesh);
      vertex_descriptor vi = target(halfedge(facets[i],pmesh),pmesh);
      vertex_descriptor vj = target(next(halfedge(facets[i],pmesh),pmesh),pmesh);
      vertex_descriptor vk = target(prev(halfedge(facets[i],pmesh),pmesh),pmesh);
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
        halfedge_descriptor h = Euler::add_center_vertex(halfedge(facets[i],pmesh),pmesh);
        put(ppmap, target(h,pmesh), c);
          scale_attribute[target(h,pmesh)] = sac;
          *vertex_out++ = target(h,pmesh);

          // collect 2 new facets for next round 
          face_descriptor h1 = face(opposite(next(h,pmesh),pmesh),pmesh);
          face_descriptor h2 = face(opposite(h,pmesh),pmesh);
          facets.push_back(h1); facets.push_back(h2);
          *facet_out++ = h1;    *facet_out++ = h2;
          // relax edges of the  patching mesh 
          halfedge_descriptor e_ij = prev(h,pmesh);
          halfedge_descriptor e_ik = next(opposite(h,pmesh),pmesh);
          halfedge_descriptor e_jk = prev(opposite(next(h,pmesh),pmesh),pmesh);

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

  bool relax(const std::vector<face_descriptor>& facets, 
             const std::set<halfedge_descriptor>& border_edges)
  {
    int flips = 0;
    std::list<halfedge_descriptor> interior_edges;
    std::set<halfedge_descriptor> included_map; 

    for(typename std::vector<face_descriptor>::const_iterator it = facets.begin(); it!= facets.end(); ++it) {
      Halfedge_around_face_circulator<PolygonMesh>  circ(halfedge(*it,pmesh),pmesh), done(circ);
      do {
        halfedge_descriptor h = *circ;        
        if(border_edges.find(h) == border_edges.end()){
          // do not remove included_map and use if(&*h < &*oh) { interior_edges.push_back(h) } 
          // which will change the order of edges from run to run
          halfedge_descriptor oh = opposite(h,pmesh);
          halfedge_descriptor h_rep = (h < oh) ? h : oh; // AF: was &*h < &*oh
          if(included_map.insert(h_rep).second) {
            interior_edges.push_back(h_rep);
          }
        }
      } while(++circ != done);
    }

    CGAL_TRACE_STREAM << "Test " << interior_edges.size() << " edges " << std::endl;
    //do not just use std::set (included_map) for iteration, the order effects the output (we like to make it deterministic)
    for(typename std::list<halfedge_descriptor>::iterator it = interior_edges.begin(); it != interior_edges.end();++it) {
      if(relax(*it)) {
        ++flips;
      }
    }

    CGAL_TRACE_STREAM << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

  double average_length(vertex_descriptor vh,
                        const std::set<face_descriptor>& interior_map, 
                        bool accept_internal_facets)
  {
    const Point_3& vp = ppmap[vh]; 
    Halfedge_around_target_circulator<PolygonMesh> circ(halfedge(vh,pmesh),pmesh), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      face_descriptor f(face(*circ,pmesh)), f_op(face(opposite(*circ,pmesh),pmesh));

      if(!accept_internal_facets) {
        if(interior_map.find(f) != interior_map.end() && interior_map.find(f_op) != interior_map.end())
        { continue; } // which means current edge is an interior edge and should not be included in scale attribute calculation
      }

      const Point_3& vq = ppmap[target(opposite(*circ,pmesh),pmesh)];
      sum += std::sqrt(CGAL::squared_distance(vp, vq));
      ++deg;
    } while(++circ != done);

    CGAL_assertion(deg != 0); // this might happen when accept_internal_facets = false but there is
    return sum/deg;
  }

  void calculate_scale_attribute(const std::vector<face_descriptor>& facets, 
                                 const std::set<face_descriptor>& interior_map,
                                 std::map<vertex_descriptor, double>& scale_attribute,
                                 bool accept_internal_facets) 
  {
    for(typename std::vector<face_descriptor>::const_iterator f_it = facets.begin(); f_it != facets.end(); ++f_it) {
      Halfedge_around_face_circulator<PolygonMesh> circ(halfedge(*f_it,pmesh),pmesh), done(circ);
      do {
        vertex_descriptor v = target(*circ,pmesh);
        std::pair<typename std::map<vertex_descriptor, double>::iterator, bool> v_insert 
          = scale_attribute.insert(std::make_pair(v, 0));
        if(!v_insert.second) { continue; } // already calculated
        v_insert.first->second = average_length(v, interior_map, accept_internal_facets);
      } while(++circ != done);
    }
  }

  bool contain_internal_facets(const std::vector<face_descriptor>& facets,
                               const std::set<face_descriptor>& interior_map) const
  {
    for(typename std::vector<face_descriptor>::const_iterator f_it = facets.begin(); f_it != facets.end(); ++f_it) {
      Halfedge_around_face_circulator<PolygonMesh> circ(halfedge(*f_it,pmesh),pmesh), done(circ);
      do {
        vertex_descriptor v = target(*circ,pmesh);
        Halfedge_around_target_circulator<PolygonMesh> circ_v(*circ,pmesh), done_v(circ_v);
        bool internal_v = true;
        do {
          face_descriptor f(face(*circ,pmesh)), f_op(face(opposite(*circ_v,pmesh),pmesh));

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
  Refine_Polyhedron_3(PolygonMesh& pmesh)
    : pmesh(pmesh), ppmap(get(vertex_point, pmesh))
  {}

  template<class FacetRange, class FacetOutputIterator, class VertexOutputIterator>
  void refine(FacetRange faces,
              FacetOutputIterator& facet_out,
              VertexOutputIterator& vertex_out,
              double alpha)
  {
    std::vector<face_descriptor> facets(boost::begin(faces), boost::end(faces));
      // do not use just std::set, the order effects the output (for the same input we want to get same output)
    std::set<face_descriptor> interior_map(facets.begin(), facets.end());

    // store boundary edges - to be used in relax 
    std::set<halfedge_descriptor> border_edges;
    for (typename std::vector<face_descriptor>::const_iterator it = facets.begin(); it != facets.end(); ++it){
      Halfedge_around_face_circulator<PolygonMesh>  circ(halfedge(*it,pmesh),pmesh), done(circ);
      do {
        if(interior_map.find(face(opposite(*circ,pmesh),pmesh)) == interior_map.end()) {
          border_edges.insert(*circ);
        }
      } while(++circ != done);
    }

    // check whether there is any need to accept internal facets
    bool accept_internal_facets = contain_internal_facets(facets, interior_map);
    std::map<vertex_descriptor, double> scale_attribute;
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

}; //end class Refine_Polyhedron_3

}//namespace internal

}//namespace Polygon_mesh_processing

}//namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_POLYHEDRON_3_H
