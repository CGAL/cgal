#ifndef CGAL_INTERSECTION_OF_PLANE_POLYHEDRA_3_USING_AABB_2_H
#define CGAL_INTERSECTION_OF_PLANE_POLYHEDRA_3_USING_AABB_2_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_segment_primitive.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>

#include <CGAL/intersection_of_Polyhedra_3.h>

#include <vector>
#include <map>
#include <deque>

#include <boost/tuple/tuple.hpp>
#include <boost/optional.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {

template<class Polyhedron, class PK>
class AABB_to_polyline_imp_2
{
private:
  typedef typename Polyhedron::Traits::Kernel K;
  typedef AABB_const_polyhedron_edge_primitive<K, Polyhedron>     AABB_primitive;
  typedef AABB_traits<K, AABB_primitive>                          AABB_traits_;
  typedef AABB_tree<AABB_traits_>                                  AABB_tree_;

  typedef typename AABB_tree_::Object_and_primitive_id             Object_and_primitive_id;
  typedef typename AABB_tree_::Primitive_id                        Primitive_id;

  typedef typename PK::Plane_3    Plane;
  typedef typename PK::Segment_3  Segment;
  typedef typename PK::Point_3    Point;

  typedef typename Polyhedron::Edge_const_iterator   Edge_const_iterator;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Polyhedron::Facet_const_handle    Facet_const_handle;
  typedef typename Polyhedron::Vertex_const_handle   Vertex_const_handle;
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
  // to unite halfedges under an "edge"
  struct Edge_comparator {
    bool operator()(Halfedge_const_handle h1, Halfedge_const_handle h2) const {
      return (std::min)(&*h1, &*(h1->opposite())) 
        < (std::min)(&*h2, &*(h2->opposite()));
    }
  };

  ////////////////////////////////////////////////
  // to represent a intersection point with plane
  struct Node {
    Point point;
  };

  struct Node_edge { 
    bool is_processed;
    Node_edge() : is_processed(false) 
    { }
  };

  // out_edges is setS for preventing parallel edges automatically
  typedef boost::adjacency_list<
    boost::setS, boost::vecS, boost::undirectedS,
    Node, Node_edge>
  Node_graph;

  typedef typename Node_graph::vertex_descriptor Vertex_g;
  typedef typename Node_graph::vertex_iterator   Vertex_iterator_g;
  typedef typename Node_graph::edge_descriptor   Edge_g;
  typedef typename Node_graph::out_edge_iterator Out_edge_iterator_g;
  ////////////////////////////////////////////////

  // Custom type to store intersections of an edge
  // there is two nodes (intersection points) when an edge is coplanar to plane
  struct Node_pair {
    Vertex_g v1, v2;
    int vertex_count;

    Node_pair() : vertex_count(0){ }
    Node_pair(Vertex_g v1, Vertex_g v2) : v1(v1), v2(v2), vertex_count(2) { }

    void put(Vertex_g v) {
      CGAL_assertion(vertex_count <= 2);

      (vertex_count == 0 ? v1 : v2) = v;
      ++vertex_count;
    }
  };

  typedef std::map<Halfedge_const_handle, Node_pair, Edge_comparator> Edge_node_map;
  typedef typename Edge_node_map::iterator Edge_node_map_iterator;

  enum Intersection_type { BOUNDARY, INTERVAL, PLANAR };

  typedef std::map<Halfedge_const_handle, Intersection_type, Edge_comparator> Edge_intersection_map;
  typedef typename Edge_intersection_map::iterator Edge_intersection_map_iterator;

  // member variables //
  AABB_tree_ tree;
  Node_graph node_graph;

  bool is_planar(Halfedge_const_handle hf, const Plane& plane) {
    const Point& s = hf->vertex()->point();
    const Point& t = hf->opposite()->vertex()->point();

    Oriented_side s_os = plane.oriented_side(s);
    Oriented_side t_os = plane.oriented_side(t);

    return t_os == ON_ORIENTED_BOUNDARY && s_os == ON_ORIENTED_BOUNDARY;
  }

  boost::tuple<Point, Intersection_type, Vertex_const_handle>
  halfedge_intersection(Halfedge_const_handle hf, const Plane& plane)
  {
    boost::tuple<Point, Intersection_type, Vertex_const_handle> ret;

    const Point& s = hf->vertex()->point();
    const Point& t = hf->opposite()->vertex()->point();

    Oriented_side s_os = plane.oriented_side(s);
    Oriented_side t_os = plane.oriented_side(t);

    // in case of planar just Intersection_type is not empty in tuple
    if(t_os == ON_ORIENTED_BOUNDARY && s_os == ON_ORIENTED_BOUNDARY) {
      ret.get<1>() = PLANAR;
      return ret;
    }
    // in case of boundary Intersection_type and vertex are not empty in tuple
    if(t_os == ON_ORIENTED_BOUNDARY || s_os == ON_ORIENTED_BOUNDARY) { 
      ret.get<1>() = BOUNDARY;
      ret.get<2>() = t_os == ON_ORIENTED_BOUNDARY ? hf->opposite()->vertex() :
                                                    hf->vertex();
      return ret;
    }

    CGAL_assertion(t_os != s_os); // should one positive one negative
    // in case of interval Intersection_type and point are not empty in tuple
    ret.get<1>() = INTERVAL;
    Object intersection = CGAL::intersection(plane, Segment(s,t));

    if(const Point* i_point  = object_cast<Point>(&intersection)) { 
      ret.get<0>() = *i_point;
    }
    else if(const Segment* i_segment = object_cast<Segment>(&intersection)) {
      // is it possible to predicate not-planar but construct segment ?
      CGAL_assertion(false); // not sure what to do in such case, will look back
      return Point(0,0,0);
    }
    else {
      // prediction indicates intersection but construction returns nothing, returning closest point
      CGAL_warning(false && "on interval case - no intersection found");
      ret.get<0>() = squared_distance(plane, s) < squared_distance(plane, t) ? s : t; 
    }
    return ret;
  }

  void add_one_ring(Vertex_const_handle v, Vertex_g v_g, Edge_node_map& edge_node_map) {
    Halfedge_around_vertex_const_circulator hf_c = v->vertex_begin();
    do {
      edge_node_map[hf_c].put(v_g);
    } 
    while(++hf_c != v->vertex_begin());
  }
  
  void add_edge_around_hf(Halfedge_const_handle hf, Edge_node_map& edge_node_map) {
    Edge_node_map_iterator hf_node_it = edge_node_map.find(hf);
    CGAL_assertion(hf_node_it != edge_node_map.end());

    Node_pair& hf_node_pair = hf_node_it->second;
    CGAL_assertion(hf_node_pair.vertex_count == 1);
    for(int i = 0; i < 2; ++i) // loop for hf and hf->opposite
    {
      if(i == 1) { hf = hf->opposite(); }
      if(hf->is_border_edge()) { continue;} 

      for(int i = 0; i < 2; ++i) 
      {
        Halfedge_const_handle facet_hf = i == 0 ? hf->next() : hf->prev();
        Edge_node_map_iterator facet_hf_it = edge_node_map.find(facet_hf);

        if(facet_hf_it != edge_node_map.end()) {
          Node_pair& facet_hf_node_pair = facet_hf_it->second;
          CGAL_assertion(facet_hf_node_pair.vertex_count == 1); // == 2 if hf is planar and that cant happen 
          boost::add_edge(hf_node_pair.v1, facet_hf_node_pair.v1, node_graph); 
        }
      }
    } // for hf and hf->opposite
  }
  
  void add_edge_around_vertex(Halfedge_const_handle hf, Edge_node_map& edge_node_map, bool both_end){
    Edge_node_map_iterator hf_node_it = edge_node_map.find(hf);
    CGAL_assertion(hf_node_it != edge_node_map.end());

    // going to be completed
    CGAL_assertion(false);
  }

  template<class OutputIterator>
  void intersect_plane(const Plane& plane, OutputIterator out) 
  {
    node_graph.clear();

    // find out intersecting halfedges (note that tree contains edges only with custom comparator)
    std::vector<Halfedge_const_handle> intersected_edges;
    tree.all_intersected_primitives(plane, std::back_inserter(intersected_edges));

    // create node graph from segments
    // each node is associated with multiple edges
    Edge_node_map edge_node_map;
    Edge_intersection_map edge_intersection_map;
    for(typename std::vector<Halfedge_const_handle>::iterator it = intersected_edges.begin();
      it != intersected_edges.end(); ++it)
    {
      Halfedge_const_handle hf = *it;
      Node_pair& assoc_nodes = edge_node_map[hf];
      CGAL_assertion(assoc_nodes.vertex_count < 3); // every Node_pair can at most contain 2 nodes

      if(assoc_nodes.vertex_count == 2) // means we associate nodes previously via other halfedges
      { 
        CGAL_assertion(is_planar(hf, plane));
        edge_intersection_map[hf] = PLANAR;
        continue; 
      } 
      if(assoc_nodes.vertex_count == 1) // there is a previous association
      {
        // we can construct one more node if halfedge is planar
        if(is_planar(hf, plane)) {
          edge_intersection_map[hf] = PLANAR;

          Vertex_g v = boost::add_vertex(node_graph);
          assoc_nodes.put(v);

          Vertex_const_handle v1 = hf->vertex();
          Vertex_const_handle v2 = hf->opposite()->vertex();

          const Point& included_loc = node_graph[assoc_nodes.v1].point;
          CGAL_assertion(included_loc == v1->point() || included_loc == v2->point());

          Vertex_const_handle not_included_v = included_loc == v1->point() ? v2 : v1;
          node_graph[v].point = not_included_v->point();

          // update edge_node_map for neighbor halfedges 
          add_one_ring(not_included_v, v, edge_node_map);
        }
      }
      // else case 
      if(assoc_nodes.vertex_count == 0) {
        boost::tuple<Point, Intersection_type, Vertex_const_handle> intersection 
          = halfedge_intersection(hf, plane);
        edge_intersection_map[hf] = intersection.get<1>();

        if(intersection.get<1>() == INTERVAL) 
        {
          Vertex_g v = boost::add_vertex(node_graph);
          assoc_nodes.put(v);
          node_graph[v].point = intersection.get<0>();
          // no related hf to edge node, done
        }
        else if(intersection.get<1>() == BOUNDARY) 
        {
          Vertex_g v = boost::add_vertex(node_graph);
          assoc_nodes.put(v);
          Vertex_const_handle v_h = intersection.get<2>();
          node_graph[v].point = v_h->point();
          
          // update edge_node_map for neighbor halfedges 
          add_one_ring(v_h, v, edge_node_map);
        }
        else 
        {
          CGAL_assertion(intersection.get<1>() == PLANAR);

          Vertex_g v1 = boost::add_vertex(node_graph);
          Vertex_g v2 = boost::add_vertex(node_graph);
          assoc_nodes.put(v1);
          assoc_nodes.put(v2);
          node_graph[v1].point = hf->vertex()->point();
          node_graph[v2].point = hf->opposite()->vertex()->point();

          // update edge_node_map for neighbor halfedges 
          add_one_ring(hf->vertex(), v1, edge_node_map);
          add_one_ring(hf->opposite()->vertex(), v2, edge_node_map);
        }
      } // if(assoc_nodes.vertex_count == 0)...
    } // for(typename std::vector<Halfedge_const_handle>::iterator it = intersected_edges.begin()...
    
    // introduce node connectivity
    std::cout << "intersected edges: " << intersected_edges.size() << std::endl;
    
    for(typename std::vector<Halfedge_const_handle>::iterator it = intersected_edges.begin();
      it != intersected_edges.end(); ++it)
    {
      Halfedge_const_handle hf = *it;
      Edge_intersection_map_iterator intersection_it = edge_intersection_map.find(hf);
      CGAL_assertion(intersection_it != edge_intersection_map.end());

      Intersection_type intersection = intersection_it->second;
      if(intersection == INTERVAL) {
        add_edge_around_hf(hf, edge_node_map);
      }
      else {
        // this part is missing currently, going to be added
        CGAL_assertion(false);
      }
      //else if(intersection == BOUNDARY) {
      //  //add_edge_around_vertex()
      //}
    }
    construct_polylines(out);
  }

  // find a new node to advance, if no proper node exists return v
  Vertex_g find_next(Vertex_g v) {
    Vertex_g v_f;

    Out_edge_iterator_g ei_b, ei_e;
    bool found = false;
    for(boost::tie(ei_b, ei_e) = boost::out_edges(v, node_graph); ei_b != ei_e; ++ei_b) {
      if(node_graph[*ei_b].is_processed) { continue; } // do not go over same edge twice
      found = true;
      CGAL_assertion(boost::source(*ei_b, node_graph) == v);
      v_f = boost::target(*ei_b, node_graph);
      node_graph[*ei_b].is_processed = true;
      break;
    }

    if(found) {
      CGAL_assertion(v_f != v);
      return v_f;
    }
    return v;
  }

  template<class OutputIterator>
  void construct_polylines(OutputIterator out) {
    std::cout << "size of vertices: " << boost::num_vertices(node_graph) << std::endl;
    std::cout << "size of edges: " << boost::num_edges(node_graph) << std::endl;

    Vertex_iterator_g v_b, v_e;
    for(boost::tie(v_b, v_e) = boost::vertices(node_graph); v_b != v_e; ++v_b)
    {
      Vertex_g v = *v_b;
      Out_edge_iterator_g e_b, e_e;
      boost::tie(e_b, e_e) = boost::out_edges(v, node_graph);
      // isolated intersection, construct a polyline with one point
      if(e_b == e_e) {
        std::vector<Point> polyline;
        polyline.push_back(node_graph[v].point);
        *out++ = polyline;
        continue;
      }

      Vertex_g v_head = v;
      // there are edges, construct one or more polylines
      for( ;e_b != e_e; ++e_b) // keep current vertex active while there is 
      {
        if(node_graph[*e_b].is_processed) { continue; } // we previously put it to polylines

        // construct new polyline
        std::vector<Point> polyline;

        bool loop_ok = true;
        while(true) {
          polyline.push_back(node_graph[v].point);
          Vertex_g v_new = find_next(v);

          if(v_new == v) { // intersection at boundary
            if(loop_ok) { // this is the first border, so go to the other direction
              loop_ok = false;
              v_new = find_next(v_head); // try to go other direction;
              if(v_new == v_head) { break; }
              std::reverse(polyline.begin(), polyline.end());
              continue;
            }
            else { break; } // completed non-loop polyline break
          }

          if(v == v_head) { // loop is completed
            polyline.push_back(node_graph[v].point); // put first point again
            break;
          }
        } // while( true )...
        *out++ = polyline;

      } // for( ;e_b != e_e;...
    } // for(boost::tie(v_b, v_e) = boost::vertices(node_graph)...
  }

public:
  template<class InputIterator, class OutputIterator>
  void operator()(Polyhedron& mesh, 
    InputIterator plane_begin,
    InputIterator plane_end,
    OutputIterator out) 
  {
    tree.rebuild(mesh.edges_begin(), mesh.edges_end());
    for(; plane_begin != plane_end; ++plane_begin) {
      intersect_plane(*plane_begin, out);
    }
  }
};


template<class PK, class Polyhedron, class OutputIterator>
void intersection_of_plane_Polyhedra_3_using_AABB_2(Polyhedron& mesh, 
  const Point_3<PK>& center, 
  const Vector_3<PK>& base1,
  const Vector_3<PK>& base2,
  OutputIterator out
  )
{
  Plane_3<PK> plane(center, center + base1, center + base2);
  intersection_of_plane_Polyhedra_3_using_AABB_2<PK>(mesh, plane, out);
}

template<class PK, class Polyhedron, class OutputIterator>
void intersection_of_plane_Polyhedra_3_using_AABB_2(Polyhedron& mesh, 
  const Plane_3<PK>& plane,
  OutputIterator out
  )
{
  intersection_of_plane_Polyhedra_3_using_AABB_2<PK>(mesh, &plane, &plane + 1, out);
}

template<class PK, class Polyhedron, class InputIterator, class OutputIterator>
void intersection_of_plane_Polyhedra_3_using_AABB_2(Polyhedron& mesh, 
  InputIterator plane_begin,
  InputIterator plane_end,
  OutputIterator out) 
{
  AABB_to_polyline_imp_2<Polyhedron, PK>()(mesh, plane_begin, plane_end, out);
}


}// end of namespace CGAL
#endif //CGAL_INTERSECTION_OF_PLANE_POLYHEDRA_3_USING_AABB_2_H

