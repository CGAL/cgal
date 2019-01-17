#ifndef TYPES_H_
#define TYPES_H_

#include "enriched_polyhedron.h"

#include <CGAL/Search_traits.h>                   // Kd tree
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/AABB_tree.h>                       // AABB tree
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Enriched_polyhedron<Kernel, CGAL::Enriched_item> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Facet_primitive;
typedef CGAL::AABB_traits<Kernel, Facet_primitive> Facet_traits;
typedef CGAL::AABB_tree<Facet_traits> Facet_tree;

// The following cannot be complied under debug mode
/*typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron> Edge_primitive;
typedef CGAL::AABB_traits<Kernel, Edge_primitive> Edge_traits;
typedef CGAL::AABB_tree<Edge_traits> Edge_tree;

struct Vertex_property_map {      // definition of the property map of a vertex
  typedef Point value_type;
  typedef const value_type &reference;
  typedef const Vertex_handle key_type;
  typedef boost::readable_property_map_tag category;
};
inline Vertex_property_map::reference get(Vertex_property_map, 
                                          Vertex_property_map::key_type p) {
  return p->point();
}
typedef CGAL::Search_traits_3<Kernel> Vertex_traits_base;
typedef CGAL::Search_traits_adapter<Vertex_handle, Vertex_property_map, Vertex_traits_base> Vertex_traits;
typedef CGAL::Orthogonal_k_neighbor_search<Vertex_traits> K_neighbor_search;
typedef K_neighbor_search::Tree Vertex_tree;
typedef K_neighbor_search::Distance Vertex_distance;*/

// the following make the Enriched_polyhedron compatible with Polyhedron_3 when using AABB
namespace boost {
  template<>
  struct graph_traits<Polyhedron> :
    public graph_traits<CGAL::Polyhedron_3<Kernel, CGAL::Enriched_item> >{};

  template<>
  struct graph_traits<Polyhedron const> :
    public graph_traits<CGAL::Polyhedron_3<Kernel, CGAL::Enriched_item> const>{};

  template<class Tag>
  struct property_map<Polyhedron, Tag> :
    public property_map<CGAL::Polyhedron_3<Kernel, CGAL::Enriched_item>, Tag>{};
}

namespace CGAL {
  template<>
  struct halfedge_graph_traits<Polyhedron> :
    public halfedge_graph_traits<CGAL::Polyhedron_3<Kernel, CGAL::Enriched_item> >{};
}

#endif // TYPES_H_