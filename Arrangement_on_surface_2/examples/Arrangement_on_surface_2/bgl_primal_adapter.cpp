//! \file examples/Arrangement_on_surface_2/bgl_primal_adapter.cpp
// Adapting an arrangement to a BGL graph.

#include <vector>

#include <CGAL/config.h>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <CGAL/graph_traits_Arrangement_2.h>
#include <CGAL/Arr_vertex_index_map.h>
#include <CGAL/property_map.h>

#include "arr_exact_construction_segments.h"
#include "Edge_length.h"

typedef CGAL::Arr_vertex_index_map<Arrangement>  Vertex_index_map;
typedef Edge_length<Arrangement>                 My_edge_length;

int main() {
  // Construct an arrangement of seven intersecting line segments.
  // We keep a handle for the vertex v0 that corresponds to the point (1,1).
  Point p1(1, 1), p2(1, 4), p3(2, 2), p4(3, 7), p5(4, 4), p6(7, 1), p7(9, 3);
  Arrangement arr;
  Segment s(p1, p6);
  Arrangement::Halfedge_handle e = insert_non_intersecting_curve(arr, s);
  Arrangement::Vertex_handle v0 = e->source();
  insert(arr, Segment(p1, p4));  insert(arr, Segment(p2, p6));
  insert(arr, Segment(p3, p7));  insert(arr, Segment(p3, p5));
  insert(arr, Segment(p6, p7));  insert(arr, Segment(p4, p7));

  // Create a mapping of the arrangement vertices to indices.
  Vertex_index_map index_map(arr);

  // Create a property map based on std::vector to keep the result distances.
  boost::vector_property_map<Number_type, Vertex_index_map>
    dist_map(static_cast<unsigned int>(arr.number_of_vertices()), index_map);

  // Perform Dijkstra's algorithm from the vertex v0.
  My_edge_length edge_length;
  boost::dijkstra_shortest_paths(arr, v0, boost::vertex_index_map(index_map).
                                 weight_map(edge_length).distance_map(dist_map).
                                 distance_zero(Number_type(0)).
                                 distance_inf(Number_type(1000)));

  // Print the distance of each vertex from v0.
  std::cout << "The graph distances of the arrangement vertices from ("
            << v0->point() << ") :\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    std::cout << "(" << vit->point() << ") at distance "
              << CGAL::to_double(dist_map[vit]) << std::endl;

  return 0;
}
