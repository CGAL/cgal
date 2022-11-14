//! \file examples/Arrangement_on_surface_2/bgl_primal_adapter.cpp
// Adapting an arrangement to a BGL graph.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/graph_traits_Arrangement_2.h>
#include <CGAL/Arr_vertex_index_map.h>

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <CGAL/property_map.h>

typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;

// A functor used to compute the length of an edge.
class Edge_length_func
{
public:

  // Boost property type definitions:
  typedef boost::readable_property_map_tag        category;
  typedef double                                  value_type;
  typedef value_type                              reference;
  typedef Arrangement_2::Halfedge_handle          key_type;

  double operator()(Arrangement_2::Halfedge_handle e) const
  {
    const double     x1 = CGAL::to_double (e->source()->point().x());
    const double     y1 = CGAL::to_double (e->source()->point().y());
    const double     x2 = CGAL::to_double (e->target()->point().x());
    const double     y2 = CGAL::to_double (e->target()->point().y());
    const double     diff_x = x2 - x1;
    const double     diff_y = y2 - y1;

    return std::sqrt(diff_x*diff_x + diff_y*diff_y);
  }
};

double get(const Edge_length_func& edge_length, Arrangement_2::Halfedge_handle e)
{
  return edge_length(e);
}

int main()
{
  Arrangement_2   arr;

  // Construct an arrangement of seven intersecting line segments.
  // We keep a handle for the vertex v_0 that corresponds to the point (1,1).
  Arrangement_2::Halfedge_handle  e =
    insert_non_intersecting_curve (arr, Segment_2 (Point_2 (1, 1),
                                                   Point_2 (7, 1)));
  Arrangement_2::Vertex_handle    v0 = e->source();
  insert (arr, Segment_2 (Point_2 (1, 1), Point_2 (3, 7)));
  insert (arr, Segment_2 (Point_2 (1, 4), Point_2 (7, 1)));
  insert (arr, Segment_2 (Point_2 (2, 2), Point_2 (9, 3)));
  insert (arr, Segment_2 (Point_2 (2, 2), Point_2 (4, 4)));
  insert (arr, Segment_2 (Point_2 (7, 1), Point_2 (9, 3)));
  insert (arr, Segment_2 (Point_2 (3, 7), Point_2 (9, 3)));

  // Create a mapping of the arrangement vertices to indices.
  CGAL::Arr_vertex_index_map<Arrangement_2>        index_map(arr);

  // Perform Dijkstra's algorithm from the vertex v0.
  Edge_length_func                                      edge_length;

  boost::vector_property_map<double, CGAL::Arr_vertex_index_map<Arrangement_2> > dist_map(static_cast<unsigned int>(arr.number_of_vertices()), index_map);
  boost::dijkstra_shortest_paths(arr, v0,
                                 boost::vertex_index_map(index_map).
                                 weight_map(edge_length).
                                 distance_map(dist_map));

  // Print the results:
  Arrangement_2::Vertex_iterator      vit;

  std::cout << "The distances of the arrangement vertices from ("
            << v0->point() << ") :" << std::endl;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    std::cout << "(" << vit->point() << ") at distance "
              << dist_map[vit] << std::endl;

  return 0;
}
