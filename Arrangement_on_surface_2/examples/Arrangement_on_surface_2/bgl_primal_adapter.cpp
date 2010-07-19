//! \file examples/Arrangement_on_surface_2/bgl_primal_adapter.cpp
// Adapting an arrangement to a BGL graph.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/graph_traits_Arrangement_2.h>
#include <CGAL/Arr_vertex_index_map.h>

#include <climits>
#include <boost/graph/dijkstra_shortest_paths.hpp>

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

double get(Edge_length_func edge_length, Arrangement_2::Halfedge_handle e)
{
  return edge_length(e);
}

/* The folowing is a workaround for a bug in the BGL upto and including version
 * 103400.
 *
 * Unfortunately some of the calls to the get() function below from the BGL
 * code are qualified with the boost namespace, while others are not. For The
 * qualified calls the compiler naturally looks for the definition of the
 * function in boost namespace. For the other calls it searches the CGAL
 * namespace according to ADL (Koenig Lookup), as the type of the 1st
 * parameter is in CGAL namespace.
 *
 * One way to get around it is to provide 2 similar functions that do the
 * same thing. One in CGAL namespace provided in CGAL/Arr_vertex_map.h, and
 * the other in boost namespace below. The signature of the latter is slightly
 * changed to avoid redefinition. The type of its 1st parameter is defined in
 * boost namespace, and is a simple derivation of the 1st parameter of the
 * CGAL::get() function.
 */

namespace boost {

template <typename Arrangement_2>
class Arr_vertex_index_map_boost :
    public CGAL::Arr_vertex_index_map<Arrangement_2>
{
 public:
  typedef CGAL::Arr_vertex_index_map<Arrangement_2>     Base;
  /*! Default constructor. */
  Arr_vertex_index_map_boost() : Base() {}

  /*! Constructor from CGAL index map. */
  Arr_vertex_index_map_boost(Base & other) :
    CGAL::Arr_vertex_index_map<Arrangement_2>(other)
  {}
};

/*!
 * Get the index property-map function. Provided so that boost is able to
 * access the Arr_vertex_index_map above.
 * \param index_map The index map.
 * \param v A vertex handle.
 * \return The vertex index.
 */
template<class Arrangement>
unsigned int
get(const boost::Arr_vertex_index_map_boost<Arrangement> & index_map,
    typename Arrangement::Vertex_handle v) 
{
  const CGAL::Arr_vertex_index_map<Arrangement> & index_map_tmp =
    static_cast<const CGAL::Arr_vertex_index_map<Arrangement> &>(index_map);
  return CGAL::get<Arrangement>(index_map_tmp, v);
}

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
  CGAL::Arr_vertex_index_map<Arrangement_2>        index_map_tmp(arr);
  boost::Arr_vertex_index_map_boost<Arrangement_2> index_map(index_map_tmp);
  
  // Perform Dijkstra's algorithm from the vertex v0.
  Edge_length_func                                      edge_length;
  CGAL::Arr_vertex_property_map<Arrangement_2, double>  dist_map(index_map);
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
