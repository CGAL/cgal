// Adapting the dual of an arrangement to a BGL graph.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>

#include <climits>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <CGAL/graph_traits_Dual_Arrangement_2.h>
#include <CGAL/Arr_face_index_map.h>

#include "arr_print.h"

typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2,
                                     unsigned int>      Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>             Arrangement_2;
typedef CGAL::Dual<Arrangement_2>                       Dual_arrangement_2;

// A BFS visitor class that associates each vertex with its discover time.
// In our case graph vertices represent arrangement faces.
template <class IndexMap>
class Discover_time_bfs_visitor : public boost::default_bfs_visitor
{
private:

  const IndexMap  *index_map;      // Mapping vertices to indices.
  unsigned int     time;           // The current time stamp.

public:

  // Constructor.
  Discover_time_bfs_visitor (const IndexMap& imap) :
    index_map (&imap),
    time (0)
  {}

  // Write the discover time for a given vertex.
  template <typename Vertex, typename Graph>
  void discover_vertex (Vertex u, const Graph& )
  {
    u->set_data (time);
    time++;
  }
};

int main ()
{
  Arrangement_2   arr;

  // Construct an arrangement of seven intersecting line segments.
  insert (arr, Segment_2 (Point_2 (1, 1), Point_2 (7, 1)));
  insert (arr, Segment_2 (Point_2 (1, 1), Point_2 (3, 7)));
  insert (arr, Segment_2 (Point_2 (1, 4), Point_2 (7, 1)));
  insert (arr, Segment_2 (Point_2 (2, 2), Point_2 (9, 3)));
  insert (arr, Segment_2 (Point_2 (2, 2), Point_2 (4, 4)));
  insert (arr, Segment_2 (Point_2 (7, 1), Point_2 (9, 3)));
  insert (arr, Segment_2 (Point_2 (3, 7), Point_2 (9, 3)));

  // Create a mapping of the arrangement faces to indices.
  CGAL::Arr_face_index_map<Arrangement_2>      index_map (arr);

  // Perform breadth-first search from the unbounded face, and use the BFS
  // visitor to associate each arrangement face with its discover time.
  Discover_time_bfs_visitor<CGAL::Arr_face_index_map<Arrangement_2> >
                                               bfs_visitor (index_map);
  Arrangement_2::Face_handle                   uf = arr.unbounded_face();

  boost::breadth_first_search (Dual_arrangement_2 (arr), uf,
                               boost::vertex_index_map (index_map).
                               visitor (bfs_visitor));

  // Print the results:
  Arrangement_2::Face_iterator      fit;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
  {
    std::cout << "Discover time " << fit->data() << " for ";
    if (fit != uf)
    {
      std::cout << "face ";
      print_ccb<Arrangement_2> (fit->outer_ccb());
    }
    else
      std::cout << "the unbounded face." << std::endl;
  }

  return (0);
}
