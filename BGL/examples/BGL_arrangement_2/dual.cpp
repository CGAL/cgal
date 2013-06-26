//! \file examples/Arrangement_on_surface_2/bgl_dual_adapter.cpp
// Adapting the dual of an arrangement to a BGL graph.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/graph_traits_Dual_Arrangement_2.h>
#include <CGAL/Arr_face_index_map.h>

#include <climits>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>

#include "arr_print.h"

// A property map that reads/writes the information to/from the extended 
// face.
template <typename Arrangement, class Type> class Extended_face_property_map {
public:
  typedef typename Arrangement::Face_handle       Face_handle;

  // Boost property type definitions.
  typedef boost::read_write_property_map_tag      category;
  typedef Type                                    value_type;
  typedef value_type&                             reference;
  typedef Face_handle                             key_type;

  // The get function is required by the property map concept.
  friend reference get(const Extended_face_property_map&, key_type key)
  { return key->data(); }

  // The put function is required by the property map concept.
  friend void put(const Extended_face_property_map&,
                  key_type key, value_type val)
  { key->set_data(val); }
};

typedef CGAL::Cartesian<Number_type>                         Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                   Traits_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, unsigned int> Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>                  Ex_arrangement;
typedef CGAL::Dual<Ex_arrangement>                           Dual_arrangement;
typedef CGAL::Arr_face_index_map<Ex_arrangement>             Face_index_map;
typedef Extended_face_property_map<Ex_arrangement,unsigned int>
                                                             Face_property_map;
typedef Kernel::Point_2                                      Point_2;
typedef Kernel::Segment_2                                    Segment_2;

int main()
{
  // Construct an arrangement of seven intersecting line segments.
  Point_2 p1(1, 1), p2(1, 4), p3(2, 2), p4(3, 7), p5(4, 4), p6(7, 1), p7(9, 3);
  Ex_arrangement  arr;
  insert(arr, Segment_2(p1, p6));
  insert(arr, Segment_2(p1, p4));  insert(arr, Segment_2(p2, p6));
  insert(arr, Segment_2(p3, p7));  insert(arr, Segment_2(p3, p5));
  insert(arr, Segment_2(p6, p7));  insert(arr, Segment_2(p4, p7));

  // Create a mapping of the arrangement faces to indices.
  Face_index_map  index_map(arr);

  // Perform breadth-first search from the unbounded face, using the event
  // visitor to associate each arrangement face with its discover time.
  unsigned int    time = 0;
  boost::breadth_first_search(Dual_arrangement(arr), arr.unbounded_face(),
                              boost::vertex_index_map(index_map).visitor
                              (boost::make_bfs_visitor
                               (stamp_times(Face_property_map(), time,
                                            boost::on_discover_vertex()))));

  // Print the discover time of each arrangement face.
  Ex_arrangement::Face_iterator  fit;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    std::cout << "Discover time " << fit->data() << " for ";
    if (fit != arr.unbounded_face()) {
      std::cout << "face ";
      print_ccb<Ex_arrangement>(fit->outer_ccb());
    }
    else std::cout << "the unbounded face." << std::endl;
  }
  return 0;
}
