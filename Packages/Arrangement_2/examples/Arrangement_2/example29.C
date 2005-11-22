// file: examples/Arrangement_2/example22.C


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>

#include <CGAL/Arr_bgl_dual_adaptor.h>
#include <CGAL/Object_color_map.h>

#include <map>

typedef int                                             Number_type;
typedef CGAL::Simple_cartesian<Number_type>             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      Segment_2;

typedef CGAL::Arr_face_extended_dcel<Traits, boost::default_color_type>
                                                        Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arrangement_2;
typedef Arrangement_2::Vertex_handle                    Vertex_handle;
typedef Arrangement_2::Halfedge_handle                  Halfedge_handle;
typedef Arrangement_2::Face_handle                      Face_handle;
typedef Arrangement_2::Face_iterator                    Face_iterator;


// compare struct for Face_handle.
struct less_face_handle
{
  bool operator()(Face_iterator s1, Face_iterator s2) const
  {
    return reinterpret_cast<long>(&(*s1)) < reinterpret_cast<long>(&(*s2));
  }
};


int main ()
{
  Arrangement_2   arr;
 
  // create arrangment in the form of 3*3
  insert_curve(arr, Segment_2 (Point_2(0, 0), Point_2(0, 3)));
  insert_curve(arr, Segment_2 (Point_2(1, 0), Point_2(1, 3)));
  insert_curve(arr, Segment_2 (Point_2(2, 0), Point_2(2, 3)));
  insert_curve(arr, Segment_2 (Point_2(3, 0), Point_2(3, 3)));
  
  insert_curve(arr, Segment_2 (Point_2(0, 0), Point_2(3, 0)));
  insert_curve(arr, Segment_2 (Point_2(0, 1), Point_2(3, 1)));
  insert_curve(arr, Segment_2 (Point_2(0, 2), Point_2(3, 2)));
  insert_curve(arr, Segment_2 (Point_2(0, 3), Point_2(3, 3)));
  
  // create a map to number the faces
  typedef std::map<Face_handle, long, less_face_handle> Face_map;
  Face_map face_2_index_map;
  long i = 0;
  for (Face_iterator it = arr.faces_begin(); it != arr.faces_end(); ++it) {
    face_2_index_map.insert(std::make_pair(it, i));
    ++i;
  }

  // run the breadth first search.
  std::list<long> face_visiting_order;
  boost::associative_property_map<Face_map> numbering_map(face_2_index_map);
  CGAL::Object_color_map<Face_handle> my_map;
  
  boost::breadth_first_search(make_arr_bgl_dual_adaptor(arr),
                              arr.faces_begin(), 
                              boost::color_map(my_map).visitor
                              (boost::make_bfs_visitor
                               (write_property(numbering_map,
                                               std::back_inserter
                                               (face_visiting_order),
                                               boost::on_discover_vertex()))));

  // print out the order of the visit
  std::cout << "The order that the algorithm has visited the faces is:"
            << std::endl;
  std::copy(face_visiting_order.begin(), face_visiting_order.end(),
            std::ostream_iterator<long>(std::cout, " "));
  std::cout << std::endl;
  return 0;
}
