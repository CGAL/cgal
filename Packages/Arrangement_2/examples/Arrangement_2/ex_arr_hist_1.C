// file: examples/Arrangement_2/ex_arr_hist_1.C

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>

typedef CGAL::Gmpq                                    Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Segment_2;
typedef CGAL::Arrangement_with_history_2<Traits_2>    Arr_with_hist_2;
typedef Arr_with_hist_2::Curve_handle                 Curve_handle;

int main ()
{
  Arr_with_hist_2   arr;

  Segment_2         s1 (Point_2 (1, 3), Point_2 (4, 3));
  Curve_handle      c1 = arr.insert (s1);
  Segment_2         s2 (Point_2 (3, 1), Point_2 (3, 5));
  Curve_handle      c2 = arr.insert (s2);
  Segment_2         s3 (Point_2 (2, 3), Point_2 (5, 3));
  Curve_handle      c3 = arr.insert (s3);

  // Print the arrangement.
  Arr_with_hist_2::Vertex_iterator  vit;

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    std::cout << "(" << vit->point() << ")" << std::endl;

  // Print the arrangement edges.
  Arr_with_hist_2::Edge_iterator       eit;
  Arr_with_hist_2::Origin_iterator     oit;

  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
  {
    std::cout << "[" << eit->curve() << "]. Origin: ";
    for (oit = eit->curve().data_begin();
	 oit != eit->curve().data_end(); ++oit)
    {
      std::cout << " [" << **oit << "]";
    }
    std::cout << std::endl;
  }

  return (0);
}

