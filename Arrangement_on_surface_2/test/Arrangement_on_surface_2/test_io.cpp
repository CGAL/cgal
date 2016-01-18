#include <fstream>
#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef Segment_traits_2::Point_2                       Point_2;
typedef Segment_traits_2::X_monotone_curve_2            Segment_2;
typedef CGAL::Arrangement_2<Segment_traits_2>           Segment_arrangement_2;

typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Polyline_traits_2;
typedef Polyline_traits_2::Curve_2                      Polyline_2;
typedef CGAL::Arrangement_2<Polyline_traits_2>          Polyline_arrangement_2;

typedef CGAL::Arr_polyline_traits_2<Polyline_traits_2>  Polycurve_traits_2;
typedef Polycurve_traits_2::Curve_2                     Polycurve_2;
typedef CGAL::Arrangement_2<Polycurve_traits_2>         Polycurve_arrangement_2;

template <typename Arrangement>
bool test_io(Arrangement& arr, const char* filename = "arr.txt")
{
  std::cout << "Writing an arrangement of size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Write the arrangement to a file.
  std::ofstream out_file(filename);
  out_file << arr;
  out_file.close();

  // Read the arrangement from the file.
  Arrangement arr2;
  std::ifstream in_file(filename);
  in_file >> arr2;
  in_file.close();

  std::cout << "Read an arrangement of size:" << std::endl
            << "   V = " << arr2.number_of_vertices()
            << ",  E = " << arr2.number_of_edges()
            << ",  F = " << arr2.number_of_faces() << std::endl;

  return ((arr.number_of_vertices() == arr2.number_of_vertices()) &&
          (arr.number_of_halfedges() == arr2.number_of_halfedges()) &&
          (arr.number_of_faces() == arr2.number_of_faces()));
}

int main()
{
  // Segment arrangement
  std::list<Segment_2> segments1;
  segments1.push_back(Segment_2(Point_2(0, 0), Point_2(1, 0)));
  segments1.push_back(Segment_2(Point_2(1, 0), Point_2(2, 0)));
  segments1.push_back(Segment_2(Point_2(2, 0), Point_2(3, 0)));

  std::list<Segment_2> segments2;
  segments2.push_back(Segment_2(Point_2(3, 0), Point_2(2, 1)));
  segments2.push_back(Segment_2(Point_2(2, 1), Point_2(1, 2)));
  segments2.push_back(Segment_2(Point_2(1, 2), Point_2(0, 3)));

  std::list<Segment_2> segments3;
  segments3.push_back(Segment_2(Point_2(0, 3), Point_2(0, 2)));
  segments3.push_back(Segment_2(Point_2(0, 2), Point_2(0, 1)));
  segments3.push_back(Segment_2(Point_2(0, 1), Point_2(0, 0)));

  Segment_arrangement_2 segment_arr;
  CGAL::insert(segment_arr, segments1.begin(), segments1.end());
  CGAL::insert(segment_arr, segments2.begin(), segments2.end());
  CGAL::insert(segment_arr, segments3.begin(), segments3.end());
  if (! test_io(segment_arr, "arr_segments.dat")) {
    std::cerr << "IO for segments tailed!" << std::endl;
    return 1;
  }

  // Polyline arrangement
  Polyline_traits_2 polyline_traits;
  Polyline_traits_2::Construct_curve_2 polyline_ctr =
    polyline_traits.construct_curve_2_object();
  Polyline_2 polyline1 = polyline_ctr(segments1.begin(), segments1.end());
  Polyline_2 polyline2 = polyline_ctr(segments2.begin(), segments2.end());
  Polyline_2 polyline3 = polyline_ctr(segments3.begin(), segments3.end());
  std::list<Polyline_2> polylines;
  polylines.push_back(polyline1);
  polylines.push_back(polyline2);
  polylines.push_back(polyline3);

  Polyline_arrangement_2 polyline_arr(&polyline_traits);
  CGAL::insert(polyline_arr, polylines.begin(), polylines.end());
  if (! test_io(polyline_arr, "arr_polylines.dat")) return 1;

  // Poly-polyline arrangement
  Polycurve_traits_2 polycurve_traits;
  Polycurve_traits_2::Construct_curve_2 polycurve_ctr =
    polycurve_traits.construct_curve_2_object();
  Polycurve_2 polycurve = polycurve_ctr(polylines.begin(), polylines.end());
  Polycurve_arrangement_2 polycurve_arr(&polycurve_traits);
  CGAL::insert(polycurve_arr, polycurve);
  if (! test_io(polycurve_arr, "arr_polycurves.dat")) return 1;

  return 0;
}
