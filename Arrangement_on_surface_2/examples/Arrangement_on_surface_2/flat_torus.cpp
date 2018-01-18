#include <string>
#include <cstring>
#include <vector>

// #include <CGAL/Gmpq.h>
// #include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_flat_torus_traits_2.h>
#include <CGAL/Arr_toroidal_topology_traits_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_batched_point_location.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/IO/Arr_iostream.h>

// typedef CGAL::Gmpq                                           Number_type;
// typedef CGAL::Cartesian<Number_type>                         Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel    Kernel;
// TODO
typedef CGAL::Arr_flat_torus_traits_2<Kernel>                Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::Curve_2                               Curve_2;
typedef Geom_traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_toroidal_topology_traits_2<Geom_traits_2>  Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;

typedef Arrangement_2::Halfedge_handle                       Halfedge_handle;
typedef Arrangement_2::Halfedge_iterator                     Halfedge_iterator;
typedef Arrangement_2::Vertex_const_handle                   Vertex_const_handle;

bool test_one_file(std::ifstream& in_file, bool /* verbose */)
{
  Geom_traits_2 traits;

  size_t i;

  // Read the points:
  size_t num_of_points;
  in_file >> num_of_points;
  std::cout << "#pts: " << num_of_points << std::endl;
  std::vector<Point_2> points(num_of_points);
  for (i = 0; i < num_of_points; ++i) traits.read(in_file, points[i]);

  // Read the x-monotone curves:
  size_t num_of_xcurves;
  in_file >> num_of_xcurves;
  std::cout << "#xcvs: " << num_of_xcurves << std::endl;
  std::vector<X_monotone_curve_2> xcurves(num_of_xcurves);
  for (i = 0; i < num_of_xcurves; ++i) traits.read(in_file, xcurves[i]);

  // Read the isolated points.
  size_t num_of_isolated_points;
  in_file >> num_of_isolated_points;
  std::cout << "#ipts: " << num_of_isolated_points << std::endl;
  std::vector<unsigned int> isolated_points(num_of_isolated_points);
  for (i = 0; i < num_of_isolated_points; ++i)
    in_file >> isolated_points[i];

  // Read the number of edges to remove.
  size_t num_edges_to_remove;
  in_file >> num_edges_to_remove;
  std::cout << "#e2r: " << num_edges_to_remove << std::endl;

  // Read the number of cells left.
  Arrangement_2::Size num_vertices_left, num_edges_left, num_faces_left;
  in_file >> num_vertices_left >> num_edges_left >> num_faces_left;

  std::cout << "#v: " << num_vertices_left << std::endl;
  std::cout << "#e: " << num_edges_left << std::endl;
  std::cout << "#f: " << num_faces_left << std::endl;

  Arrangement_2 arr(&traits);
  std::vector<Halfedge_handle> halfedges;
  std::vector< X_monotone_curve_2 >::const_iterator xit;

#if 0
  // Insert the curves incrementally.
  for (xit = xcurves.begin(); xit != xcurves.end(); ++xit) {
    std::cout << "inserting " << *xit << " ... ";
    std::cout.flush();
    Halfedge_handle he = CGAL::insert_non_intersecting_curve(arr, *xit);
    halfedges.push_back(he);
    std::cout << "inserted" << std::endl;
  }
#else
  // Insert the curves aggregately.
  std::cout << "inserting " << " ... ";
  std::cout.flush();
  CGAL::insert(arr, xcurves.begin(), xcurves.end());
  std::cout << "inserted" << std::endl;
  for (Halfedge_iterator hit = arr.halfedges_begin(); hit != arr.halfedges_end(); hit++) {
    halfedges.push_back(hit);
  }
#endif

  xcurves.clear();

  // Insert the isolated points.
  if (isolated_points.size() != 0) {
    std::vector<unsigned int>::const_iterator pit;
    for (pit = isolated_points.begin(); pit != isolated_points.end(); ++pit) {
      Point_2 point(points[*pit]);
      CGAL::insert_point(arr, point);
    }
  }
  isolated_points.clear();
  points.clear();

  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  {
    std::cout << "Faces:" << std::endl;
    Arrangement_2::Face_const_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
      std::cout << "  Face: "
                << &(*fit)
                << std::endl;
      std::cout << "  Outer CCBs: "
                << std::distance(fit->outer_ccbs_begin(), fit->outer_ccbs_end())
                << std::endl;
      std::cout << "  Inner CCBs: "
                << std::distance(fit->inner_ccbs_begin(), fit->inner_ccbs_end())
                << std::endl;
      std::cout << std::endl;
    }
  }

  // test some overlay
  Arrangement_2 arr2, out;
  CGAL::overlay(arr, arr2, out);

  // test batched point location
  std::list<Point_2> qpoints;

  typedef CGAL::Arr_point_location_result<Arrangement_2>  Point_location_result;
  typedef std::pair<Point_2, Point_location_result::Type> Query_result;

  std::list<Query_result>  qresults;
  locate(arr, qpoints.begin(), qpoints.end(), std::back_inserter(qresults));

  // test vertical decomposition
  typedef std::pair< Vertex_const_handle, std::pair< CGAL::Object, CGAL::Object> > Decompose_result;

  std::list<Decompose_result> dresults;
  decompose(arr, std::back_inserter(dresults));

  // Remove the halfedges.
  if (num_edges_to_remove-- > 0) {
    // Remove the first halfedge inserted. Then, remove the rest.
    std::vector<Halfedge_handle>::const_iterator hit = halfedges.begin();
    std::cout << "removing (" << (*hit)->source()->point()
              << ") => (" << (*hit)->target()->point()
              << ") ... ";
    std::cout.flush();
    CGAL::remove_edge(arr, *hit);
    std::cout << "removed" << std::endl;
    halfedges.clear();

    while (num_edges_to_remove-- && (arr.number_of_edges() > 0)) {
      Halfedge_iterator hi = arr.halfedges_begin();
      std::cout << "removing (" << hi->source()->point()
                << ") => (" << hi->target()->point()
                << ") ... ";
      std::cout.flush();
      CGAL::remove_edge(arr, hi);
      std::cout << "removed" << std::endl;
    }
  }

  // Verify the resulting arrangement.
  Arrangement_2::Size num_vertices = arr.number_of_vertices();
  Arrangement_2::Size num_edges = arr.number_of_edges();
  Arrangement_2::Size num_faces = arr.number_of_faces();

  if ((num_vertices != num_vertices_left) ||
      (num_edges != num_edges_left) ||
      (num_faces != num_faces_left))
  {
    std::cerr << " Failed. The number of arrangement cells is incorrect:"
              << std::endl
              << "   V = " << arr.number_of_vertices()
              << ", E = " << arr.number_of_edges()
              << ", F = " << arr.number_of_faces()
              << std::endl;
    arr.clear();
    return false;
  }

  arr.clear();
  return true;
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    return -1;
  }
  bool verbose = false;
  if ((argc > 2) && (std::strncmp(argv[2], "-v", 2) == 0))
    verbose = true;
  int success = 0;
  for (int i = 1; i < argc; ++i) {
    std::string str(argv[i]);
    if (str.empty()) continue;

    std::string::iterator itr = str.end();
    --itr;
    while (itr != str.begin()) {
      std::string::iterator tmp = itr;
      --tmp;
      if (*itr == 't')  break;

      str.erase(itr);
      itr = tmp;
    }
    if (str.size() <= 1) continue;
    std::ifstream inp(str.c_str());
    if (!inp.is_open()) {
      std::cerr << "Failed to open " << str << std::endl;
      return -1;
    }
    if (! test_one_file(inp, verbose)) {
      inp.close();
      std::cerr << str << ": ERROR" << std::endl;
      success = -1;
    }
    else std::cout <<str << ": succeeded" << std::endl;
    inp.close();
  }

  return success;
}
