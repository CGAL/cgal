#define CGAL_IDENTIFICATION_XY  2

#define CGAL_ARR_GEODESIC_ARC_ON_SPHERE_DETAILS 1

#include <string>
#include <vector>
#include <fstream>

#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Exact_rational                                 Number_type;
typedef CGAL::Cartesian<Number_type>                         Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>    Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::Curve_2                               Curve_2;
typedef Geom_traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;

typedef Arrangement_2::Vertex_iterator                       Vertex_iterator;
typedef Arrangement_2::Halfedge_handle                       Halfedge_handle;
typedef Arrangement_2::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;

bool test_one_file(std::ifstream& in_file, bool /* verbose */)
{
  unsigned int i;

  // Read the points:
  unsigned int num_of_points;
  in_file >> num_of_points;
  std::vector<Point_2> points(num_of_points);
  for (i = 0; i < num_of_points; ++i)
    in_file >> points[i];

  // Read the curves:
  unsigned int num_of_curves;
  in_file >> num_of_curves;
  std::vector<std::pair<unsigned int, unsigned int> > curves(num_of_curves);
  for (i = 0; i < num_of_curves; ++i) {
    unsigned int j, k;
    in_file >> j >> k;
    curves[i] = std::pair<unsigned int, unsigned int>(j, k);
  }

  // Read the isolated points.
  unsigned int num_of_isolated_points;
  in_file >> num_of_isolated_points;
  std::vector<unsigned int> isolated_points(num_of_isolated_points);
  for (i = 0; i < num_of_isolated_points; ++i)
    in_file >> isolated_points[i];

  // Read the number of vertices to remove.
  unsigned int num_vertices_to_remove;
  in_file >> num_vertices_to_remove;

  // Read the number of cells left.
  unsigned int num_vertices_left, num_edges_left, num_faces_left;
  in_file >> num_vertices_left >> num_edges_left >> num_faces_left;
  // std::cout << "Expected number of cells left:"
  //             << "V = " << num_vertices_left
  //             << ", E = " << num_edges_left
  //             << ", F = " << num_faces_left
  //             << std::endl;

  // Insert the curves incrementally.
  Arrangement_2 arr;
  std::vector<std::pair<unsigned int, unsigned int> >::const_iterator cit;
  for (cit = curves.begin(); cit != curves.end(); ++cit) {
    X_monotone_curve_2 xcv(points[cit->first], points[cit->second]);
    std::cout << "inserting " << xcv << " ... ";
    std::cout.flush();
    CGAL::insert_non_intersecting_curve(arr, xcv);
    std::cout << "inserted" << std::endl;
  }
  curves.clear();

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

  // Remove the vertices.

  const Geom_traits_2* traits = arr.geometry_traits();
  Vertex_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    if (num_vertices_to_remove-- == 0) break;
    if (vit->degree() != 2) continue;
    Halfedge_around_vertex_circulator hit = vit->incident_halfedges();
    if (!traits->are_mergeable_2_object()(hit->curve(), hit->next()->curve()))
      continue;
    X_monotone_curve_2 c;
    traits->merge_2_object()(hit->curve(), hit->next()->curve(), c);
    arr.merge_edge(hit, hit->next(), c);
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
