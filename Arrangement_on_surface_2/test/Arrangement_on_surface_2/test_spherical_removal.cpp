// #define CGAL_IDENTIFICATION_XY  2

#include <string>
#include <vector>

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Gmpq                                           Number_type;
typedef CGAL::Cartesian<Number_type>                         Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>    Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::Curve_2                               Curve_2;
typedef Geom_traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;

typedef Arrangement_2::Halfedge_handle                       Halfedge_handle;

bool test_one_file(std::ifstream& in_file, bool verbose)
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

  // Read the number of edges to remove.
  unsigned int num_edges_to_remove;
  in_file >> num_edges_to_remove;
  
  // Insert the curves incrementally.
  Arrangement_2 arr;
  std::vector<Halfedge_handle> halfedges;
  std::vector<std::pair<unsigned int, unsigned int> >::const_iterator cit;
  for (cit = curves.begin(); cit != curves.end(); ++cit) {
    X_monotone_curve_2 xcv(points[cit->first], points[cit->second]);
    std::cout << "inserting " << xcv << " ... ";
    std::cout.flush();
    Halfedge_handle he = CGAL::insert_non_intersecting_curve(arr, xcv);
    halfedges.push_back(he);
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
  
  // Remove the halfedges.
  std::vector<Halfedge_handle>::const_iterator hit;
  for (hit = halfedges.begin(); hit != halfedges.end(); ++hit) {
    if (num_edges_to_remove-- == 0) break;
    std::cout << "removing " << (*hit)->curve() << " ... ";
    std::cout.flush();
    CGAL::remove_edge(arr, *hit);
    std::cout << "removed" << std::endl;
  }
  halfedges.clear();
  
  // Verify that the arrangement is empty.
  unsigned int num_vertices = arr.number_of_vertices();
  unsigned int num_edges = arr.number_of_edges();
  unsigned int num_faces = arr.number_of_faces();
  arr.clear();
  
  if ((num_vertices != num_of_points) || num_edges || (num_faces != 1)) {
    std::cerr << " Failed. The arrangement is not empty:" << std::endl
              << "   V = " << arr.number_of_vertices()
              << ", E = " << arr.number_of_edges() 
              << ", F = " << arr.number_of_faces() 
              << std::endl;
    return false;
  }
  
  return true;
}

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    return -1;
  }
  bool verbose = false;
  if (argc > 2 && std::strncmp(argv[2], "-v", 2) == 0)
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
