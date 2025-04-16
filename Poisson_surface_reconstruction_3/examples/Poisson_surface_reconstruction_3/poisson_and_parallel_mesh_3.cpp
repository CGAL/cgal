#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Poisson_reconstruction_function.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <boost/iterator/transform_iterator.hpp>

#include <iostream>
#include <fstream>
#include <type_traits>
#include <vector>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

namespace params = CGAL::parameters;

template<typename Concurrency_tag, typename PointSet>
void poisson_reconstruction(const PointSet& points, const char* output)
{
  typedef CGAL::Labeled_mesh_domain_3<Kernel> Mesh_domain;
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  // Poisson options
  FT sm_angle = 20.0; // Min triangle angle in degrees.
  FT sm_radius = 1.; // Max triangle size w.r.t. point set average spacing.
  FT sm_distance = 0.25; // Surface Approximation error w.r.t. point set average spacing.

  CGAL::Real_timer time;
  time.start();

  CGAL::Real_timer total_time;
  total_time.start();

  // Creates implicit function from the read points using the default solver.

  // Note: this method requires an iterator over points
  // + property maps to access each point's position and normal.
  Poisson_reconstruction_function function(points.begin(), points.end(),
                                           Point_map(), Normal_map());

  // Computes the Poisson indicator function f()
  // at each vertex of the triangulation.
  if(!function.compute_implicit_function())
  {
    std::cerr << "compute_implicit_function() failed." << std::endl;
    return;
  }

  time.stop();
  std::cout << "compute_implicit_function() : " << time.time() << " seconds." << std::endl;
  time.reset();
  time.start();

  // Computes average spacing
  FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>
    (points, 6 /* knn = 1 ring */, params::point_map(Point_map()));

  time.stop();
  std::cout << "Average spacing : " << time.time() << " seconds." << std::endl;
  time.reset();
  time.start();

  // Gets one point inside the implicit surface
  // and computes implicit function bounding sphere radius.
  const Sphere bsphere = function.bounding_sphere();
  FT radius = std::sqrt(bsphere.squared_radius());

  // Defines the implicit surface: requires defining a
  // conservative bounding sphere centered at inner point.
  FT sm_sphere_radius = 5.0 * radius;
  FT sm_dichotomy_error = sm_distance * average_spacing / 1000.0; // Dichotomy error must be << sm_distance
  std::cout << "dichotomy error = " << sm_dichotomy_error << std::endl;
  std::cout << "sm_dichotomy_error / sm_sphere_radius = " << sm_dichotomy_error / sm_sphere_radius << std::endl;

  time.stop();
  std::cout << "Surface created in " << time.time() << " seconds." << std::endl;
  time.reset();
  time.start();

  // Defines surface mesh generation criteria
  Mesh_criteria criteria(params::facet_angle = sm_angle,
                         params::facet_size = sm_radius * average_spacing,
                         params::facet_distance = sm_distance * average_spacing);

  Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(function, bsphere,
    params::relative_error_bound(sm_dichotomy_error / sm_sphere_radius));

  // Generates surface mesh with manifold option
  std::cout << "Start meshing...";
  std::cout.flush();
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      params::no_exude()
                                             .no_perturb()
                                             .manifold_with_boundary());

  time.stop();
  std::cout << "\nTet mesh created in " << time.time() << " seconds." << std::endl;
  time.reset();
  time.start();

  const auto& tr = c3t3.triangulation();
  if(tr.number_of_vertices() == 0)
  {
    std::cerr << "Triangulation empty!" << std::endl;
    return;
  }

  // saves reconstructed surface mesh
  Polyhedron output_mesh;
  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, output_mesh);

  time.stop();
  std::cout << "Surface extracted in " << time.time() << " seconds." << std::endl;
  time.reset();
  time.start();

  total_time.stop();
  std::cout << "Total time : " << total_time.time() << " seconds." << std::endl;

  CGAL::IO::write_polygon_mesh(output, output_mesh, params::stream_precision(17));
  std::cout << "File written " << output << std::endl;
}

int main(int argc, char* argv[])
{
  const std::string file = (argc > 1) ? std::string(argv[1])
                                      : CGAL::data_file_path("points_3/kitten.xyz");

  // Reads the point set file in points[].
  // Note: read_points() requires an iterator over points
  // + property maps to access each point's position and normal.
  PointList points;
  if(!CGAL::IO::read_points(file, std::back_inserter(points),
                            params::point_map(Point_map())
                                   .normal_map(Normal_map())))
  {
    std::cerr << "Error: cannot read file input file!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "File " << file << " has been read, " << points.size() << " points." << std::endl;

  std::cout << "\n\n### Sequential mode ###" << std::endl;
  poisson_reconstruction<CGAL::Sequential_tag>(points, "out_sequential.off");

#ifdef CGAL_LINKED_WITH_TBB
  std::cout << "\n\n### Parallel mode ###" << std::endl;
  poisson_reconstruction<CGAL::Parallel_tag>(points, "out_parallel.off");
#endif

  return EXIT_SUCCESS;
}
