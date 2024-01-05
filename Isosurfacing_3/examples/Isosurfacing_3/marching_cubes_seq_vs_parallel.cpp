#include <CGAL/Simple_cartesian.h>
#include <CGAL/Isosurfacing_3/Implicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Timer.h>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;
using Point_range = std::vector<Point>;
using Triangle_range = std::vector<std::vector<std::size_t> >;

// Sphere function = Euclidean distance function to the origin
auto sphere_function = [&](const Point& p) -> FT
{
	return std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
};

int main(int, char**)
{
  // box domain and spacing vector
  const CGAL::Bbox_3 bbox{ -1.0, -1.0, -1.0,  1.0, 1.0, 1.0 };
  const FT spacing = 0.01;
  const Vector vec_spacing(spacing, spacing, spacing);

  // create domain with sphere function
  auto domain = CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain<Kernel>
	  (bbox, vec_spacing, sphere_function);

  // points and triangles for the output indexed mesh
  Point_range points;
  Triangle_range triangles;

  // run marching cubes sequential
  std::cout << "marching cubes sequential...";
  const FT isovalue = 0.8;
  CGAL::Timer timer;
  timer.start();
  CGAL::Isosurfacing::marching_cubes<CGAL::Sequential_tag>(domain, isovalue, points, triangles);
  timer.stop();
  std::cout << "done (" << timer.time() << "s, " << triangles.size() << " triangles)" << std::endl;
  CGAL::IO::write_OFF("output-seq.off", points, triangles);

  // clear points and triangles
  points.clear();
  triangles.clear();

  // run marching cubes parallel
  std::cout << "marching cubes parallel...";
  timer.start();
  CGAL::Isosurfacing::marching_cubes<CGAL::Parallel_tag>(domain, isovalue, points, triangles);
  timer.stop();
  std::cout << "done (" << timer.time() << "s, " << triangles.size() << " triangles)" << std::endl;
  CGAL::IO::write_OFF("output-parallel.off", points, triangles);

  return EXIT_SUCCESS;
}
