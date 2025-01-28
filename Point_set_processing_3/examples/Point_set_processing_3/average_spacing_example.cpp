#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/IO/read_points.h>

#include <vector>
#include <fstream>
#include <boost/tuple/tuple.hpp>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

// Data type := index, followed by the point, followed by three integers that
// define the Red Green Blue color of the point.
typedef boost::tuple<int, Point, int, int, int> IndexedPointWithColorTuple;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("points_3/sphere_20k.xyz");

  // Reads a file in points.
  // As the point is the second element of the tuple (that is with index 1)
  // we use a property map that accesses the 1st element of the tuple.

  std::vector<IndexedPointWithColorTuple> points;
  if (!CGAL::IO::read_points(fname, std::back_inserter(points),
                             CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<1, IndexedPointWithColorTuple>())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize index and RGB color fields in tuple.
  // As the index and RGB color are respectively the first and third-fifth elements
  // of the tuple we use a get function from the property map that accesses the 0
  // and 2-4th elements of the tuple.
  for(unsigned int i = 0; i < points.size(); i++)
  {
    points[i].get<0>() = i; // set index value of tuple to i

    points[i].get<2>() = 0; // set RGB color to black
    points[i].get<3>() = 0;
    points[i].get<4>() = 0;
  }

  // Computes average spacing.
  const unsigned int nb_neighbors = 6; // 1 ring
  FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
                         points, nb_neighbors,
                         CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<1,IndexedPointWithColorTuple>()));

  std::cout << "Average spacing: " << average_spacing << std::endl;

  return EXIT_SUCCESS;
}

