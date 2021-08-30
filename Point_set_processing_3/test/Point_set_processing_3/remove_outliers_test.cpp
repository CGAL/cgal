// remove_outliers_test.cpp

//----------------------------------------------------------
// Test the outlier removal methods:
// For each input point set, remove outliers.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// remove_outliers_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/property_map.h>

// This package
#include <CGAL/remove_outliers.h>
#include <CGAL/IO/read_points.h>

#include <deque>
#include <string>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<float> Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

// Removes outliers
template <class PointContainer, class PointMap>
void test_avg_knn_sq_distance(PointContainer& points, // input point set
                              unsigned int nb_neighbors_remove_outliers, // K-nearest neighbors
                              double removed_percentage,
                              PointMap point_map) // percentage of points to remove
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Removes outliers wrt average squared distance to k nearest neighbors (remove "
            << removed_percentage << "%, k="
            << nb_neighbors_remove_outliers << ")...\n";

  // Removes outliers using erase-remove idiom
  points.erase(CGAL::remove_outliers<CGAL::Parallel_if_available_tag>
               (points, nb_neighbors_remove_outliers,
                CGAL::parameters::threshold_percent(removed_percentage).
                                  point_map(point_map)),
               points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  PointContainer(points).swap(points);


  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "ok: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Outlier removal test" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if(argc < 2)
  {
      std::cerr << "For each input point set, remove outliers.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  // Outlier Removal options
  const double removed_percentage = 5.0 /* % */; // percentage of outliers to remove
  const unsigned int nb_neighbors_remove_outliers = 24; // K-nearest neighbors

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for(int i=1; i<argc; i++)
  {
    std::cerr << std::endl;

    //***************************************
    // Loads point set
    //***************************************

    // Reads the point set file in points[].
    std::deque<Point> points;
    std::cerr << "Open " << argv[i] << " for reading..." << std::endl;

    // If XYZ file format:
    if(CGAL::IO::read_points(argv[i], std::back_inserter(points)))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
    }
    else
    {
      std::cerr << "Error: cannot read file " << argv[i] << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    //***************************************
    // Test
    //***************************************

    test_avg_knn_sq_distance(points, nb_neighbors_remove_outliers, removed_percentage,
                             CGAL::Identity_property_map<Point>());

    struct A{};
    std::vector< std::pair<Point, A> > points_bis;
    points_bis.reserve(points.size());
    for (const Point& p : points)
      points_bis.push_back( std::make_pair(p, A()) );
    test_avg_knn_sq_distance(points_bis, nb_neighbors_remove_outliers, removed_percentage,
                             CGAL::First_of_pair_property_map<std::pair<Point,A>>());
  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}
