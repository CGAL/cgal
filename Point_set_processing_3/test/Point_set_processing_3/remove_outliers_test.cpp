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

// This package
#include <CGAL/remove_outliers.h>
#include <CGAL/IO/read_xyz_points.h>

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
void test_avg_knn_sq_distance(std::deque<Point>& points, // input point set
                              unsigned int nb_neighbors_remove_outliers, // K-nearest neighbors
                              double removed_percentage) // percentage of points to remove
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Removes outliers wrt average squared distance to k nearest neighbors (remove "
            << removed_percentage << "%, k="
            << nb_neighbors_remove_outliers << ")...\n";

  // Removes outliers using erase-remove idiom
  points.erase(CGAL::remove_outliers(points.begin(), points.end(), 
                                     nb_neighbors_remove_outliers, removed_percentage),
               points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  std::deque<Point>(points).swap(points);


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

    // File name is:
    std::string input_filename  = argv[i];

    // Reads the point set file in points[].
    std::deque<Point> points;
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;

    // If XYZ file format:
    std::ifstream stream(input_filename.c_str());
    if(stream &&
       CGAL::read_xyz_points(stream, std::back_inserter(points)))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    //***************************************
    // Test
    //***************************************

    test_avg_knn_sq_distance(points, nb_neighbors_remove_outliers, removed_percentage);

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

