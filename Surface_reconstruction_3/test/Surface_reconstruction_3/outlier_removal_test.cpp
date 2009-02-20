// outlier_removal_test.cpp

//----------------------------------------------------------
// Test the outlier removal methods:
// For each input point set, remove outliers.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// outlier_removal_test points1.xyz points2.xyz...


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/boost/graph/properties.h>

// This package
#include <CGAL/outlier_removal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

// STL
#include <deque>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_avg_knn_sq_distance(std::deque<Point>& points, // input point set
                              double nb_neighbors_outlier_removal, // number of neighbors
                              double threshold_percent_avg_knn_sq_dst) // percentage of points to remove
{
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(points.size()) * nb_neighbors_outlier_removal / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > points.size()-1)
    nb_neighbors = points.size()-1;

  std::cerr << "Remove outliers wrt average squared distance to knn (remove "
            << threshold_percent_avg_knn_sq_dst << "%, knn="
            << nb_neighbors_outlier_removal << "%=" << nb_neighbors << ")...\n";

  std::deque<Point> output;
  CGAL::outlier_removal_3(
                  points.begin(), points.end(),
                  std::back_inserter(output),
                  nb_neighbors,
                  threshold_percent_avg_knn_sq_dst);

  // mutating version of the same function
  points.erase(CGAL::outlier_removal_3(
                  points.begin(), points.end(),
                  nb_neighbors,
                  threshold_percent_avg_knn_sq_dst),
               points.end());

  long memory = CGAL::Memory_sizer().virtual_size();
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
  const double threshold_percent_avg_knn_sq_dst = 5.0 /* % */; // percentage of outliers to remove
  const double nb_neighbors_outlier_removal = 0.05 /* % */; // K-nearest neighbors (outlier removal)

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for(int i=1; i<argc; i++)
  {
    std::cerr << std::endl;

    //***************************************
    // Load point set
    //***************************************

    // File name is:
    std::string input_filename  = argv[i];

    std::deque<Point> points;

    std::cerr << "Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
                                             std::back_inserter(points),
                                             false /*skip normals*/))
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
      
    test_avg_knn_sq_distance(points, nb_neighbors_outlier_removal, threshold_percent_avg_knn_sq_dst);

  } // for each input file

  std::cerr << std::endl;

  // Return accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

