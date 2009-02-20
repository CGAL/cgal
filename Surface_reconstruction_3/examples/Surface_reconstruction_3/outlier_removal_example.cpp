// outlier_removal_example.cpp

//----------------------------------------------------------
// Outlier Removal method.
// Read a point set and remove outliers.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// outlier_removal_example point_set.xyz


// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_avg_knn_sq_distance(std::deque<Point>& points, // input point set
                              double nb_neighbors_outlier_removal, // number of neighbors
                              double threshold_percent_avg_knn_sq_dst) // percentage of points to remove
{
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

  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Outlier Removal" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 1)
    {
      std::cerr << "Read a point set and remove outliers.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " point_set.xyz" << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
    }

    // Outlier Removal options
    const double threshold_percent_avg_knn_sq_dst = 5.0 /* % */; // percentage of outliers to remove
    const double nb_neighbors_outlier_removal = 0.05 /* % */; // K-nearest neighbors (outlier removal)

    // File name is:
    std::string input_filename = argv[1];

    //***************************************
    // Load point set
    //***************************************

    std::deque<Point> points;

    // Read the point set file in points[]
    std::cerr << "Open " << input_filename << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
                                             std::back_inserter(points), 
                                             false /*skip normals*/))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Test
    //***************************************
      
    test_avg_knn_sq_distance(points, nb_neighbors_outlier_removal, threshold_percent_avg_knn_sq_dst);

    return EXIT_SUCCESS;
}

