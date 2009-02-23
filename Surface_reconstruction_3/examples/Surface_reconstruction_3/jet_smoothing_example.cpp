// jet_smoothing_example.cpp

//----------------------------------------------------------
// Smoothing method.
// Read a point set and smooth it.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// jet_smoothing_example point_set.xyz

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// This package
#include <CGAL/jet_smoothing_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

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

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_smooth_jet_fitting(std::deque<Point>& points, // input point set
                             double nb_neighbors_smooth_jet_fitting) // number of neighbors
{
  // percentage -> number of neighbors
  int nb_neighbors = int(double(points.size()) * nb_neighbors_smooth_jet_fitting / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > points.size()-1)
    nb_neighbors = points.size()-1;

  std::cerr << "Smooth Point Set (knn="
            << nb_neighbors_smooth_jet_fitting << "%=" << nb_neighbors << ")...\n";

  std::deque<Point> output;
  CGAL::jet_smoothing_3(points.begin(), points.end(),
                        std::back_inserter(output),
                        nb_neighbors);

  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Smoothing" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 1)
    {
      std::cerr << "Read a point set and smoth it.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " point_set.xyz" << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
    }

    // Smoothing options
    const double nb_neighbors_smooth_jet_fitting = 0.1 /* % */; // K-nearest neighbors (smooth points by Jet Fitting)

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

    test_smooth_jet_fitting(points, nb_neighbors_smooth_jet_fitting);

    return EXIT_SUCCESS;
}

