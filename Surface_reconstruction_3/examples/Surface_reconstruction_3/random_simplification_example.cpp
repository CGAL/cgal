// random_simplification_example.cpp

//----------------------------------------------------------
// Random Point Set Simplification method.
// Read a point set and simplify it.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// random_simplification_example point_set.xyz

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// This package
#include <CGAL/random_simplification_3.h>
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


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_random_simplification(std::deque<Point>& points, // input point set
                                double random_simplification_percentage) // percentage of points to remove
{
  std::cerr << "Random point cloud simplification (remove "
            << random_simplification_percentage << "%)...\n";

  std::deque<Point> output;
  CGAL::random_simplification_3(
                  points.begin(), points.end(),
                  std::back_inserter(output),
                  random_simplification_percentage);

  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Random Point Set Simplification" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 1)
    {
      std::cerr << "Read a point set and simplify it.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " point_set.xyz" << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
    }

    // Random Point Set Simplification options
    const double random_simplification_percentage = 50.0 /* % */; // percentage of outliers to remove

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

    test_random_simplification(points, random_simplification_percentage);

    return EXIT_SUCCESS;
}

