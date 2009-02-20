// average_spacing_example.cpp

//----------------------------------------------------------
// Average Spacing method.
// Read a point set and compute its average spacing.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// average_spacing_example point_set.xyz


// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// This package
#include <CGAL/average_spacing_3.h>
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

void test_average_spacing(std::deque<Point>& points, // input point set
                          unsigned int nb_neighbors_avg_spacing) // number of neighbors
{
  std::cerr << "Compute average spacing to knn (knn="<< nb_neighbors_avg_spacing << ")... ";

  typedef std::deque<Point>::iterator Iterator;
  std::cerr << CGAL::average_spacing_3<Iterator,FT>(points.begin(), points.end(),
                                                    nb_neighbors_avg_spacing) << "\n";
  
  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Average Spacing" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 1)
    {
      std::cerr << "Read a point set and compute its average spacing.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " point_set.xyz" << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
    }

    // Average Spacing options
    const unsigned int nb_neighbors_avg_spacing = 7; // K-nearest neighbors = 1 ring

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
      
    test_average_spacing(points,nb_neighbors_avg_spacing);

    return EXIT_SUCCESS;
}

