// surface_reconstruction_read_write_xyz_example.cpp

//----------------------------------------------------------
// Read/write .xyz point set.
// Input format is .xyz.
// Output format is .xyz.
//----------------------------------------------------------
// surface_reconstruction_read_write_xyz_example in_point_set.xyz out_point_set.xyz

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This package
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>
#include <CGAL/IO/surface_reconstruction_write_xyz.h>

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
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::deque<Point_with_normal> PointList;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Read/write point set" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 2)
    {
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " in_point_set.xyz out_point_set.xyz" << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "Output file format is .xyz.\n";
      return EXIT_FAILURE;
    }

    // File names are:
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];

    //***************************************
    // Load point set
    //***************************************

    PointList points;

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
    // Save point set
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;
    if(CGAL::surface_reconstruction_write_xyz(output_filename.c_str(),
                                              points.begin(), points.end()))
    {
      std::cerr << "ok" << std::endl;
    }
    else
    {
        std::cerr << "Error: cannot write file " << output_filename << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

