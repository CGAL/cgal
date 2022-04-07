// smoothing_test.cpp

//----------------------------------------------------------
// Test the smoothing methods:
// For each input point set, smooth it.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// smoothing_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/IO/read_points.h>

#include <deque>
#include <string>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_smooth_jet_fitting(std::deque<Point>& points,// input point set
                             unsigned int nb_neighbors_smooth_jet_fitting) // number of neighbors
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Smoothes Point Set (k=" << nb_neighbors_smooth_jet_fitting <<  ")...\n";

  CGAL::jet_smooth_point_set<Concurrency_tag>
    (points, nb_neighbors_smooth_jet_fitting);

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
  std::cerr << "Smoothing test" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if(argc < 2)
  {
      std::cerr << "For each input point set, smooth it.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  // Smoothing options
  const unsigned int nb_neighbors_smooth_jet_fitting = 24; // K-nearest neighbors (smooth points by Jet Fitting)

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

    test_smooth_jet_fitting(points, nb_neighbors_smooth_jet_fitting);

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

