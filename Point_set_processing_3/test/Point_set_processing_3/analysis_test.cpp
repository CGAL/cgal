// analysis_test.cpp

//----------------------------------------------------------
// Test the analysis methods:
// For each input point set, compute the average spacing.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// analysis_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/compute_average_spacing.h>
#include <CGAL/IO/read_xyz_points.h>

#include <deque>
#include <cstdlib>
#include <fstream>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<float> Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

// Computes average spacing.
void test_average_spacing(std::deque<Point>& points, // input point set
                          unsigned int nb_neighbors) // number of neighbors
{
  std::cerr << "Computes average spacing to k nearest neighbors (k="<< nb_neighbors << ")... ";
  CGAL::Timer task_timer; task_timer.start();

  FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(points.begin(), points.end(), nb_neighbors);
  std::cout << average_spacing << std::endl;

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
  std::cerr << "Analysis test" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if(argc < 2)
  {
      std::cerr << "For each input point set, compute the average spacing.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  // Average Spacing options
  const unsigned int nb_neighbors = 6; // K-nearest neighbors = 1 ring (average spacing)

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

    test_average_spacing(points, nb_neighbors);

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

