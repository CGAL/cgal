// wlop_simplify_and_regularize_test.cpp

//----------------------------------------------------------
// Test the wlop simplify and regularize method:
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// wlop_simplify_and_regularize_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/IO/read_xyz_points.h>

#include <deque>
#include <cstdlib>
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
void test_wlop_simplify_and_regularize(std::vector<Point>& points, // input point set                            
                              double retain_percentage, // percentage of points to remove
                              double neighbor_radius, // neighborhood size
                              unsigned int iter_number, // iteration number
                              bool need_compute_density)

{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Run WLOP simplify and regularize, (retain_percentage: "
            << retain_percentage << "%, neighbor_radius="
            << neighbor_radius << ")...\n";

  // Make room for sample points
  std::vector<Point> points_sampled;
  points_sampled.resize(points.size() * (retain_percentage / 100.));

  // Run algorithm 
  std::vector<Point>::const_iterator sample_points_begin =
    CGAL::wlop_simplify_and_regularize_point_set(
            points.begin(), 
            points.end(), 
            retain_percentage, 
            neighbor_radius,
            iter_number,
            need_compute_density);

  // Copy results to sample points
  std::copy(sample_points_begin,
            static_cast<std::vector<Point>::const_iterator>(points.end()),
            points_sampled.begin());

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
  std::cerr << "WLOP simplify and regularize" << std::endl;

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

  //Algorithm parameters
  const double retain_percentage = 10;   // percentage of points to retain.
  const double neighbor_radius = 0.25;   // neighbors size.
  const unsigned int iter_number = 30;     // number of iterations.
  const bool need_compute_density = true;  // if needed to compute density.

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
    std::vector<Point> points;
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

    test_wlop_simplify_and_regularize(points, 
                                      retain_percentage,
                                      neighbor_radius,
                                      iter_number,
                                      need_compute_density);

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

