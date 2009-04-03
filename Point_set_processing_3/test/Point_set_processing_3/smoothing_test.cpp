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
#include <CGAL/jet_smoothing_3.h>
#include <CGAL/IO/read_xyz_point_set.h>

#include <deque>
#include <cstdlib>
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


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_smooth_jet_fitting(std::deque<Point>& points,// input point set
                             double nb_neighbors_smooth_jet_fitting) // number of neighbors
{
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(points.size()) * nb_neighbors_smooth_jet_fitting / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > points.size()-1)
    nb_neighbors = points.size()-1;

  std::cerr << "Smooth Point Set (k="
            << nb_neighbors_smooth_jet_fitting << "%=" << nb_neighbors << ")...\n";

  std::deque<Point> output;
  CGAL::jet_smoothing_3(points.begin(), points.end(),
                        std::back_inserter(output),
                        nb_neighbors);

  // mutating version of the same function
  CGAL::jet_smoothing_3(points.begin(), points.end(),
                        nb_neighbors);

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
  const double nb_neighbors_smooth_jet_fitting = 0.1 /* % */; // K-nearest neighbors (smooth points by Jet Fitting)

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

    // Read the point set file in points[].
    std::deque<Point> points;
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;

    // If XYZ file format:
    std::ifstream stream(input_filename.c_str());
    if(stream &&
       CGAL::read_xyz_point_set(stream,
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

    test_smooth_jet_fitting(points, nb_neighbors_smooth_jet_fitting);

  } // for each input file

  std::cerr << std::endl;

  // Return accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

