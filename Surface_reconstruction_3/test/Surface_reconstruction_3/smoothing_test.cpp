// smoothing_test.cpp

// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the smoothing methods
// Input files are .xyz
// No output
//----------------------------------------------------------
// smoothing_test points1.xyz points2.xyz...


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/smooth_jet_fitting_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

// STL stuff
#include <deque>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;


// ----------------------------------------------------------------------------
// Private functions
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

  std::cerr << "Smooth Point Set (knn="
            << nb_neighbors_smooth_jet_fitting << "%=" << nb_neighbors << ")...\n";

  std::deque<Point> output;
  CGAL::smooth_jet_fitting_3(points.begin(), points.end(),
                             std::back_inserter(output),
                             nb_neighbors);

  // mutating version of the same function
  CGAL::smooth_jet_fitting_3(points.begin(), points.end(),
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
  std::cerr << "No output" << std::endl;

  // decode parameters
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
    return EXIT_FAILURE;
  }

  // Smoothing options
  const double nb_neighbors_smooth_jet_fitting = 0.05 /* % */; // K-nearest neighbors (smooth points by Jet Fitting)
                                                               // LS: was 20

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for(int i=1; i<argc; i++)
  {
    std::cerr << std::endl;

    // Load point set
    std::deque<Point> points;
    std::cerr << "Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(argv[i],
                                             std::back_inserter(points),
                                             false /*skip normals*/))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;

      test_smooth_jet_fitting(points, nb_neighbors_smooth_jet_fitting);
    }
    else
    {
      std::cerr << "Error: cannot read file " << argv[i] << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
    }
  } // for each input file

  std::cerr << std::endl;

  // Return accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}
