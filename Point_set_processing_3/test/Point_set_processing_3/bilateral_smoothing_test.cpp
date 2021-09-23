// bilateral_smoothing_test.cpp

//----------------------------------------------------------
// Test the smoothing methods:
// For each input point set, smooth it.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// bilateral_smoothing_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/IO/read_points.h>

#include <deque>
#include <cstdlib>
#include <fstream>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

template<typename Concurrency_tag>
void test_bilateral_smoothing(std::deque<PointVectorPair>& points,// input point set
                             unsigned int nb_neighbors, // number of neighbors
                             double sharpness_sigma)
{
  CGAL::Real_timer task_timer; task_timer.start();

  for (int i = 0; i < 3; i++)
  {
      CGAL::bilateral_smooth_point_set <Concurrency_tag>(
        points, nb_neighbors,
        CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
        normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()).
        sharpness_angle (sharpness_sigma));
  }

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
  std::cerr << "Bilateral smoothing test" << std::endl;

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
  const unsigned int nb_neighbors = 50; // K-nearest neighbors
  const double sharpness_sigma = 25; // control sharpness

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
    std::deque<PointVectorPair> points;
    std::cerr << "Opening " << argv[i] << " for reading..." << std::endl;

    // If XYZ file format:
    if(CGAL::IO::read_points(argv[i],
                             std::back_inserter(points),
                             CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                              .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())))
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

#ifdef CGAL_LINKED_WITH_TBB
    std::deque<PointVectorPair> points2(points);
#endif // CGAL_LINKED_WITH_TBB

    test_bilateral_smoothing<CGAL::Sequential_tag>(
      points, nb_neighbors, sharpness_sigma);

#ifdef CGAL_LINKED_WITH_TBB
    test_bilateral_smoothing<CGAL::Parallel_tag>(
      points2, nb_neighbors, sharpness_sigma);
#endif // CGAL_LINKED_WITH_TBB

  } // for each input file

  std::cerr << std::endl;
  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

