// edge_aware_upsample_test.cpp

//----------------------------------------------------------
// Test the edge aware up-sample test method:
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// edge_aware_upsample_test points1.xyz points2.xyz...

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/IO/read_points.h>

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
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

void test_edge_aware_upsample(std::vector<PointVectorPair>& points, // input point set
                              double sharpness_sigma, //control sharpness
                              double edge_sensitivity,// more points will up-sample on edge
                              double neighbor_radius,  // initial neighbors size.
                              unsigned int times_of_output_points)

{
  CGAL::Real_timer task_timer; task_timer.start();
  std::cerr << "Running edge aware up-sample, (sharpness_sigma: "
            << sharpness_sigma << "%, number_of_output_points="
            << points.size() * times_of_output_points << ")...\n";

   //Run algorithm
  CGAL::edge_aware_upsample_point_set<Concurrency_tag>(
    points, std::back_inserter(points),
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
    normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()).
    sharpness_angle(sharpness_sigma).
    edge_sensitivity(edge_sensitivity).
    neighbor_radius(neighbor_radius).
    number_of_output_points(points.size() * times_of_output_points));


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
  std::cerr << "Edge aware up-sample" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if(argc < 2)
  {
      std::cerr << "Upsample each input point set.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  //Algorithm parameters
  const double sharpness_sigma = 25;   //control sharpness of the result.
  const double edge_sensitivity = 0;    // more points will up-sample on edge.
  const double neighbor_radius = 0.2;      // initial neighbors size.
  const unsigned int times_of_output_points = 4;

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
    std::vector<PointVectorPair> points;
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

    test_edge_aware_upsample(points,
                             sharpness_sigma,
                             edge_sensitivity,
                             neighbor_radius,
                             times_of_output_points);

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

