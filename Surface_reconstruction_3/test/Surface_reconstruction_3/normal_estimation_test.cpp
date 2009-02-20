// normal_estimation_test.cpp

//----------------------------------------------------------
// Test the normal estimation methods:
// For each input point set, compute and orient its normals.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// normal_estimation_test points1.xyz points2.xyz...


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/boost/graph/properties.h>

// This package
#include <CGAL/pca_normal_estimation.h>
#include <CGAL/jet_normal_estimation.h>
#include <CGAL/mst_normal_orientation.h>
#include <CGAL/Orientable_normal_3.h>
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
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Orientable_normal_3<Kernel> Orientable_normal; // normal vector + orientation

typedef std::deque<Point> PointList;


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

// Compute normals direction by Principal Component Analysis
// @return true on success.
bool test_pca_normal_estimation(const PointList& points, // input point set
                                std::deque<Orientable_normal>& computed_normals, // normals to estimate
                                double nb_neighbors_pca_normals /* % */) // number of neighbors
{
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(points.size()) * nb_neighbors_pca_normals / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > points.size()-1)
    nb_neighbors = points.size()-1;

  std::cerr << "Estimate Normals Direction by PCA (knn="
            << nb_neighbors_pca_normals << "%=" << nb_neighbors <<")...\n";

  CGAL::pca_normal_estimation(points.begin(), points.end(),
                              std::back_inserter(computed_normals),
                              nb_neighbors);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;

  return true;
}

// Compute normals direction by Jet Fitting
// @return true on success.
bool test_jet_normal_estimation(const PointList& points, // input point set
                                std::deque<Orientable_normal>& computed_normals, // normals to estimate
                                double nb_neighbors_jet_fitting_normals /* % */) // number of neighbors
{
  CGAL::Timer task_timer; task_timer.start();

  // percentage -> number of neighbors
  int nb_neighbors = int(double(points.size()) * nb_neighbors_jet_fitting_normals / 100.0);
  if (nb_neighbors < 7)
    nb_neighbors = 7;
  if ((unsigned int)nb_neighbors > points.size()-1)
    nb_neighbors = points.size()-1;

  std::cerr << "Estimate Normals Direction by Jet Fitting (knn="
            << nb_neighbors_jet_fitting_normals << "%=" << nb_neighbors <<")...\n";

  CGAL::jet_normal_estimation(points.begin(), points.end(),
                              std::back_inserter(computed_normals),
                              nb_neighbors);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;

  return true;
}

// Count non-oriented normals.
bool verify_normal_orientation(const PointList& points, // input point set
                                const std::deque<Orientable_normal>& computed_normals) // oriented normals
{
  bool success = true;
  
  // Count non-oriented normals
  int unoriented_normals = 0;
  std::deque<Orientable_normal>::const_iterator n;
  for (n = computed_normals.begin(); n != computed_normals.end(); n++)
  {
    // Check unit vector
    Vector v = *n;
    double norm = std::sqrt( v*v );
    assert(norm > 0.99 || norm < 1.01);

    // Check orientation
    if ( ! n->is_oriented() )
    {
      unoriented_normals++;
    }
  }
  if (unoriented_normals > 0)
  {
    std::cerr << "Error: " << unoriented_normals << " normals are unoriented\n";
    success = false;
  }

  return success;
}

// Test Hoppe92 normal orientation using a Minimum Spanning Tree.
// @return true on success.
bool test_mst_normal_orientation(const PointList& points, // input point set
                                 std::deque<Orientable_normal>& computed_normals, // normals to orient
                                 unsigned int nb_neighbors_mst) // number of neighbors
{
  std::cerr << "Orient Normals with a Minimum Spanning Tree (knn="<< nb_neighbors_mst << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  // Mark all normals as unoriented
  std::deque<Orientable_normal>::iterator n;
  for (n = computed_normals.begin(); n != computed_normals.end(); n++)
    n->set_orientation(false);
  
  // mst_normal_orientation() requires an iterator over points
  // + property maps to access each point's index, position and normal.
  // We use the points index as iterator.
  boost::identity_property_map index_id; // identity
  CGAL::mst_normal_orientation(
         (std::size_t)0, points.size(), // use the points index as iterator
         index_id, // index -> index property map = identity
         boost::make_iterator_property_map(points.begin(), index_id), // index -> position prop. map
         boost::make_iterator_property_map(computed_normals.begin(), index_id), // index -> normal prop. map
         nb_neighbors_mst);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;

  // Check the accuracy of normal orientation.
  // If original normals are available, compare with them.
  return verify_normal_orientation(points, computed_normals);
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Normal estimation test" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if(argc < 2)
  {
      std::cerr << "For each input point set, compute and orient its normals.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  // Normals Computing options
  double nb_neighbors_pca_normals = 0.15 /* % */; // K-nearest neighbors (estimate normals by PCA)
  double nb_neighbors_jet_fitting_normals = 0.1 /* % */; // K-nearest neighbors (estimate normals by Jet Fitting)
  unsigned int nb_neighbors_mst = 18; // K-nearest neighbors = 3 rings (orient normals by MST)

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

    PointList points;

    std::cerr << "Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
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
    // Compute normals (PCA + MST)
    //***************************************

    std::deque<Orientable_normal> computed_normals;
    bool success = false;

    // Estimate normals direction.
    success = test_pca_normal_estimation(points, computed_normals, nb_neighbors_pca_normals);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

    // Orient normals.
    // Check normal orientation.
    success = test_mst_normal_orientation(points, computed_normals, nb_neighbors_mst);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

    //***************************************
    // Compute normals (jet fitting + MST)
    //***************************************

    computed_normals.clear();

    // Estimate normals direction
    success = test_jet_normal_estimation(points, computed_normals, nb_neighbors_jet_fitting_normals);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

    // Orient normals
    // Check normal orientation.
    success = test_mst_normal_orientation(points, computed_normals, nb_neighbors_mst);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

  } // for each input file

  std::cerr << std::endl;

  // Return accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

