// normal_estimation_test.cpp

//----------------------------------------------------------
// Test the normal estimation methods:
// For each input point set, compute and orient its normals.
// If an input mesh has normals, print the normal deviation.
// Input file formats are .off, .xyz and .pwn.
// No output.
//----------------------------------------------------------
// normal_estimation_test points1.xyz points2.xyz...

// With iterator debugging this testsuite takes to long and the process gets killed
//#define _HAS_ITERATOR_DEBUGGING 0

// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

// This package
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>

#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<double> Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal; // position + normal vector
typedef std::vector<Point_with_normal> PointList;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

// Check the accuracy of normals direction estimation.
// If original normals are available, compare with them and count normals with large deviation.
// @return true on success.
bool verify_normal_direction(const PointList& points, // input points + computed normals
                             const std::vector<Vector>& original_normals) // may be empty
{
  bool success = true;

  bool points_have_original_normals = ! original_normals.empty();
  if (points_have_original_normals)
  {
    assert(points.size() == original_normals.size());

    std::cerr << "Compare with original normals:" << std::endl;

    double min_normal_deviation = DBL_MAX; // deviation / original normal
    double max_normal_deviation = DBL_MIN;
    double avg_normal_deviation = 0;
    int invalid_normals = 0; // #normals with large deviation
    PointList::const_iterator p;
    std::vector<Vector>::const_iterator n;
    for (p = points.begin(), n = original_normals.begin(); p != points.end(); p++, n++)
    {
      // Computes normal deviation.
      Vector v1 = *n; // original normal
      double norm1 = std::sqrt( v1*v1 );
      assert(norm1 != 0.0);
      Vector v2 = p->normal(); // computed normal
      double norm2 = std::sqrt( v2*v2 );
      assert(norm2 != 0.0);
      double cos_normal_deviation = (v1*v2)/(norm1*norm2);
      if (cos_normal_deviation < 0)
      {
        cos_normal_deviation = -cos_normal_deviation;
      }
      double normal_deviation = std::acos(cos_normal_deviation);

      // statistics about normals deviation
      min_normal_deviation = (std::min)(min_normal_deviation, normal_deviation);
      max_normal_deviation = (std::max)(max_normal_deviation, normal_deviation);
      avg_normal_deviation += normal_deviation;

      // count normal if large deviation
      bool valid = (normal_deviation <= CGAL_PI/3.); // valid if deviation <= 60 degrees
      if ( ! valid )
      {
        invalid_normals++;
      }
    }
    avg_normal_deviation /= double(points.size());

    std::cerr << "  Min normal deviation=" << min_normal_deviation*180.0/CGAL_PI << " degrees\n";
    std::cerr << "  Max normal deviation=" << max_normal_deviation*180.0/CGAL_PI << " degrees\n";
    std::cerr << "  Avg normal deviation=" << avg_normal_deviation*180.0/CGAL_PI << " degrees\n";
    if (invalid_normals > 0)
    {
      std::cerr << "  Error: " << invalid_normals << " normals have a deviation > 60 degrees\n";
      success = false;
    }
  }

  return success;
}

// Computes normals direction by Principal Component Analysis
// @return true on success.
bool run_pca_estimate_normals(PointList& points, // input points + output normals
                              unsigned int nb_neighbors_pca_normals, // number of neighbors
                              const std::vector<Vector>& original_normals) // may be empty
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Estimates Normals Direction by PCA (k="
            << nb_neighbors_pca_normals << ")...\n";

  CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
                             CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()), 
                             nb_neighbors_pca_normals);

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;

  // Check the accuracy of normals direction estimation.
  // If original normals are available, compare with them.
  return verify_normal_direction(points, original_normals);
}

// Computes normals direction by Jet Fitting
// @return true on success.
bool run_jet_estimate_normals(PointList& points, // input points + output normals
                              unsigned int nb_neighbors_jet_fitting_normals, // number of neighbors
                              const std::vector<Vector>& original_normals) // may be empty
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Estimates Normals Direction by Jet Fitting (k="
            << nb_neighbors_jet_fitting_normals << ")...\n";

  CGAL::jet_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
                             CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()), 
                             nb_neighbors_jet_fitting_normals);

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;

  // Check the accuracy of normals direction estimation.
  // If original normals are available, compare with them.
  return verify_normal_direction(points, original_normals);
}

// Check the accuracy of normal orientation.
// Count non-oriented normals.
// If original normals are available, compare with them and count flipped normals.
bool verify_normal_orientation(const PointList& points, // input points + computed normals
                               PointList::const_iterator unoriented_points_begin, // first pt w/ unoriented normal
                               const std::vector<Vector>& original_normals) // may be empty
{
  bool success = true;

  // Count non-oriented normals
  int unoriented_normals = 0;
  for (PointList::const_iterator p = unoriented_points_begin ; p != points.end(); p++)
  {
      unoriented_normals++;
  }
  if (unoriented_normals > 0)
  {
    std::cerr << "Error: " << unoriented_normals << " normals are unoriented\n";
    success = false;
  }

  // Compare oriented normals with original ones and count flipped normals
  bool points_have_original_normals = ! original_normals.empty();
  if (points_have_original_normals)
  {
    assert(points.size() == original_normals.size());

    std::cerr << "Compare with original normals:" << std::endl;

    int flipped_normals = 0; // #normals with wrong orientation
    PointList::const_iterator p;
    std::vector<Vector>::const_iterator n;
    for (p = points.begin(), n = original_normals.begin(); p != unoriented_points_begin; p++, n++)
    {
      Vector v1 = *n; // original normal
      double norm1 = std::sqrt( v1*v1 );
      assert(norm1 != 0.0);
      Vector v2 = p->normal(); // computed normal
      double norm2 = std::sqrt( v2*v2 );
      assert(norm2 != 0.0);
      double cos_normal_deviation = (v1*v2)/(norm1*norm2);
      if (cos_normal_deviation < 0) // if flipped
      {
        flipped_normals++;
      }
    }

    if (flipped_normals == 0)
      std::cerr << "  ok\n";
    else
      std::cerr << "  Error: " << flipped_normals << " normal(s) are flipped\n";
  }

  return success;
}

// Hoppe92 normal orientation using a Minimum Spanning Tree.
// @return true on success.
bool run_mst_orient_normals(PointList& points, // input points + input/output normals
                            unsigned int nb_neighbors_mst, // number of neighbors
                            const std::vector<Vector>& original_normals) // may be empty
{
#if (BOOST_VERSION / 100) == 1054
  std::cerr <<
    "In run_mst_orient_normals():\n"
    "NOTICE: This function is incompatible with Boost 1.54, "
    "and will not be tested. See the following bug:\n"
    "  https://svn.boost.org/trac/boost/ticket/9012\n";
  return true;
#endif // Boost version is 1.54

  std::cerr << "Orients Normals with a Minimum Spanning Tree (k="<< nb_neighbors_mst << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  PointList::iterator unoriented_points_begin = 
    CGAL::mst_orient_normals(points.begin(), points.end(),
    CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()), 
                             nb_neighbors_mst);

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;

  // Note: we do *not* delete points with unoriented normals in this test.
  // Instead, we check the accuracy of normal orientation and,
  // if original normals are available, compare with them.
  return verify_normal_orientation(points, unoriented_points_begin, original_normals);
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
      std::cerr << "If an input mesh has normals, print the normals deviation.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz..." << std::endl;
      std::cerr << "Input file formats are .off, .xyz and .pwn.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  // Normals Computing options
  unsigned int nb_neighbors_pca_normals = 18; // K-nearest neighbors = 3 rings (estimate normals by PCA)
  unsigned int nb_neighbors_jet_fitting_normals = 18; // K-nearest neighbors (estimate normals by Jet Fitting)
  unsigned int nb_neighbors_mst = 18; // K-nearest neighbors (orient normals by MST)

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
    PointList points;
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;

    // If OFF file format
    bool success = false;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      std::ifstream stream(input_filename.c_str());
      success = stream && 
                CGAL::read_off_points_and_normals(stream,
                                                  std::back_inserter(points),
                                                  CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()) 
                                                  );
    }
    // If XYZ file format
    else if (extension == ".xyz" || extension == ".XYZ" ||
             extension == ".pwn" || extension == ".PWN")
    {
      std::ifstream stream(input_filename.c_str());
      success = stream && 
                CGAL::read_xyz_points_and_normals(stream,
                                                  std::back_inserter(points),
                                                  CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())
                                                  );
    }
    if (success)
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
    // Check requirements
    //***************************************

    if (points.size() == 0)
    {
      std::cerr << "Error: empty file" << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    //***************************************
    // Copy original normals
    //***************************************

    std::vector<Vector> original_normals; 
    bool points_have_original_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( points_have_original_normals )
    {
      for (PointList::iterator p = points.begin() ; p != points.end(); p++)
        original_normals.push_back(p->normal());
    }

    //***************************************
    // Computes normals (PCA + MST)
    //***************************************

    // Estimates normals direction.
    success = run_pca_estimate_normals(points, nb_neighbors_pca_normals, original_normals);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

    // Orients normals.
    success = run_mst_orient_normals(points, nb_neighbors_mst, original_normals);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

    //***************************************
    // Computes normals (jet fitting + MST)
    //***************************************

    // Estimates normals direction
    success = run_jet_estimate_normals(points, nb_neighbors_jet_fitting_normals, original_normals);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

    // Orients normals
    success = run_mst_orient_normals(points, nb_neighbors_mst, original_normals);
    if ( ! success )
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}
