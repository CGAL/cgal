// normal_estimation_test.cpp

// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the normal estimation methods
// Input files are .xyz
// No output
//----------------------------------------------------------
// normal_estimation_test points1.xyz points2.xyz...


#include <CGAL/basic.h> // include basic.h before testing #defines

#ifdef CGAL_USE_LAPACK


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/properties.h>

// This package
#include <CGAL/estimate_normals_pca_3.h>
#include <CGAL/estimate_normals_jet_fitting_3.h>
#include <CGAL/orient_normals_minimum_spanning_tree_3.h>
#include <CGAL/Oriented_normal_3.h>
#include <CGAL/Vector_index_property_map.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

// STL stuff
#include <list>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iterator>

// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Oriented_normal_3<Kernel> Normal;

// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

void test_pca(const std::vector<Point>& points, // input point set
              std::vector<Normal>& normals, // computed normals
              unsigned int k) // number of neighbors
{
  std::cerr << "  Estimate normals using KNN and point-based PCA...";
  CGAL::estimate_normals_pca_3(points.begin(),points.end(),std::back_inserter(normals),k);
  std::cerr << "ok" << std::endl;
}

void test_jet_fitting(const std::vector<Point>& points, // input point set
              std::vector<Normal>& normals, // computed normals
              unsigned int k) // number of neighbors)
{
  std::cerr << "  Estimate normals using KNN and jet fitting...";
  CGAL::estimate_normals_jet_fitting_3(points.begin(),points.end(),std::back_inserter(normals),k);
  std::cerr << "ok" << std::endl;
}

void test_orient_normals_MST(
              const std::vector<Point>& points, // input point set
              std::vector<Normal>& normals, // normals to orient
              unsigned int k) // number of neighbors
{
  std::cerr << "  Orient normals using a minimum spanning tree...";

  // orient_normals_minimum_spanning_tree_3() requires an iterator over points
  // + property maps to access each point's index, position and normal.
  // We use the points index as iterator.
  boost::identity_property_map index_id; // identity
  CGAL::orient_normals_minimum_spanning_tree_3(
         (std::size_t)0, points.size(), // use the points index as iterator
         index_id, // index -> index prop. map = identity
         boost::make_iterator_property_map(points.begin(), index_id), // index -> position prop. map
         boost::make_iterator_property_map(normals.begin(), index_id), // index -> normal prop. map
         k);

  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Normal estimation test" << std::endl;
  std::cerr << "No output" << std::endl;

  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " file1.xyz file2.xyz ..." << std::endl;
    return EXIT_FAILURE;
  }

  // for each input file
  const unsigned int k = 10; // # neighbors
  for(int i=1; i<argc; i++)
  {
    std::vector<Point> points;
    std::cerr << "  Open " << argv[i] << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(argv[i], std::back_inserter(points)))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
      
      std::vector<Normal> normals_pca, normals_jet_fitting;
      test_pca(points, normals_pca, k);
      test_jet_fitting(points, normals_jet_fitting, k);
      test_orient_normals_MST(points, normals_jet_fitting, k);
    }
    else
      std::cerr << "  Unable to open file " << argv[i] << std::endl;
  }
  return EXIT_SUCCESS;
}


#else // CGAL_USE_LAPACK


#include <iostream>
#include <cstdlib>

// ----------------------------------------------------------------------------
// Empty main() if LAPACK is not installed
// ----------------------------------------------------------------------------

int main()
{
    std::cerr << "Skip test as LAPACK is not installed" << std::endl;
    return EXIT_SUCCESS;
}

#endif // CGAL_USE_LAPACK

