// pca_normal_estimation_example.cpp

//----------------------------------------------------------
// Normal estimation:
// Read a point set, compute its normals using Principal Component Analysis, then orient them.
// Input format is .xyz.
// No output.
//----------------------------------------------------------
// pca_normal_estimation_example point_set.xyz

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/properties.h>

// This package
#include <CGAL/pca_normal_estimation.h>
#include <CGAL/mst_normal_orientation.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

#include <deque>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Orientable_normal_3<Kernel> Orientable_normal; // normal vector + orientation
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal; // position + normal vector
typedef std::deque<Point_with_normal> PointList;


// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

// Compute normals direction by Principal Component Analysis
void test_pca_normal_estimation(const PointList& points, // input point set
                                std::deque<Orientable_normal>& computed_normals, // normals to estimate
                                double nb_neighbors_pca_normals /* % */) // number of neighbors
{
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

  std::cerr << "ok" << std::endl;
}

// Test Hoppe92 normal orientation using a Minimum Spanning Tree.
void test_mst_normal_orientation(const PointList& points, // input point set
                                 std::deque<Orientable_normal>& computed_normals, // normals to orient
                                 unsigned int nb_neighbors_mst) // number of neighbors
{
  std::cerr << "Orient Normals with a Minimum Spanning Tree (knn="<< nb_neighbors_mst << ")...\n";

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

  std::cerr << "ok" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Normal estimation and orientation" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 1)
    {
      std::cerr << "Read a point set, compute its normals using Principal Component Analysis, then orient them.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " point_set.xyz" << std::endl;
      std::cerr << "Input file format is .xyz.\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
    }

    // Normals Computing options
    double nb_neighbors_pca_normals = 0.15 /* % */; // K-nearest neighbors (estimate normals by PCA)
    unsigned int nb_neighbors_mst = 18; // K-nearest neighbors = 3 rings (orient normals by MST)

    // File name is:
    std::string input_filename = argv[1];

    //***************************************
    // Load point set
    //***************************************

    PointList points;

    // Read the point set file in points[]
    std::cerr << "Open " << input_filename << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
                                             std::back_inserter(points),
                                             false /*skip normals*/))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Compute normals (PCA + MST)
    //***************************************

    std::deque<Orientable_normal> computed_normals;

    // Estimate normals direction.
    test_pca_normal_estimation(points, computed_normals, nb_neighbors_pca_normals);

    // Orient normals.
    test_mst_normal_orientation(points, computed_normals, nb_neighbors_mst);

    return EXIT_SUCCESS;
}

