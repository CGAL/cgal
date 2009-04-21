// normal_estimation.cpp

//----------------------------------------------------------
// Normal estimation:
// Read a point set or a mesh's set of vertices, compute and orient its normals,
// and save the point set.
// Input file formats are .off (mesh) and .xyz or .pwn (point set).
// Output file format is .xyz or .pwn (point set).
//----------------------------------------------------------
// normal_estimation file_in file_out [options]

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// This package
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>

#include "compute_normal.h"

#include <deque>
#include <cstdlib>
#include <fstream>


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
// Private functions
// ----------------------------------------------------------------------------

// Compute normals direction by Principal Component Analysis
void run_pca_estimate_normals(const PointList& points, // input point set
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

  std::cerr << "Estimate Normals Direction by PCA (k="
            << nb_neighbors_pca_normals << "%=" << nb_neighbors <<")...\n";

  CGAL::pca_estimate_normals(points.begin(), points.end(),
                             std::back_inserter(computed_normals),
                             nb_neighbors);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}

// Compute normals direction by Jet Fitting
void run_jet_estimate_normals(const PointList& points, // input point set
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

  std::cerr << "Estimate Normals Direction by Jet Fitting (k="
            << nb_neighbors_jet_fitting_normals << "%=" << nb_neighbors <<")...\n";

  CGAL::jet_estimate_normals(points.begin(), points.end(),
                             std::back_inserter(computed_normals),
                             nb_neighbors);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}

// Test Hoppe92 normal orientation using a Minimum Spanning Tree.
void run_mst_orient_normals(const PointList& points, // input point set
                            std::deque<Orientable_normal>& computed_normals, // normals to orient
                            unsigned int nb_neighbors_mst) // number of neighbors
{
  std::cerr << "Orient Normals with a Minimum Spanning Tree (k="<< nb_neighbors_mst << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  // Mark all normals as unoriented
  std::deque<Orientable_normal>::iterator n;
  for (n = computed_normals.begin(); n != computed_normals.end(); n++)
    n->set_orientation(false);

  // mst_orient_normals() requires an iterator over points
  // + property maps to access each point's index, position and normal.
  // We use the points index as iterator.
  boost::identity_property_map index_id; // identity
  CGAL::mst_orient_normals(
         (std::size_t)0, points.size(), // use the points index as iterator
         index_id, // index -> index property map = identity
         boost::make_iterator_property_map(points.begin(), index_id), // index -> position prop. map
         boost::make_iterator_property_map(computed_normals.begin(), index_id), // index -> normal prop. map
         nb_neighbors_mst);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Normal estimation" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 < 2)
    {
      std::cerr << "Read a point set or a mesh's set of vertices, compute and orient its normals,\n";
      std::cerr << "and save the point set.\n";
      std::cerr << "If the input mesh has normals, print the normals deviation.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .xyz or .pwn (point set).\n";
      std::cerr << "Options:\n";
      std::cerr << "  -estimate plane|quadric          Estimate normals direction\n";
      std::cerr << "  using a tangent plane or quadric (default=quadric)\n";
      std::cerr << "  -nb_neighbors_pca <int>          Number of neighbors\n";
      std::cerr << "  to compute tangent plane (default=0.15% of points)\n";
      std::cerr << "  -nb_neighbors_jet_fitting <int>  Number of neighbors\n";
      std::cerr << "  to compute quadric (default=default=0.1% of points)\n";
      std::cerr << "  -orient MST                      Orient normals\n";
      std::cerr << "  using a Minimum Spanning Tree (default=MST)\n";
      std::cerr << "  -nb_neighbors_mst <int>          Number of neighbors\n";
      std::cerr << "  to compute the MST (default=18)\n";
      return EXIT_FAILURE;
    }

    // Normals Computing options
    double nb_neighbors_pca_normals = 0.15 /* % */; // K-nearest neighbors (estimate normals by PCA)
    double nb_neighbors_jet_fitting_normals = 0.1 /* % */; // K-nearest neighbors (estimate normals by Jet Fitting)
    unsigned int nb_neighbors_mst = 18; // K-nearest neighbors = 3 rings (orient normals by MST)
    std::string estimate = "quadric"; // estimate normals by jet fitting
    std::string orient = "MST"; // orient normals using a Minimum Spanning Tree

    // decode parameters
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];
    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i])=="-estimate") {
        estimate = argv[++i];
        if (estimate != "plane" && estimate != "quadric")
          std::cerr << "invalid option " << argv[i] << "\n";
      }
      else if (std::string(argv[i])=="-nb_neighbors_pca") {
        nb_neighbors_pca_normals = atof(argv[++i]);
      }
      else if (std::string(argv[i])=="-nb_neighbors_jet_fitting") {
        nb_neighbors_jet_fitting_normals = atof(argv[++i]);
      }
      else if (std::string(argv[i])=="-orient") {
        orient = argv[++i];
        if (orient != "MST")
          std::cerr << "invalid option " << argv[i] << "\n";
      }
      else if (std::string(argv[i])=="-nb_neighbors_mst") {
        nb_neighbors_mst = atoi(argv[++i]);
      }
      else {
        std::cerr << "invalid option " << argv[i] << "\n";
      }
    }

    // Accumulated errors
    int accumulated_fatal_err = EXIT_SUCCESS;

    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Load mesh/point set
    //***************************************

    PointList points;

    // If OFF file format
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      // Read the mesh file in a polyhedron
      std::ifstream stream(input_filename.c_str());
      typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }

      // Convert Polyhedron vertices to point set.
      // Compute vertices' normals from connectivity.
      Polyhedron::Vertex_const_iterator v;
      for (v = input_mesh.vertices_begin(); v != input_mesh.vertices_end(); v++)
      {
        const Point& p = v->point();
        Vector n = compute_vertex_normal<Polyhedron::Vertex,Kernel>(*v);
        points.push_back(Point_with_normal(p,n));
      }
    }
    // If XYZ file format
    else if (extension == ".xyz" || extension == ".XYZ" ||
             extension == ".pwn" || extension == ".PWN")
    {
      // Read the point set file in points[]
      std::ifstream stream(input_filename.c_str());
      if(!stream || 
         !CGAL::read_xyz_point_set(stream,
                                   std::back_inserter(points)))
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Print status
    int nb_vertices = points.size();
    std::cerr << "Read file " << input_filename << ": " << nb_vertices << " vertices, "
                                                        << task_timer.time() << " seconds"
                                                        << std::endl;
    task_timer.reset();

    //***************************************
    // Check requirements
    //***************************************

    if (nb_vertices == 0)
    {
      std::cerr << "Error: empty file" << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Compute normals
    //***************************************

    std::deque<Orientable_normal> computed_normals;

    // Estimate normals direction.
    if (estimate == "plane")
      run_pca_estimate_normals(points, computed_normals, nb_neighbors_pca_normals);
    else if (estimate == "quadric")
      run_jet_estimate_normals(points, computed_normals, nb_neighbors_jet_fitting_normals);

    // Orient normals.
    if (orient == "MST")
      run_mst_orient_normals(points, computed_normals, nb_neighbors_mst);

    //***************************************
    // Save the point set
    //***************************************

    // Replace old normals by new ones
    PointList::iterator p;
    std::deque<Orientable_normal>::iterator n;
    for (p = points.begin(), n = computed_normals.begin(); p != points.end(); p++, n++)
      p->normal() = *n;

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    // If XYZ file format
    /*std::string*/ extension = output_filename.substr(output_filename.find_last_of('.'));
    if (extension == ".xyz" || extension == ".XYZ" ||
        extension == ".pwn" || extension == ".PWN")
    {
      std::ofstream stream(output_filename.c_str());
      if (!stream || 
          !CGAL::write_xyz_point_set(stream,
                                     points.begin(), points.end()) )
      {
        std::cerr << "Error: cannot write file " << output_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
        std::cerr << "Error: cannot write file " << output_filename << std::endl;
        return EXIT_FAILURE;
    }

    // Return accumulated fatal error
    std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
    return accumulated_fatal_err;
}

