// normal_estimation.cpp

//----------------------------------------------------------
// Normal estimation:
// Reads a point set, compute and orient its normals,
// and save the point set.
// Input file formats are .off, .xyz and .pwn.
// Output file formats are .xyz and .pwn.
//----------------------------------------------------------
// normal_estimation file_in file_out [options]

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>

// This package
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <utility> // defines std::pair
#include <vector>
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

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// Computes normals direction by Principal Component Analysis
void run_pca_estimate_normals(PointList& points, // input points + output normals
                              unsigned int nb_neighbors_pca_normals) // number of neighbors
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Estimates Normals Direction by PCA (k="
            << nb_neighbors_pca_normals << ")...\n";

  // Estimates normals direction.
  // Note: pca_estimate_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  CGAL::pca_estimate_normals(points.begin(), points.end(),
                             CGAL::First_of_pair_property_map<PointVectorPair>(),
                             CGAL::Second_of_pair_property_map<PointVectorPair>(),
                             nb_neighbors_pca_normals);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}

// Computes normals direction by Jet Fitting
void run_jet_estimate_normals(PointList& points, // input points + output normals
                              unsigned int nb_neighbors_jet_fitting_normals) // number of neighbors
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Estimates Normals Direction by Jet Fitting (k="
            << nb_neighbors_jet_fitting_normals << ")...\n";

  // Estimates normals direction.
  // Note: jet_estimate_normals() requires an iterator over points
  // + property maps to access each point's position and normal.
  CGAL::jet_estimate_normals(points.begin(), points.end(),
                             CGAL::First_of_pair_property_map<PointVectorPair>(),
                             CGAL::Second_of_pair_property_map<PointVectorPair>(),
                             nb_neighbors_jet_fitting_normals);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}

// Hoppe92 normal orientation using a Minimum Spanning Tree.
void run_mst_orient_normals(PointList& points, // input points + input/output normals
                            unsigned int nb_neighbors_mst) // number of neighbors
{
  std::cerr << "Orients Normals with a Minimum Spanning Tree (k="<< nb_neighbors_mst << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  // Orients normals.
  // Note: mst_orient_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  PointList::iterator unoriented_points_begin =
    CGAL::mst_orient_normals(points.begin(), points.end(),
                             CGAL::First_of_pair_property_map<PointVectorPair>(),
                             CGAL::Second_of_pair_property_map<PointVectorPair>(),
                             nb_neighbors_mst);

  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented normals.
  points.erase(unoriented_points_begin, points.end());

  // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
  PointList(points).swap(points);

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
      std::cerr << "Reads a point set, compute and orient its normals,\n";
      std::cerr << "and save the point set.\n";
      std::cerr << "If the input mesh has normals, print the normals deviation.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off, .xyz and .pwn.\n";
      std::cerr << "Output file formats are .xyz and .pwn.\n";
      std::cerr << "Options:\n";
      std::cerr << "  -estimate plane|quadric          Estimates normals direction\n";
      std::cerr << "  using a tangent plane or quadric (default=quadric)\n";
      std::cerr << "  -nb_neighbors_pca <int>          Number of neighbors\n";
      std::cerr << "  to compute tangent plane (default=18)\n";
      std::cerr << "  -nb_neighbors_jet_fitting <int>  Number of neighbors\n";
      std::cerr << "  to compute quadric (default=18)\n";
      std::cerr << "  -orient MST                      Orient normals\n";
      std::cerr << "  using a Minimum Spanning Tree (default=MST)\n";
      std::cerr << "  -nb_neighbors_mst <int>          Number of neighbors\n";
      std::cerr << "  to compute the MST (default=18)\n";
      return EXIT_FAILURE;
    }

    // Normals Computing options
    unsigned int nb_neighbors_pca_normals = 18; // K-nearest neighbors = 3 rings (estimate normals by PCA)
    unsigned int nb_neighbors_jet_fitting_normals = 18; // K-nearest neighbors (estimate normals by Jet Fitting)
    unsigned int nb_neighbors_mst = 18; // K-nearest neighbors (orient normals by MST)
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
        nb_neighbors_pca_normals = atoi(argv[++i]);
      }
      else if (std::string(argv[i])=="-nb_neighbors_jet_fitting") {
        nb_neighbors_jet_fitting_normals = atoi(argv[++i]);
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
    // Loads point set
    //***************************************

    // Reads a .off or .xyz point set file in points[].
    PointList points;
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;

    // If OFF file format
    bool success = false;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      std::ifstream stream(input_filename.c_str());
      success = stream &&
                CGAL::read_off_points(stream,
                                      std::back_inserter(points),
                                      CGAL::First_of_pair_property_map<PointVectorPair>());
    }
    // If XYZ file format
    else if (extension == ".xyz" || extension == ".XYZ" ||
             extension == ".pwn" || extension == ".PWN")
    {
      std::ifstream stream(input_filename.c_str());
      success = stream &&
                CGAL::read_xyz_points(stream,
                                      std::back_inserter(points),
                                      CGAL::First_of_pair_property_map<PointVectorPair>());
    }
    if (!success)
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Prints status
    int nb_points = points.size();
    std::cerr << "Reads file " << input_filename << ": " << nb_points << " points, "
                                                         << task_timer.time() << " seconds"
                                                         << std::endl;
    task_timer.reset();

    //***************************************
    // Check requirements
    //***************************************

    if (nb_points == 0)
    {
      std::cerr << "Error: empty file" << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Computes normals
    //***************************************

    // Estimates normals direction.
    if (estimate == "plane")
      run_pca_estimate_normals(points, nb_neighbors_pca_normals);
    else if (estimate == "quadric")
      run_jet_estimate_normals(points, nb_neighbors_jet_fitting_normals);

    // Orient normals.
    if (orient == "MST")
      run_mst_orient_normals(points, nb_neighbors_mst);

    //***************************************
    // Saves the point set
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    // If XYZ file format
    /*std::string*/ extension = output_filename.substr(output_filename.find_last_of('.'));
    if (extension == ".xyz" || extension == ".XYZ" ||
        extension == ".pwn" || extension == ".PWN")
    {
      std::ofstream stream(output_filename.c_str());
      if (!stream ||
          !CGAL::write_xyz_points_and_normals(stream,
                                              points.begin(), points.end(),
                                              CGAL::First_of_pair_property_map<PointVectorPair>(),
                                              CGAL::Second_of_pair_property_map<PointVectorPair>()))
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

    // Returns accumulated fatal error
    std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
    return accumulated_fatal_err;
}

