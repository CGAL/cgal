// normal_estimation.cpp

//----------------------------------------------------------
// Normal estimation:
// Reads a point set, compute and orient its normals,
// and save the point set.
// Input file formats are .off, .xyz and .pwn.
// Output file formats are .xyz and .pwn.
//----------------------------------------------------------
// normal_estimation file_in file_out [options]

// This package
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>

#include <utility> // defines std::pair
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

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

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

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
  CGAL::pca_estimate_normals<Concurrency_tag>(points,
     nb_neighbors_pca_normals,
     CGAL::parameters::point_map (CGAL::First_of_pair_property_map<PointVectorPair>()).
     normal_map (CGAL::Second_of_pair_property_map<PointVectorPair>()));

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
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
  CGAL::jet_estimate_normals<Concurrency_tag>
    (points,
     nb_neighbors_jet_fitting_normals,
     CGAL::parameters::point_map (CGAL::First_of_pair_property_map<PointVectorPair>()).
     normal_map (CGAL::Second_of_pair_property_map<PointVectorPair>()));


  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}

// Compute normals direction using the VCM
void run_vcm_estimate_normals(PointList &points, // input points + output normals
                              double R, // radius of the offset
                              double r) { // radius used during the convolution
    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Estimates Normals Direction using VCM (R="
        << R << " and r=" << r << ")...\n";

  // Estimates normals direction.
  // Note: vcm_estimate_normals() requires an iterator over points
  // + property maps to access each point's position and normal.
    CGAL::vcm_estimate_normals(points, R, r,
                               CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                               normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
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
    CGAL::mst_orient_normals(points,
                             nb_neighbors_mst,
                             CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                             normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));

  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented normals.
  points.erase(unoriented_points_begin, points.end());

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
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
      std::cerr << "  -estimate plane|quadric|vcm          Estimates normals direction\n";
      std::cerr << "  using a tangent plane or quadric or vcm (default=quadric)\n";
      std::cerr << "  -nb_neighbors_pca <int>              Number of neighbors\n";
      std::cerr << "  to compute tangent plane (default=18)\n";
      std::cerr << "  -nb_neighbors_jet_fitting <int>      Number of neighbors\n";
      std::cerr << "  to compute quadric (default=18)\n";
      std::cerr << "  -offset_radius_vcm <double>           Offset radius\n";
      std::cerr << "  to compute VCM (default=0.1)\n";
      std::cerr << "  -convolve_radius_vcm <double>         Convolve radius\n";
      std::cerr << "  to compute VCM (default=0)\n";
      std::cerr << "  -orient MST                          Orient normals\n";
      std::cerr << "  using a Minimum Spanning Tree (default=MST)\n";
      std::cerr << "  -nb_neighbors_mst <int>              Number of neighbors\n";
      std::cerr << "  to compute the MST (default=18)\n";
      std::cerr << "Running with " << argv[0] << "data/ChineseDragon-10kv.off ChineseDragon-10kv.pwn"
                                   << " -nb_neighbors_jet_fitting 10 -nb_neighbors_mst 10\n";
    }

    // Normals Computing options
    unsigned int nb_neighbors_pca_normals = 18; // K-nearest neighbors = 3 rings (estimate normals by PCA)
    unsigned int nb_neighbors_jet_fitting_normals = 18; // K-nearest neighbors (estimate normals by Jet Fitting)
    unsigned int nb_neighbors_mst = 18; // K-nearest neighbors (orient normals by MST)
    double offset_radius_vcm = 0.1; // Offset radius (estimate normals by VCM)
    double convolve_radius_vcm = 0; // Convolve radius (estimate normals by VCM)
    std::string estimate = "quadric"; // estimate normals by jet fitting
    std::string orient = "MST"; // orient normals using a Minimum Spanning Tree

    // decode parameters
    std::string input_filename  = argc == 1 ?  CGAL::data_file_path("meshes/ChineseDragon-10kv.off") : argv[1];
    std::string output_filename = argc == 1 ?  "ChineseDragon-10kv.pwn" : argv[2];
    if (argc==1)
    {
      nb_neighbors_jet_fitting_normals = 10;
      nb_neighbors_mst = 10;
    }

    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i])=="-estimate") {
        estimate = argv[++i];
        if (estimate != "plane" && estimate != "quadric" && estimate != "vcm")
          std::cerr << "invalid option " << argv[i] << "\n";
      }
      else if (std::string(argv[i])=="-nb_neighbors_pca") {
        nb_neighbors_pca_normals = atoi(argv[++i]);
      }
      else if (std::string(argv[i])=="-nb_neighbors_jet_fitting") {
        nb_neighbors_jet_fitting_normals = atoi(argv[++i]);
      }
      else if (std::string(argv[i])=="-offset_radius_vcm") {
          offset_radius_vcm = atof(argv[++i]);
      }
      else if (std::string(argv[i])=="-convolve_radius_vcm") {
          convolve_radius_vcm = atof(argv[++i]);
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

    // Reads a point set file in points[].
    PointList points;
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;

    if(!CGAL::IO::read_points(input_filename.c_str(), std::back_inserter(points),
                              CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())))
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Prints status
    std::size_t nb_points = points.size();
    std::cerr << "Reads file " << input_filename << ": " << nb_points << " points, "
                                                         << task_timer.time() << " seconds"
                                                         << std::endl;
    task_timer.reset();

    //***************************************
    // Check requirements
    //***************************************

    if(nb_points == 0)
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
    else if (estimate == "vcm")
      run_vcm_estimate_normals(points, offset_radius_vcm, convolve_radius_vcm);

    // Orient normals.
    if (orient == "MST")
      run_mst_orient_normals(points, nb_neighbors_mst);

    //***************************************
    // Saves the point set
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    if(!CGAL::IO::write_points(output_filename, points,
                               CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
                                                .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())
                                                .stream_precision(17)))
    {
      std::cerr << "Error: cannot write file " << output_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Returns accumulated fatal error
    std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
    return accumulated_fatal_err;
}

