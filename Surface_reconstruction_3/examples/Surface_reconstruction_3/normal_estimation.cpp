// normal_estimation.cpp


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Normal estimation:
// Read a point set or a mesh's set of vertices, compute and orient its normals,
// and save the point set.
// If the input mesh has normals, print the normals deviation.
//----------------------------------------------------------
// normal_estimation file_in file_out [options]


#include <CGAL/basic.h> // include basic.h before testing #defines

#ifdef CGAL_USE_LAPACK


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// This package
#include <CGAL/estimate_normals_pca_3.h>
#include <CGAL/estimate_normals_jet_fitting_3.h>
#include <CGAL/orient_normals_minimum_spanning_tree_3.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>
#include <CGAL/IO/surface_reconstruction_read_pwn.h>
#include <CGAL/IO/surface_reconstruction_write_xyz.h>
#include <CGAL/IO/surface_reconstruction_write_pwn.h>

// This test
#include "enriched_polyhedron.h"

// STL stuff
#include <deque>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <math.h>
#ifndef M_PI
  #define M_PI       3.14159265358979323846
#endif

// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Simple_cartesian<float> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Orientable_normal_3<Kernel> Orientable_normal; // normal vector + orientation
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal; // position + normal vector

typedef std::deque<Point_with_normal> PointList;

// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

void estimate_normals_pca(const PointList& points, // input point set
                          std::deque<Orientable_normal>& normals, // computed normals
                          unsigned int k) // number of neighbors
{
  std::cerr << "Estimate normals using KNN and point-based PCA (knn="<<k << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  CGAL::estimate_normals_pca_3(points.begin(),points.end(),std::back_inserter(normals),k);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "ok: " << task_timer.time() << " seconds, "
                      << (memory>>20) << " Mb allocated"
                      << std::endl;
}

void estimate_normals_jet_fitting(const PointList& points, // input point set
                                  std::deque<Orientable_normal>& normals, // computed normals
                                  unsigned int k) // number of neighbors)
{
  std::cerr << "Estimate normals using KNN and jet fitting (knn="<<k << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  CGAL::estimate_normals_jet_fitting_3(points.begin(),points.end(),std::back_inserter(normals),k);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "ok: " << task_timer.time() << " seconds, "
                      << (memory>>20) << " Mb allocated"
                      << std::endl;
}

void orient_normals_MST(const PointList& points, // input point set
                        std::deque<Orientable_normal>& normals, // normals to orient
                        unsigned int k) // number of neighbors
{
  std::cerr << "Orient normals using a minimum spanning tree (knn="<<k << ")...\n";
  CGAL::Timer task_timer; task_timer.start();

  // orient_normals_minimum_spanning_tree_3() requires an iterator over points
  // + property maps to access each point's index, position and normal.
  // We use the points index as iterator.
  boost::identity_property_map index_id; // identity
  CGAL::orient_normals_minimum_spanning_tree_3(
         (std::size_t)0, points.size(), // use the points index as iterator
         index_id, // index -> index property map = identity
         boost::make_iterator_property_map(points.begin(), index_id), // index -> position prop. map
         boost::make_iterator_property_map(normals.begin(), index_id), // index -> normal prop. map
         k);


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
    std::cerr << "Normal estimation" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc<3)
    {
      std::cerr << "Read a point set or a mesh's set of vertices, compute and orient its normals,\n";
      std::cerr << "and save the point set.\n";
      std::cerr << "If the input mesh has normals, print the normals deviation.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .xyz or .pwn (point set).\n";
      std::cerr << "Options:\n";
      std::cerr << "  -estimate plane|quadric     Estimate normals direction\n";
      std::cerr << "  using a tangent plane or quadric (default=quadric)\n";
      std::cerr << "  -estimate_neighbors <int>   Number of neighbors\n";
      std::cerr << "  to compute tangent plane or quadric (default=50)\n";
      std::cerr << "  -orient MST                 Orient normals\n";
      std::cerr << "  using a Minimum Spanning Tree (default=MST)\n";
      std::cerr << "  -orient_neighbors <int>   Number of neighbors\n";
      std::cerr << "  to compute the MST (default=10)\n";
      return EXIT_FAILURE;
    }

    // Default options
    std::string estimate = "quadric"; // estimate normals by jet fitting
    int estimate_neighbors = 50; // # neighbors
    std::string orient = "MST"; // orient normals using a Minimum Spanning Tree
    int orient_neighbors = 10; // # neighbors

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
      else if (std::string(argv[i])=="-estimate_neighbors") {
        estimate_neighbors = atoi(argv[++i]);
      } 
      else if (std::string(argv[i])=="-orient") {
        orient = argv[++i];
        if (orient != "MST")
          std::cerr << "invalid option " << argv[i] << "\n";
      } 
      else if (std::string(argv[i])=="-orient_neighbors") {
        orient_neighbors = atoi(argv[++i]);
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

    PointList pwns;

    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      // Read the mesh file in a polyhedron
      std::ifstream stream(input_filename.c_str());
      typedef Enriched_polyhedron<Kernel,Enriched_items> Polyhedron;
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }

      // Compute vertices' normals from connectivity
      input_mesh.compute_normals();

      // Convert vertices and normals to PointList
      Polyhedron::Vertex_iterator v;
      for (v = input_mesh.vertices_begin(); v != input_mesh.vertices_end(); v++)
      {
        const Point& p = v->point();
        const Vector& n = v->normal();
        pwns.push_back(Point_with_normal(p,n));
      }
    }
    else if (extension == ".xyz" || extension == ".XYZ")
    {
      // Read the point set file in pwns[]
      if(!CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
                                                std::back_inserter(pwns)))
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (extension == ".pwn" || extension == ".PWN")
    {
      // Read the point set file in pwns[]
      if(!CGAL::surface_reconstruction_read_pwn(input_filename.c_str(),
                                                std::back_inserter(pwns)))
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
    long memory = CGAL::Memory_sizer().virtual_size();
    int nb_vertices = pwns.size();
    std::cerr << "Read file " << input_filename << ": " << nb_vertices << " vertices, "
                                                        << task_timer.time() << " seconds, "
                                                        << (memory>>20) << " Mb allocated"
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
    
    // Estimate normals direction
    if (estimate == "plane")
      estimate_normals_pca(pwns, computed_normals, estimate_neighbors);
    else if (estimate == "quadric")
      estimate_normals_jet_fitting(pwns, computed_normals, estimate_neighbors);
      
    // Orient normals
    if (orient == "MST")
      orient_normals_MST(pwns, computed_normals, orient_neighbors);
    
    // Check computed normals
    int unoriented_normals = 0;
    std::deque<Orientable_normal>::iterator n;
    for (n = computed_normals.begin(); n != computed_normals.end(); n++)
    {
      // Check unit vector
      Vector v = *n;
      double norm = std::sqrt( v*v );
      assert(norm > 0.99 || norm < 1.01);
        
      // Check orientation
      if ( ! n->is_oriented() )
        unoriented_normals++;
    }
    if (unoriented_normals > 0)
    {
      std::cerr << "Error: " << unoriented_normals << " normal(s) are unoriented\n";
      accumulated_fatal_err = EXIT_FAILURE; // set error and continue
    }

    //***************************************
    // Compare with original normals
    //***************************************

    assert(pwns.begin() != pwns.end());
    bool input_points_have_normals = (pwns.begin()->normal() != CGAL::NULL_VECTOR);
    if (input_points_have_normals)
    {
      std::cerr << "Compare with original normals..." << std::endl;
    
      PointList::iterator p;
      /*std::deque<Orientable_normal>::iterator*/ n;
      double min_normal_deviation = DBL_MAX;
      double max_normal_deviation = DBL_MIN;
      double avg_normal_deviation = 0;
      int flipped_normals = 0;
      for (p = pwns.begin(), n = computed_normals.begin(); p != pwns.end(); p++, n++)
      {
        Vector v1 = p->normal(); // input normal
        double norm1 = std::sqrt( v1*v1 );
        assert(norm1 != 0.0);
        Vector v2 = *n; // computed normal
        double norm2 = std::sqrt( v2*v2 );
        assert(norm2 != 0.0);
        double cos_normal_deviation = (v1*v2)/(norm1*norm2);
        if (cos_normal_deviation < 0)
        {
          cos_normal_deviation = -cos_normal_deviation;
          if ( n->is_oriented() ) // unoriented normals are already reported
            flipped_normals++;
        }
        double normal_deviation = std::acos(cos_normal_deviation);
        
        min_normal_deviation = (std::min)(min_normal_deviation, normal_deviation);
        max_normal_deviation = (std::max)(max_normal_deviation, normal_deviation);
        avg_normal_deviation += normal_deviation;
      }
      avg_normal_deviation /= double(pwns.size());
      
      if (flipped_normals > 0)
      {
        std::cerr << "Error: " << flipped_normals << " normal(s) are flipped\n";
        accumulated_fatal_err = EXIT_FAILURE; // set error and continue
      }

      std::cerr << "Min normal deviation=" << min_normal_deviation*180.0/M_PI << " degrees\n";
      std::cerr << "Max normal deviation=" << max_normal_deviation*180.0/M_PI << " degrees\n";
      std::cerr << "Avg normal deviation=" << avg_normal_deviation*180.0/M_PI << " degrees\n";
      if (max_normal_deviation*180.0/M_PI > 10.0) // 10 degrees
        accumulated_fatal_err = EXIT_FAILURE; // set error and continue
    }

    //***************************************
    // Save the point set
    //***************************************
    
    // Replace old normals by new ones
    PointList::iterator p;
    /*std::deque<Orientable_normal>::iterator*/ n;
    for (p = pwns.begin(), n = computed_normals.begin(); p != pwns.end(); p++, n++)
      p->normal() = *n;

    std::cerr << "Write file " << output_filename << std::endl << std::endl;
    
    /*std::string*/ extension = output_filename.substr(output_filename.find_last_of('.'));
    if (extension == ".pwn" || extension == ".PWN")
    {
      if( ! CGAL::surface_reconstruction_write_pwn(output_filename.c_str(),
                                                   pwns.begin(), pwns.end()) )
      {
        std::cerr << "Error: cannot write file " << output_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (extension == ".xyz" || extension == ".XYZ")
    {
      if( ! CGAL::surface_reconstruction_write_xyz(output_filename.c_str(),
                                                   pwns.begin(), pwns.end()) )
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

