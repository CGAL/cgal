// APSS_reconstruction_test.cpp

//----------------------------------------------------------
// Test the APSS reconstruction method.
// For each input point set or mesh's set of vertices, reconstruct a surface.
// No output.
//----------------------------------------------------------
// APSS_reconstruction_test mesh1.off point_set2.xyz...

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This package
#include <CGAL/APSS_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/point_set_property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

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
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;

// polyhedron
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// APSS implicit function
typedef CGAL::APSS_reconstruction_function<Kernel> APSS_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, APSS_reconstruction_function> Surface_3;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Test the APSS reconstruction method" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  // usage
  if (argc-1 == 0)
  {
      std::cerr << "For each input point set or mesh's set of vertices, reconstruct a surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " mesh1.off point_set2.xyz..." << std::endl;
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "No output" << std::endl;
      return EXIT_FAILURE;
  }

  // APSS options
  FT sm_angle = 20.0; // Min triangle angle (degrees). 20=fast, 30 guaranties convergence (PA).
  FT sm_radius = 0.1; // Max triangle radius w.r.t. point set radius. 0.1 is fine (LR).
  FT sm_distance = 0.015; // Approximation error w.r.t. p.s.r. For APSS: 0.015=fast, 0.003=smooth (LS).
                          // Note: 1.5 * Poisson's  approximation error gives roughly the same number of triangles.
  FT smoothness = 2; // Smoothness factor. In the range 2 (clean datasets) and 8 (noisy datasets).

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for (int i = 1; i <= argc-1; i++)
  {
    CGAL::Timer task_timer; task_timer.start();

    std::cerr << std::endl;

    //***************************************
    // Loads mesh/point set
    //***************************************

    // File name is:
    std::string input_filename  = argv[i];

    PointList points;

    // If OFF file format
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      // Reads the mesh file in a polyhedron
      std::ifstream stream(input_filename.c_str());
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        accumulated_fatal_err = EXIT_FAILURE;
        continue;
      }

      // Converts Polyhedron vertices to point set.
      // Computes vertices normal from connectivity.
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
      // Reads the point set file in points[].
      // Note: read_xyz_points_and_normals() requires an iterator over points
      // + property maps to access each point's position and normal.
      // The position property map can be omitted here as we use iterators over Point_3 elements.
      std::ifstream stream(input_filename.c_str());
      if (!stream ||
          !CGAL::read_xyz_points_and_normals(
                                stream,
                                std::back_inserter(points),
                                CGAL::make_normal_vector_property_map(std::back_inserter(points))))
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        accumulated_fatal_err = EXIT_FAILURE;
        continue;
      }
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    // Prints status
    long memory = CGAL::Memory_sizer().virtual_size();
    int nb_points = points.size();
    std::cerr << "Reads file " << input_filename << ": " << nb_points << " points, "
                                                        << task_timer.time() << " seconds, "
                                                        << (memory>>20) << " Mb allocated"
                                                        << std::endl;
    task_timer.reset();

    //***************************************
    // Check requirements
    //***************************************

    if (nb_points == 0)
    {
      std::cerr << "Error: empty file" << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      // this is not a bug => do not set accumulated_fatal_err
      continue;
    }

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Computes APSS implicit function (smoothness=" << smoothness << ")...\n";

    // Creates implicit function and insert points.
    // Note: APSS_implicit_function() requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    APSS_reconstruction_function implicit_function(
                              points.begin(), points.end(),
                              CGAL::make_normal_vector_property_map(points.begin()),
                              smoothness);

    // Recover memory used by points[]
    points.clear();

    // Prints status
    /*long*/ memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Computes implicit function: " << task_timer.time() << " seconds, "
                                               << (memory>>20) << " Mb allocated"
                                               << std::endl;
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Gets point inside the implicit surface
    Point inner_point = implicit_function.get_inner_point();
    FT inner_point_value = implicit_function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    // Gets implicit function's radius
    Sphere bsphere = implicit_function.bounding_sphere();
    FT size = sqrt(bsphere.squared_radius());

    // defining the implicit surface = implicit function + bounding sphere centered at inner_point
    Point sm_sphere_center = inner_point;
    FT    sm_sphere_radius = size + std::sqrt(CGAL::squared_distance(bsphere.center(),inner_point));
    sm_sphere_radius *= 1.01; // make sure that the bounding sphere contains the surface
    Surface_3 surface(implicit_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius));

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*size,  // Max triangle radius
                                                        sm_distance*size); // Approximation error

    CGAL_TRACE_STREAM << "  make_surface_mesh(sphere center=("<<sm_sphere_center << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius/size<<" * p.s.r.,\n"
                      << "                    angle="<<sm_angle << " degrees,\n"
                      << "                    radius="<<sm_radius<<" * p.s.r.,\n"
                      << "                    distance="<<sm_distance<<" * p.s.r.,\n"
                      << "                    Manifold_with_boundary_tag)\n";

    // meshing surface
    STr tr; // 3D-Delaunay triangulation for Surface Mesher
    C2t3 surface_mesher_c2t3 (tr); // 2D-complex in 3D-Delaunay triangulation
    CGAL::make_surface_mesh(surface_mesher_c2t3, // reconstructed mesh
                            surface, // implicit surface
                            criteria, // meshing criteria
                            CGAL::Manifold_with_boundary_tag()); // require manifold mesh

    // Prints status
    /*long*/ memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices, "
                                     << (memory>>20) << " Mb allocated"
                                     << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0) {
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    // Converts to polyhedron
    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(surface_mesher_c2t3, output_mesh);

    //***************************************
    // Erases small connected components
    //***************************************

    std::cerr << "Erases small connected components...\n";

    unsigned int nb_erased_components =
      output_mesh.keep_largest_connected_components( 1 /* keep largest component only*/ );

    // Prints status
    /*long*/ memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Erases small connected components: " << task_timer.time() << " seconds, "
                                                      << nb_erased_components << " component(s) erased, "
                                                      << (memory>>20) << " Mb allocated"
                                                      << std::endl;
    task_timer.reset();

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;
}

