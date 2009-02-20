// APSS_reconstruction.cpp

//----------------------------------------------------------
// APSS reconstruction method:
// Read a point set or a mesh's set of vertices, reconstruct a surface,
// and save the surface.
//----------------------------------------------------------
// APSS_reconstruction file_in file_out [options]


// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#define CGAL_C2T3_USE_FILE_WRITER_OFF
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// This package
#include <CGAL/APSS_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>
#include <CGAL/IO/surface_reconstruction_read_pwn.h>

// This test
#include "enriched_polyhedron.h"

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
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;

typedef std::deque<Point_with_normal> PointList;

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
    std::cerr << "RECONSTRUCTION" << std::endl;
    std::cerr << "APSS reconstruction method." << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc<3)
    {
      std::cerr << "Read a point set or a mesh's set of vertices, reconstruct a surface,\n";
      std::cerr << "and save the surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .off.\n";
      std::cerr << "Options:\n";
      std::cerr << "  -sm_radius <float>     Radius upper bound (default=0.1 * point set radius)\n";
      std::cerr << "  -sm_distance <float>   Distance upper bound (default=0.003 * point set radius)\n";
      std::cerr << "  -k <int>               Number of neighbors (default=24)\n";
      std::cerr << "                           - should be greater than 7,\n";
      std::cerr << "                           - high numbers lead to smoother surfaces.\n";
      return EXIT_FAILURE;
    }

    // APSS options
    FT sm_angle_apss = 20.0; // theorical guaranty if angle >= 30, but slower
    FT sm_radius_apss = 0.1; // as suggested by LR
    FT sm_error_bound_apss = 1e-3;
    FT sm_distance_apss = 0.003; // Upper bound of distance to surface (APSS).
                                 // Note: 1.5 * Poisson's distance gives roughly the same number of triangles.
    unsigned int nb_neighbors_apss = 24; // APSS K-nearest neighbors = 4 rings, as suggested by GG

    // decode parameters
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];
    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i])=="-sm_radius")
        sm_radius_apss = atof(argv[++i]);
      else if (std::string(argv[i])=="-sm_distance")
        sm_distance_apss = atof(argv[++i]);
      else if (std::string(argv[i])=="-k")
        nb_neighbors_apss = atoi(argv[++i]);
      else
        std::cerr << "invalid option " << argv[i] << "\n";
    }

    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Load mesh/point set
    //***************************************

    PointList points;

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
        points.push_back(Point_with_normal(p,n));
      }
    }
    else if (extension == ".xyz" || extension == ".XYZ")
    {
      // Read the point set file in points[]
      if(!CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
                                                std::back_inserter(points)))
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (extension == ".pwn" || extension == ".PWN")
    {
      // Read the point set file in points[]
      if(!CGAL::surface_reconstruction_read_pwn(input_filename.c_str(),
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
    long memory = CGAL::Memory_sizer().virtual_size();
    int nb_vertices = points.size();
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

    assert(points.begin() != points.end());
    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Compute implicit function
    //***************************************

    std::cerr << "Compute APSS implicit function (knn=" << nb_neighbors_apss << ")...\n";

    // Create implicit function
    APSS_reconstruction_function apss_function(points.begin(), points.end(),
                                               nb_neighbors_apss);

    // Recover memory used by points[]
    points.clear();

    // Print status
    /*long*/ memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Compute implicit function: " << task_timer.time() << " seconds, "
                                               << (memory>>20) << " Mb allocated"
                                               << std::endl;
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    STr tr; // 3D-Delaunay triangulation for Surface Mesher
    C2t3 surface_mesher_c2t3 (tr); // 2D-complex in 3D-Delaunay triangulation

    // Get inner point
    Point inner_point = apss_function.get_inner_point();
    FT inner_point_value = apss_function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      return EXIT_FAILURE;
    }

    // Get implicit surface's radius
    Sphere bounding_sphere = apss_function.bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the surface
    Point sm_sphere_center = inner_point; // bounding sphere centered at inner_point
    FT    sm_sphere_radius = 2 * size;
    sm_sphere_radius *= 1.1; // <= the Surface Mesher fails if the sphere does not contain the surface
    Surface_3 surface(apss_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius),
                      sm_error_bound_apss*size/sm_sphere_radius); // dichotomy stops when segment < sm_error_bound_apss*size

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle_apss,  // lower bound of facets angles (degrees)
                                                        sm_radius_apss*size,  // upper bound of Delaunay balls radii
                                                        sm_distance_apss*size); // upper bound of distance to surface

    CGAL_TRACE_STREAM << "  make_surface_mesh(dichotomy error="<<sm_error_bound_apss<<" * point set radius,\n"
                      << "                    sphere center=("<<sm_sphere_center << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius/size<<" * p.s.r.,\n"
                      << "                    angle="<<sm_angle_apss << " degrees,\n"
                      << "                    radius="<<sm_radius_apss<<" * p.s.r.,\n"
                      << "                    distance="<<sm_distance_apss<<" * p.s.r.,\n"
                      << "                    Manifold_with_boundary_tag)\n";

    // meshing surface
    CGAL::make_surface_mesh(surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    // Print status
    /*long*/ memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices, "
                                     << (memory>>20) << " Mb allocated"
                                     << std::endl;
    task_timer.reset();

    //***************************************
    // save the mesh
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    std::ofstream out(output_filename.c_str());
    CGAL::output_surface_facets_to_off(out, surface_mesher_c2t3);

    return EXIT_SUCCESS;
}

