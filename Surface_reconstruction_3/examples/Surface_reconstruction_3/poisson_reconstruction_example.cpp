// poisson_reconstruction_example.cpp

//----------------------------------------------------------
// Poisson Delaunay Reconstruction method.
// Read a point set, reconstruct a surface using Poisson,
// and save the surface.
// Input format is .xyz (with normals).
// Output format is .off.
//----------------------------------------------------------
// poisson_reconstruction_example file_in.xyz file_out.off

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#define CGAL_C2T3_USE_FILE_WRITER_OFF
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// This package
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
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
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;

// Poisson's Delaunay triangulation 3 and implicit function
typedef CGAL::Reconstruction_triangulation_3<Kernel> Dt3;
typedef CGAL::Poisson_reconstruction_function<Kernel, Dt3> Poisson_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Poisson Delaunay Reconstruction method" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 != 2)
    {
      std::cerr << "Read a point set, reconstruct a surface using Poisson,\n";
      std::cerr << "and save the surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in.xyz file_out.off" << std::endl;
      std::cerr << "Input format is .xyz (with normals).\n";
      std::cerr << "Output file format is .off.\n";
      return EXIT_FAILURE;
    }

    // Poisson options
    FT sm_angle_poisson = 20.0; // theorical guaranty if angle >= 30, but slower
    FT sm_radius_poisson = 0.1;
    FT sm_error_bound_poisson = 1e-3;
    FT sm_distance_poisson = 0.002; // upper bound of distance to surface (Poisson)

    // decode parameters
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];

    //***************************************
    // Load point set
    //***************************************

    PointList points;

    // Read the point set file in points[]
    std::cerr << "Open " << input_filename << " for reading...";
    if(CGAL::surface_reconstruction_read_xyz(input_filename.c_str(),
                                             std::back_inserter(points)))
    {
      std::cerr << "ok (" << points.size() << " points)" << std::endl;
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Compute implicit function
    //***************************************

    // Create implicit function and triangulation.
    // Insert vertices and normals in triangulation.
    Dt3 dt;
    Poisson_reconstruction_function poisson_function(dt, points.begin(), points.end());

    /// Compute the Poisson indicator function f()
    /// at each vertex of the triangulation.
    if ( ! poisson_function.compute_implicit_function() )
    {
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    STr tr; // 3D-Delaunay triangulation for Surface Mesher
    C2t3 surface_mesher_c2t3 (tr); // 2D-complex in 3D-Delaunay triangulation

    // Get inner point
    Point inner_point = poisson_function.get_inner_point();

    // Get implicit surface's radius
    Sphere bounding_sphere = poisson_function.bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the surface
    Point sm_sphere_center = inner_point; // bounding sphere centered at inner_point
    FT    sm_sphere_radius = 2 * size;
    sm_sphere_radius *= 1.1; // <= the Surface Mesher fails if the sphere does not contain the surface
    Surface_3 surface(poisson_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius),
                      sm_error_bound_poisson*size/sm_sphere_radius); // dichotomy stops when segment < sm_error_bound_poisson*size

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle_poisson,  // lower bound of facets angles (degrees)
                                                        sm_radius_poisson*size,  // upper bound of Delaunay balls radii
                                                        sm_distance_poisson*size); // upper bound of distance to surface

    // meshing surface
    CGAL::make_surface_mesh(surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    // Print status
    std::cerr << "Surface meshing: " << tr.number_of_vertices() << " output vertices"
                                     << std::endl;

    //***************************************
    // save the mesh
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    std::ofstream out(output_filename.c_str());
    CGAL::output_surface_facets_to_off(out, surface_mesher_c2t3);

    return EXIT_SUCCESS;
}

