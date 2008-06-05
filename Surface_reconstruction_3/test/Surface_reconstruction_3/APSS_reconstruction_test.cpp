// APSS_reconstruction_test.cpp


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the APSS reconstruction method
// No output
// Input files are .off
//----------------------------------------------------------
// APSS_reconstruction_test mesh1.off mesh2.off...



// CGAL
#include <CGAL/basic.h> // include basic.h before testing #defines
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// Surface mesher
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This package
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/surface_reconstruction_output.h>
#include <CGAL/APSS_implicit_function.h>

// This test
#include "enriched_polyhedron.h"

// STL stuff
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;

typedef std::vector<Point_with_normal> PointList;

// APSS implicit function
typedef CGAL::APSS_implicit_function<Kernel,Point_with_normal> APSS_implicit_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Str;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Str> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, APSS_implicit_function&> Surface_3;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

//// Scale point set to [-1,1]^3
//void reshape(PointList& pwns)
//{
//    Point cmin = pwns[0];
//    Point cmax = pwns[0];
//    for (unsigned int i=1 ; i<pwns.size() ; ++i)
//    {
//        cmin = CGAL::min(cmin,pwns[i]);
//        cmax = CGAL::max(cmax,pwns[i]);
//    }
//
//    Point mid = midpoint(cmax,cmin);
//    Vector diag = cmax-cmin;
//    FT s = 2./(diag.x()>diag.y() ? (diag.x()>diag.z() ? diag.x() : diag.z()) : (diag.y()>diag.z() ? diag.y() : diag.z()));
//    for (unsigned int i=0 ; i<pwns.size() ; ++i)
//    {
//        pwns[i] = CGAL::ORIGIN + s * (pwns[i] - mid);
//    }
//}

// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "RECONSTRUCTION" << std::endl;
    std::cerr << "Test the APSS reconstruction method" << std::endl;
    std::cerr << "No output" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 == 0)
    {
        std::cerr << "Usage: " << argv[0] << " input_file1.off input_file2.obj ..." << std::endl;
        return(EXIT_FAILURE);
    }

    unsigned int nofNeighbors = 10;

    // Accumulated errors
    int accumulated_fatal_err = EXIT_SUCCESS;

    // Reconstruct each input file and accumulate errors
    for (int arg_index = 1; arg_index <= argc-1; arg_index++)
    {
        std::cerr << std::endl;

        // File name is:
        const char* input_filename  = argv[arg_index];

        //***************************************
        // Load mesh
        //***************************************

        // Read the mesh file in a polyhedron
        std::ifstream stream(input_filename);
        typedef Enriched_polyhedron<Kernel,Enriched_items> Polyhedron;
        Polyhedron input_mesh;
        CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
        if(!stream || !input_mesh.is_valid() || input_mesh.empty())
        {
            std::cerr << "FATAL ERROR: cannot read OFF file " << input_filename << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            continue;
        }

        // Compute vertices' normals from connectivity
        input_mesh.compute_normals();

        // Convert vertices and normals to PointList
        PointList pwns;
        Polyhedron::Vertex_iterator v;
        for(v = input_mesh.vertices_begin();
            v != input_mesh.vertices_end();
            v++)
        {
          const Point& p = v->point();
          const Vector& n = v->normal();
          pwns.push_back(Point_with_normal(p,n));
        }

        // Print status
        int nb_vertices = input_mesh.size_of_vertices();
        std::cerr << "Read file " << input_filename << ": "
                  << nb_vertices << " vertices"
                  << std::endl;

        //***************************************
        // Compute implicit function
        //***************************************

        CGAL::Timer task_timer;
        task_timer.start();

        //reshape(pwns); // Scale point set to [-1,1]^3

        APSS_implicit_function apss_function(pwns.begin(), pwns.end());
        apss_function.setNofNeighbors(nofNeighbors);

        //***************************************
        // Surface mesh generation
        //***************************************

        Str tr;           // 3D-Delaunay triangulation
        C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

        // Get inner point
        Point inner_point = apss_function.get_inner_point();
        FT inner_point_value = apss_function(inner_point);
        if(inner_point_value >= 0.0)
        {
            std::cerr << "FATAL ERROR: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            continue;
        }

        // Get implicit surface's size
        Sphere bounding_sphere = apss_function.bounding_sphere();
        FT size = sqrt(bounding_sphere.squared_radius());

        // defining the surface
        Point sm_sphere_center = inner_point; // bounding sphere centered at inner_point
        FT    sm_sphere_radius = 2 * size; 
        sm_sphere_radius *= 1.1; // <= the Surface Mesher fails if the sphere does not contain the surface
        Surface_3 surface(apss_function,
                          Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius));

        // defining meshing criteria
        FT sm_angle = 30.0; // theorical guaranty if angle >= 30
        FT sm_radius = 0.1; // as suggested by LR
        FT sm_distance = 0.005;
        CGAL::Surface_mesh_default_criteria_3<Str> criteria(sm_angle,  // lower bound of facets angles (degrees)
                                                            sm_radius*size,  // upper bound of Delaunay balls radii
                                                            sm_distance*size); // upper bound of distance to surface

        // meshing surface
std::cerr << "make_surface_mesh(sphere={center=("<<sm_sphere_center << "), radius="<<sm_sphere_radius << "},\n"
          << "                  criteria={angle="<<sm_angle << ", radius="<<sm_radius*size << ", distance="<<sm_distance*size << "},\n"
          << "                  Non_manifold_tag())...\n";
        CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

        // Print status
        std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                         << tr.number_of_vertices() << " vertices"
                                         << std::endl;
        task_timer.reset();

    } // for each input file

    std::cerr << std::endl;

    // Return accumulated fatal error
    std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
    return accumulated_fatal_err;
}
