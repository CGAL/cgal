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
#include <CGAL/IO/surface_reconstruction_output.h>
#include <CGAL/apss.h>

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
typedef Kernel::Sphere_3 Sphere;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Str;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Str> C2t3;
//typedef CGAL::Implicit_surface_3<Kernel, APSS_implicit_function&> Surface_3;
typedef std::vector<Point> PointList;
typedef CGAL::APSS<Kernel,PointList,PointList> Surface_3;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// Scale point set to [-1,1]^3
void reshape(PointList& points)
{
    Point cmin = points[0];
    Point cmax = points[0];
    for (unsigned int i=1 ; i<points.size() ; ++i)
    {
        cmin = CGAL::min(cmin,points[i]);
        cmax = CGAL::max(cmax,points[i]);
    }

    Point mid = midpoint(cmax,cmin);
    Vector diag = cmax-cmin;
    FT s = 2./(diag.x()>diag.y() ? (diag.x()>diag.z() ? diag.x() : diag.z()) : (diag.y()>diag.z() ? diag.y() : diag.z()));
    for (unsigned int i=0 ; i<points.size() ; ++i)
    {
        points[i] = Surface_3::mul(s,Surface_3::sub(points[i],mid));
    }
}

// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "RECONSTRUCTION" << std::endl;
    std::cerr << "Test the APSS reconstruction method" << std::endl;
    std::cerr << "No output" << std::endl;

    //--------------------------------------------------------------------------------
    // parse arguments
    //--------------------------------------------------------------------------------

    if (argc-1 == 0)
    {
        std::cerr << "Usage: " << argv[0] << " input_file1.off input_file2.obj ..." << std::endl;
        return(EXIT_FAILURE);
    }

    FT sm_radius = 0.07;
    FT sm_distance = 0.02;
    unsigned int nofNeighbors = 10;

    // Accumulated errors
    int accumulated_fatal_err = EXIT_SUCCESS;

    // Reconstruct each input file and accumulate errors
    for (int arg_index = 1; arg_index <= argc-1; arg_index++)
    {
        std::cerr << std::endl;

        // File name is:
        const char* input_filename  = argv[arg_index];

        //--------------------------------------------------------------------------------
        // input
        //--------------------------------------------------------------------------------

        CGAL::Timer task_timer;
        task_timer.start();

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

        // Convert vertices and normals to PointList/NormalList
        PointList points, normals;
        Polyhedron::Vertex_iterator v;
        for(v = input_mesh.vertices_begin();
            v != input_mesh.vertices_end();
            v++)
        {
          const Point& p = v->point();
          const Vector& n = v->normal();
          points.push_back(p);
          normals.push_back(CGAL::ORIGIN + n);
        }

        // Print status
        int nb_vertices = input_mesh.size_of_vertices();
        std::cerr << "Read file " << input_filename << ": "
                  << task_timer.time() << " seconds, "
                  << nb_vertices << " vertices"
                  << std::endl;
        task_timer.reset();

        //--------------------------------------------------------------------------------
        // create the MLS surface
        //--------------------------------------------------------------------------------

        reshape(points); // Scale point set to [-1, 1]^3

        Surface_3 surface(points, normals);

        surface.setNofNeighbors(nofNeighbors);

        //--------------------------------------------------------------------------------
        // Do the polygonization
        //--------------------------------------------------------------------------------

        Str tr;            // 3D-Delaunay triangulation
        C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

        // defining meshing criteria
        CGAL::Surface_mesh_default_criteria_3<Str> criteria( 20.,          // angular bound
                                                            sm_radius,    // radius upper bound
                                                            sm_distance); // distance upper bound

        // meshing surface
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
