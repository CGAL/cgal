// poisson_reconstruction_test.cpp


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test the Poisson Delaunay Reconstruction method
// No output
// Input files are .off
//----------------------------------------------------------
// poisson_reconstruction_test mesh1.off mesh2.off...


#include <CGAL/basic.h> // include basic.h before testing #defines

#ifdef CGAL_USE_TAUCS


// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Surface mesher
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
//#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// This package
#include <CGAL/surface_reconstruction_output.h>
#include <CGAL/Poisson_implicit_function.h>

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
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Tetrahedron_3 Tetrahedron;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;

// Poisson's Delaunay triangulation 3 and implicit function
typedef CGAL::Implicit_fct_delaunay_triangulation_3<Kernel> Dt3;
typedef CGAL::Poisson_implicit_function<Kernel, Dt3> Poisson_implicit_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Str;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Str> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_implicit_function&> Surface_3;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "RECONSTRUCTION" << std::endl;
    std::cerr << "Test the Poisson Delaunay Reconstruction method" << std::endl;
    std::cerr << "No output" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 == 0)
    {
        std::cerr << "Usage: " << argv[0] << " input_file1.off input_file2.obj ..." << std::endl;
        return(EXIT_FAILURE);
    }

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

        // Compute vertices normals
        input_mesh.compute_normals();

        // Insert vertices in triangulation
        Dt3 dt;
		    std::vector<Point_with_normal> pwns;
		    Polyhedron::Vertex_iterator v;
		    for(v = input_mesh.vertices_begin();
				    v != input_mesh.vertices_end();
				    v++)
		    {
			    const Point& p = v->point();
			    const Vector& n = v->normal();
			    pwns.push_back(Point_with_normal(p,n));
		    }
        dt.insert(pwns.begin(), pwns.end(), Dt3::INPUT);

        // Print status
        int nb_vertices = input_mesh.size_of_vertices();
        std::cerr << "Read file " << input_filename << ": "
                  << task_timer.time() << " seconds, "
                  << nb_vertices << " vertices"
                  << std::endl;
        task_timer.reset();

        //***************************************
        // Solve Poisson equation
        //***************************************

        Poisson_implicit_function poisson_function(dt);

        /// Computes the Poisson indicator function f()
        /// at each vertex of the triangulation
	      if ( ! poisson_function.compute_implicit_function() )
        {
            std::cerr << "FATAL ERROR: cannot solve Poisson equation" << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            continue;
        }

        // Print status
        int nb_vertices2 = poisson_function.triangulation().number_of_vertices();
        std::cerr << "Solve Poisson equation: "
                  << task_timer.time() << " seconds "
                  << "(added " << nb_vertices2-nb_vertices << " vertices)"
                  << std::endl;
        task_timer.reset();

        //***************************************
        // Surface mesh generation
        //***************************************

        Str tr;           // 3D-Delaunay triangulation
        C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

        // Get inner point
	      Point sink = poisson_function.sink();
	      FT f_sink = poisson_function(sink);
	      if(f_sink >= 0.0)
	      {
            std::cerr << "FATAL ERROR: unable to seed (" << f_sink << " at sink)" << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            continue;
	      }

        // Get implicit surface's size
        Sphere bounding_sphere = poisson_function.bounding_sphere();
		    FT size = sqrt(bounding_sphere.squared_radius());

        // defining the surface
	      Surface_3 surface(poisson_function,
                          Sphere(sink,4*size*size)); // bounding sphere

        // defining meshing criteria
	      FT sm_angle = 20.0; // LR: 30 is OK
	      FT sm_radius = 0.1; // as suggested by LR (was 0.01)
	      FT sm_distance = 0.001; // was 0.01
        CGAL::Surface_mesh_default_criteria_3<Str> criteria(sm_angle,  // lower bound of facets angles (degrees)
                                                            sm_radius*size,  // upper bound of Delaunay balls radii
                                                            sm_distance*size); // upper bound of distance to surface

	      // meshing surface
        make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

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


#else // CGAL_USE_TAUCS


#include <iostream>
#include <cstdlib>

// ----------------------------------------------------------------------------
// Empty main() if TAUCS is not installed
// ----------------------------------------------------------------------------

int main()
{
    std::cerr << "Skip test as TAUCS is not installed" << std::endl;
    return EXIT_SUCCESS;
}

#endif // CGAL_USE_TAUCS


