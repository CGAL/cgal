// poisson_reconstruction.cpp

//----------------------------------------------------------
// Poisson Delaunay Reconstruction method.
// Reads a point set or a mesh's set of vertices, reconstructs a surface using Poisson,
// and saves the surface.
// Output format is .off.
//----------------------------------------------------------
// poisson_reconstruction file_in file_out [options]

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// This package
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/read_xyz_points.h>

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

// Poisson implicit function
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

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
    if (argc-1 < 2)
    {
      std::cerr << "Reads a point set or a mesh's set of vertices, reconstructs a surface using Poisson,\n";
      std::cerr << "and saves the surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .off.\n";
      std::cerr << "Options:\n";
      std::cerr << "  -sm_radius <float>     Radius upper bound (default=0.1 * point set radius)\n";
      std::cerr << "  -sm_distance <float>   Distance upper bound (default=0.002 * point set radius)\n";
      return EXIT_FAILURE;
    }

    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees). 20=fast, 30 guaranties convergence.
    FT sm_radius = 0.1; // Max triangle size w.r.t. point set radius. 0.1 is fine.
    FT sm_distance = 0.002; // Approximation error w.r.t. p.s.r. For Poisson: 0.01=fast, 0.002=smooth.

    // decode parameters
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];
    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i])=="-sm_radius")
        sm_radius = atof(argv[++i]);
      else if (std::string(argv[i])=="-sm_distance")
        sm_distance = atof(argv[++i]);
      else
        std::cerr << "invalid option " << argv[i] << "\n";
    }

    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Loads mesh/point set
    //***************************************

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
        return EXIT_FAILURE;
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
                                CGAL::make_normal_of_point_with_normal_pmap(std::back_inserter(points))))
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

    // Prints status
    int nb_points = points.size();
    std::cerr << "Reads file " << input_filename << ": " << nb_points << " points, "
                                                        << task_timer.time() << " seconds"
                                                        << std::endl;
    task_timer.reset();

    //***************************************
    // Checks requirements
    //***************************************

    if (nb_points == 0)
    {
      std::cerr << "Error: empty point set" << std::endl;
      return EXIT_FAILURE;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Creates Poisson triangulation...\n";

    // Creates implicit function from the read points.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(
                              points.begin(), points.end(),
                              CGAL::make_normal_of_point_with_normal_pmap(points.begin()));

    // Recover memory used by points[]
    points.clear();

    // Prints status
    std::cerr << "Creates Poisson triangulation: " << task_timer.time() << " seconds\n";
    task_timer.reset();

    std::cerr << "Computes implicit function...\n";

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function() )
    {
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      return EXIT_FAILURE;
    }

    // Prints status
    std::cerr << "Computes implicit function: " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Gets one point inside the implicit surface
    Point inner_point = function.get_inner_point();
    FT inner_point_value = function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      return EXIT_FAILURE;
    }

    // Gets implicit function's radius
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface = implicit function + bounding sphere centered at inner_point
    FT sm_sphere_radius = radius + std::sqrt(CGAL::squared_distance(bsphere.center(),inner_point));
    sm_sphere_radius *= 1.01; // make sure that the bounding sphere contains the surface
    FT sm_dichotomy_error = sm_distance/10.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*radius,  // Max triangle size
                                                        sm_distance*radius); // Approximation error

    CGAL_TRACE_STREAM << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius<<",\n"
                      << "                    dichotomy error="<<sm_dichotomy_error<<" * sphere radius,\n"
                      << "                    angle="<<sm_angle << " degrees,\n"
                      << "                    triangle size="<<sm_radius<<" * point set radius,\n"
                      << "                    distance="<<sm_distance<<" * p.s.r.,\n"
                      << "                    Manifold_with_boundary_tag)\n"
                      << "  where point set radius="<<radius<<"\n";

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                // reconstructed mesh
                            surface,                             // implicit surface
                            criteria,                            // meshing criteria
                            CGAL::Manifold_with_boundary_tag()); // require manifold mesh

    // Prints status
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices"
                                     << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    //***************************************
    // saves reconstructed surface mesh
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    std::ofstream out(output_filename.c_str());
    CGAL::output_surface_facets_to_off(out, c2t3);

    return EXIT_SUCCESS;
}
