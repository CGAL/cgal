//----------------------------------------------------------
// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
//----------------------------------------------------------

#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Point_set_scene_item.h"

// CGAL
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This package
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>


// Poisson implicit function
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;


// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
Polyhedron* poisson_reconstruct(const Point_set& points,
                                FT sm_angle, // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
                                FT sm_radius, // Max triangle size w.r.t. point set radius. 0.1 is fine.
                                FT sm_distance) // Approximation error w.r.t. p.s.r.. For Poisson: 0.01=fast, 0.002=smooth.
{
    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Checks requirements
    //***************************************

    int nb_points = points.size();
    if (nb_points == 0)
    {
      std::cerr << "Error: empty point set" << std::endl;
      return NULL;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return NULL;
    }

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Creates Poisson triangulation...\n";

    // Creates implicit function from the point set.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(
                              points.begin(), points.end(),
                              CGAL::make_normal_of_point_with_normal_pmap(points.begin()));
    // Prints status
    std::cerr << "Creates Poisson triangulation: " << task_timer.time() << " seconds\n";
    task_timer.reset();

    std::cerr << "Computes implicit function...\n";

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function() )
    {
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      return NULL;
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
      return NULL;
    }

    // Gets implicit function's radius
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    // defining the implicit surface = implicit function + bounding sphere centered at inner_point
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
      return NULL;

    //***************************************
    // Converts to polyhedron
    //***************************************

    Polyhedron* output_mesh = new Polyhedron;
    CGAL::output_surface_facets_to_polyhedron(c2t3, *output_mesh);

    return output_mesh;
}

