//----------------------------------------------------------
// APSS reconstruction method:
// Reconstruct a surface mesh from a point set and return it as a polyhedron.
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
#include <CGAL/APSS_reconstruction_function.h>
#include <CGAL/IO/surface_reconstruction_output_surface_facets.h>
#include <CGAL/polyhedron_connected_components.h>


// APSS implicit function
typedef CGAL::APSS_reconstruction_function<Kernel> APSS_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, APSS_reconstruction_function> Surface_3;


// APSS reconstruction method:
// Reconstruct a surface mesh from a point set and return it as a polyhedron.
Polyhedron* APSS_reconstruct(const Point_set& points,
                             FT sm_angle, // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
                             FT sm_radius, // Max triangle radius w.r.t. point set radius. 0.1 is fine.
                             FT sm_distance, // Approximation error w.r.t. p.s.r.. For APSS: 0.015 = fast, 0.003 = smooth.
                             unsigned int k = 24) // #neighbors to compute APPS sphere fitting. 12 = fast, 24 = robust.
{
    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Check requirements
    //***************************************

    int nb_vertices = points.size();
    if (nb_vertices == 0)
    {
      std::cerr << "Error: empty file" << std::endl;
      return NULL;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return NULL;
    }

    //***************************************
    // Compute implicit function
    //***************************************

    std::cerr << "Compute APSS implicit function (k=" << k << ")...\n";

    // Create implicit function
    APSS_reconstruction_function implicit_function(points.begin(), points.end(),
                                                   k);

    // Print status
    std::cerr << "Compute implicit function: " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Get point inside the implicit surface
    Point inner_point = implicit_function.get_inner_point();
    FT inner_point_value = implicit_function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      return NULL;
    }

    // Get implicit function's radius
    Sphere bounding_sphere = implicit_function.bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the implicit surface = implicit function + bounding sphere centered at inner_point
    Point sm_sphere_center = inner_point;
    FT    sm_sphere_radius = size + std::sqrt(CGAL::squared_distance(bounding_sphere.center(),inner_point));
    sm_sphere_radius *= 1.01; // <= the Surface Mesher fails if the sphere does not contain the surface
    Surface_3 surface(implicit_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius));

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*size,  // Max triangle radius
                                                        sm_distance*size); // Approximation error

    // meshing surface
    STr tr; // 3D-Delaunay triangulation for Surface Mesher
    C2t3 surface_mesher_c2t3 (tr); // 2D-complex in 3D-Delaunay triangulation
    CGAL::make_surface_mesh(surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    // Print status
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices"
                                     << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0)
      return NULL;

    // Convert to polyhedron
    Polyhedron* output_mesh = new Polyhedron;
    surface_reconstruction_output_surface_facets(surface_mesher_c2t3, *output_mesh);

    //***************************************
    // Erase small connected components
    //***************************************

    std::cerr << "Erase small connected components...\n";
    
    unsigned int nb_erased_components = 
      erase_small_polyhedron_connected_components(*output_mesh);

    // Print status
    std::cerr << "Erase small connected components: " << task_timer.time() << " seconds, "
                                                      << nb_erased_components << " components erased"
                                                      << std::endl;
    task_timer.reset();

    return output_mesh;
}

