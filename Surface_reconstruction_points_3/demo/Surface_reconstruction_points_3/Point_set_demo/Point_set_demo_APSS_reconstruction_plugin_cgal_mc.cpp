//----------------------------------------------------------
// APSS reconstruction method:
// Reconstruct a surface mesh from a point set and return it as a polyhedron.
//----------------------------------------------------------

#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Point_set_scene_item.h"
#include "marching_cubes.h"

// CGAL
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This package
#include <CGAL/APSS_reconstruction_function.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>


// APSS implicit function
typedef CGAL::APSS_reconstruction_function<Kernel> APSS_reconstruction_function;

// APSS reconstruction method:
// Reconstruct a surface mesh from a point set and return it as a polyhedron.
Polyhedron* APSS_reconstruct_mc(const Point_set& points,
                                FT sm_angle, // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
                                FT sm_radius, // Max triangle radius w.r.t. point set radius. 0.1 is fine.
                                FT sm_distance, // Approximation error w.r.t. p.s.r.. For APSS: 0.015 = fast, 0.003 = smooth.
                                FT smoothness,  // smoothness factor
                                int grid_size)
{
    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Check requirements
    //***************************************

    int nb_points = points.size();
    if (nb_points == 0)
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
    // Creates implicit function
    //***************************************

    std::cerr << "Creates APSS implicit function (smoothness=" << smoothness << ")...\n";

    // Creates implicit function
    // Creates implicit function
    APSS_reconstruction_function implicit_function(points.begin(), points.end(),
                                                   CGAL::make_normal_of_point_with_normal_pmap(points.begin()),
                                                   smoothness);

    // Prints status
    std::cerr << "Creates implicit function: " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Marching cubes...\n";

    Polyhedron* output_mesh = new Polyhedron;
    marching_cubes(implicit_function, grid_size, *output_mesh);

    //***************************************
    // Erases small connected components
    //***************************************

//     std::cerr << "Erases small connected components...\n";

//     unsigned int nb_erased_components =
//       CGAL::keep_largest_connected_components(*output_mesh, 1/* keep largest component only*/);

    // Prints status
//     std::cerr << "Erases small connected components: " << task_timer.time() << " seconds, "
//                                                       << nb_erased_components << " components erased"
//                                                       << std::endl;
    task_timer.reset();

    return output_mesh;
}

