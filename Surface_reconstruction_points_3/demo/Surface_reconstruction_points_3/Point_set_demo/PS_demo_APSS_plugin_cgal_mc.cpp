//----------------------------------------------------------
// APSS reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
//----------------------------------------------------------

#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Point_set_scene_item.h"
#include "marching_cubes.h"

// CGAL
#include <CGAL/Timer.h>

// This package
#include <CGAL/APSS_reconstruction_function.h>


// APSS implicit function
typedef CGAL::APSS_reconstruction_function<Kernel> APSS_reconstruction_function;

// APSS reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
Polyhedron* APSS_reconstruct_mc(const Point_set& points,
                                FT smoothness,  // smoothness factor
                                int grid_size)
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
    // Creates implicit function
    //***************************************

    std::cerr << "Creates APSS implicit function (smoothness=" << smoothness << ")...\n";

    // Creates implicit function
    APSS_reconstruction_function function(points.begin(), points.end(),
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
    marching_cubes(function, grid_size, *output_mesh);

    // Prints status
    std::cerr << "Marching cubes: " << task_timer.time() << " seconds, "
                                    << output_mesh->size_of_vertices() << " output vertices"
                                    << std::endl;
    task_timer.reset();

    return output_mesh;
}

