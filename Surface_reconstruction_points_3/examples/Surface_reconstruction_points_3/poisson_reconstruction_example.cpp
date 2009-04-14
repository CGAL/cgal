#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#define CGAL_C2T3_USE_FILE_WRITER_OFF
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/read_xyz_point_set.h>

#include <deque>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;
typedef CGAL::Reconstruction_triangulation_3<Kernel> Dt3;
typedef CGAL::Poisson_reconstruction_function<Kernel, Dt3> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

int main(void)
{
    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
    FT sm_radius = 0.1; // Max triangle radius w.r.t. point set radius. 0.1 is fine.
    FT sm_distance = 0.01; // Approximation error w.r.t. p.s.r. For Poisson: 0.01 = fast, 0.002 = smooth.

    // Read a .xyz point set file in points[]
    PointList points;
    std::ifstream stream("data/oni.xyz");
    if (!stream || 
        !CGAL::read_xyz_point_set(stream,
                                  std::back_inserter(points)))
    {
      return EXIT_FAILURE;
    }

    // Create implicit function.
    // Create 3D-Delaunay triangulation for the implicit function and insert vertices.
    Dt3 dt;
    Poisson_reconstruction_function implicit_function(dt, points.begin(), points.end());

    // Compute the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! implicit_function.compute_implicit_function() )
      return EXIT_FAILURE;

    // Get point inside the implicit surface
    Point inner_point = implicit_function.get_inner_point();

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

    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    //***************************************
    // save the mesh
    //***************************************

    std::ofstream out("oni_poisson.off");
    CGAL::output_surface_facets_to_off(out, surface_mesher_c2t3);

    return EXIT_SUCCESS;
}

