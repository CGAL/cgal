#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#define CGAL_C2T3_USE_FILE_WRITER_OFF
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/surface_reconstruction_read_xyz.h>

#include <deque>
#include <fstream>

// types
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
    // Read a .xyz point set file in points[]
    PointList points;
    if (!CGAL::surface_reconstruction_read_xyz("data/oni.xyz",
                                               std::back_inserter(points)))
    {
      return EXIT_FAILURE;
    }

    // Poisson options
    FT sm_angle_poisson = 20.0; // Theorical guaranty if angle >= 30, but slower
    FT sm_radius_poisson = 0.1; // Upper bound of Delaunay balls radii. 0.1 is fine (LR).
    FT sm_distance_poisson = 0.004; // Upper bound of distance to surface (Poisson). 0.004 = fast, 0.002 = smooth.

    // Create implicit function.
    // Create 3D-Delaunay triangulation for the implicit function and insert vertices.
    Dt3 dt;
    Poisson_reconstruction_function poisson_function(dt, points.begin(), points.end());

    /// Compute the Poisson indicator function f() at each vertex of the triangulation.
    poisson_function.compute_implicit_function();

    // Get point inside the implicit surface
    Point inner_point = poisson_function.get_inner_point();

    // Get implicit surface's radius
    Sphere bounding_sphere = poisson_function.bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the surface
    Point sm_sphere_center = inner_point; // bounding sphere centered at inner_point
    FT    sm_sphere_radius = 2 * size;
    sm_sphere_radius *= 1.1; // <= the Surface Mesher fails if the sphere does not contain the surface
    Surface_3 surface(poisson_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius));

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle_poisson,
                                                        sm_radius_poisson*size,
                                                        sm_distance_poisson*size);

    // meshing surface
    STr tr; // 3D-Delaunay triangulation for Surface Mesher
    C2t3 surface_mesher_c2t3 (tr); // 2D-complex in 3D-Delaunay triangulation
    CGAL::make_surface_mesh(surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    //***************************************
    // save the mesh
    //***************************************

    std::ofstream out("oni_poisson.off");
    CGAL::output_surface_facets_to_off(out, surface_mesher_c2t3);

    return EXIT_SUCCESS;
}

