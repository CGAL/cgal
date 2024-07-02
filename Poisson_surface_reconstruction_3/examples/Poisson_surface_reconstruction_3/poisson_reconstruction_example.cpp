#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Poisson_reconstruction_function.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Polygon_mesh_processing/distance.h>

#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Labeled_mesh_domain_3<Kernel> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

int main(void)
{
    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle in degrees.
    FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.

    // Reads the point set file in points[].
    // Note: read_points() requires an iterator over points
    // + property maps to access each point's position and normal.
    PointList points;
    if(!CGAL::IO::read_points(CGAL::data_file_path("points_3/kitten.xyz"), std::back_inserter(points),
                          CGAL::parameters::point_map(Point_map())
                                           .normal_map (Normal_map())))
    {
      std::cerr << "Error: cannot read file input file!" << std::endl;
      return EXIT_FAILURE;
    }

    // Creates implicit function from the read points using the default solver.

    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    Poisson_reconstruction_function function(points.begin(), points.end(), Point_map(), Normal_map());

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function() )
      return EXIT_FAILURE;

    // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, 6 /* knn = 1 ring */,
       CGAL::parameters::point_map (Point_map()));

    //Computes implicit function bounding sphere radius.
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance

    // Defines surface mesh generation criteria
    Mesh_criteria criteria(CGAL::parameters::facet_angle = sm_angle,
                           CGAL::parameters::facet_size = sm_radius*average_spacing,
                           CGAL::parameters::facet_distance = sm_distance*average_spacing);

    // Defines mesh domain
    Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(function, bsphere,
        CGAL::parameters::relative_error_bound(sm_dichotomy_error / sm_sphere_radius));

    // Generates mesh with manifold option
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude().no_perturb()
                                        .manifold_with_boundary());

    const Tr& tr = c3t3.triangulation();
    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    // saves reconstructed surface mesh
    std::ofstream out("kitten_poisson-20-30-0.375.off");
    Polyhedron output_mesh;
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, output_mesh);
    out << output_mesh;


    /// [PMP_distance_snippet]
    // computes the approximation error of the reconstruction
    double max_dist =
      CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set
      (output_mesh,
       CGAL::make_range (boost::make_transform_iterator
                         (points.begin(), CGAL::Property_map_to_unary_function<Point_map>()),
                         boost::make_transform_iterator
                         (points.end(), CGAL::Property_map_to_unary_function<Point_map>())),
       4000);
    std::cout << "Max distance to point_set: " << max_dist << std::endl;
    /// [PMP_distance_snippet]

    return EXIT_SUCCESS;
}
