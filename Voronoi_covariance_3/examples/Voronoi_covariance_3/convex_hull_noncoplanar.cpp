#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>
#include <CGAL/Voronoi_covariance_3/halfspaces_intersection.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Plane_3 Plane;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

typedef CGAL::Voronoi_covariance_3::Convex_hull_traits_dual_3<K>    Hull_traits_dual_3;
typedef CGAL::Polyhedron_3<Hull_traits_dual_3>                      Polyhedron_dual_3;

#include <CGAL/convex_hull_3.h>
#include "include/to_dual.h"

#include <CGAL/point_generators_3.h>

template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);
    v = v / sqrt(v.squared_length());

    return plane;
}

// Write into an OFF file to visualize with GeomView
template <typename K, typename Polyhedron>
void convertToOFF (std::string const& filename, Polyhedron& P) {
    typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
    typedef typename Polyhedron::Halfedge_handle                   Halfedge_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

    std::ofstream file(filename.c_str());

    // 0. number of vertices / number of facets / number of edges
    file << "OFF" << std::endl << P.size_of_vertices() << ' '
        << P.size_of_facets() << " 0" << std::endl;

    // 1. vertices definition
    std::copy( P.points_begin(), P.points_end(),
               std::ostream_iterator<typename K::Point_3>( file, "\n"));

    // 2. facets definition
    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        file << std::endl;
    }

    file.close();
}

#include <cstdlib>

int main (int argc, char *argv[]) {
    // define dual polyhedrons to hold convex hulls
    Polyhedron_dual_3 dual_poly_sphere;
    Polyhedron_dual_3 dual_poly_cube;

    // primal polyhedrons
    Polyhedron_3 poly_cube;
    Polyhedron_3 poly_sphere;

    // traits
    // TODO
    Point origin(0, 0, 0);
    Hull_traits_dual_3 dual_traits(origin);

    // Cube
    // IMPORTANT: d <= 0
    std::list<Plane> planes;
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));
    planes.push_back(Plane(0, 1, 0, -1));
    planes.push_back(Plane(0, -1, 0, -1));
    planes.push_back(Plane(0, 0, 1, -1));
    planes.push_back(Plane(0, 0, -1, -1));

    // Translated cube
    /* planes.push_back(Plane(1, 0, 0, -1)); */
    /* planes.push_back(Plane(-1, 0, 0, -1)); */
    /* planes.push_back(Plane(0, 1, 0, -1)); */
    /* planes.push_back(Plane(0, -0.5, 0, -1)); */
    /* planes.push_back(Plane(0, 0, 1, -1)); */
    /* planes.push_back(Plane(0, 0, -1, -1)); */

    // Random points on a sphere
    std::list<Plane> sphere_planes;
    int N;
    if (argc > 1) {
        N = atoi(argv[1]);
    } else {
        N = 200;
    }

    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        sphere_planes.push_back(tangent_plane<K>(*g++));
    }

    // Compute dual convex hulls
    CGAL::convex_hull_3(planes.begin(), planes.end(), dual_poly_cube, dual_traits);
    /* CGAL::convex_hull_3(sphere_planes.begin(), sphere_planes.end(), dual_poly_sphere, dual_traits); */

    // Print the dual polyhedron in an OFF file
    convert_dual_OFF<K>("dual_cube_convex_hull.off", dual_poly_cube);
    /* convert_dual_OFF<K>("dual_sphere_convex_hull.off", dual_poly_sphere); */

    // Compute associated primal polyhedrons
    typedef CGAL::Voronoi_covariance_3::internal::Build_primal_polyhedron<K, Polyhedron_dual_3, Polyhedron_3>
        Build_primal_polyhedron;
    Build_primal_polyhedron bpp_cube(dual_poly_cube, origin);
    /* Build_primal_polyhedron bpp_sphere(dual_poly_sphere, origin); */

    poly_cube.delegate(bpp_cube);
    /* poly_sphere.delegate(bpp_sphere); */

    // Print the primal polyhedrons into an OFF file
    convertToOFF<K, Polyhedron_3>("primal_cube.off", poly_cube);
    /* convertToOFF<K, Polyhedron_3>("primal_sphere.off", poly_sphere); */

    return 0;
}

