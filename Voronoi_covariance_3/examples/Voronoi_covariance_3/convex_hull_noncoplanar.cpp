#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>

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

template <typename K, class Polyhedron_dual, class Polyhedron>
class Build_primal_polyhedron :
    public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> {
    typedef typename Polyhedron::HalfedgeDS HDS;
    const Polyhedron_dual & _dual;
    typedef typename K::Point_3 Primal_point_3;
    Primal_point_3 origin;

    public:
    Build_primal_polyhedron (const Polyhedron_dual & dual, Primal_point_3 const& o = Primal_point_3(0, 0, 0)) : _dual (dual), origin(o)
        {}

    // Compute the primal point associated to a triple of dual planes

    void operator () (HDS &hds)
    {
        typedef typename K::RT RT;
        typedef typename K::Point_3 Point_3;

        // Typedefs for dual
        typedef typename Polyhedron_dual::Facet Facet;
        typedef typename Polyhedron_dual::Vertex Vertex;
        typedef typename Vertex::Point_3 Plane;
        typedef typename Polyhedron_dual::Facet_const_handle
            Facet_const_handle;
        typedef typename Polyhedron_dual::Facet_const_iterator
            Facet_const_iterator;
        typedef typename Polyhedron_dual::Vertex_const_iterator
            Vertex_const_iterator;

        // Typedefs for primal
        typename CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        B.begin_surface(_dual.size_of_facets(),
                        _dual.size_of_vertices(),
                        _dual.size_of_vertices());

        std::map <Facet_const_handle, size_t> primal_vertices;
        size_t n = 0;

        std::cout << "Origin : " << origin << std::endl;

        // First, computing the primal vertices
        for (Facet_const_iterator it = _dual.facets_begin();
             it != _dual.facets_end(); ++it, ++n) {
            typename Facet::Halfedge_const_handle h = it->halfedge();
            // Build the dual plane corresponding to the current facet
            Plane p1 = h->vertex()->point();
            Plane p2 = h->next()->vertex()->point();
            Plane p3 = h->next()->next()->vertex()->point();

            RT dp1 = p1.d() + origin.x() * p1.a()
                + origin.y() * p1.b() + origin.z() * p1.c();
            RT dp2 = p2.d() + origin.x() * p2.a()
                + origin.y() * p2.b() + origin.z() * p2.c();
            RT dp3 = p3.d() + origin.x() * p3.a()
                + origin.y() * p3.b() + origin.z() * p3.c();

            // Normal to the dual plane
            RT alpha = (dp1 * p2.b() - dp2 * p1.b()) *
                (dp1 * p3.c() - dp3 * p1.c()) -
                (dp1 * p2.c() - dp2 * p1.c()) *
                (dp1 * p3.b() - dp3 * p1.b());

            RT beta  = (dp1 * p2.c() - dp2 * p1.c()) *
                (dp1 * p3.a() - dp3 * p1.a()) -
                (dp1 * p2.a() - dp2 * p1.a()) *
                (dp1 * p3.c() - dp3 * p1.c());

            RT gamma = (dp1 * p2.a() - dp2 * p1.a()) *
                (dp1 * p3.b() - dp3 * p1.b()) -
                (dp1 * p2.b() - dp2 * p1.b()) *
                (dp1 * p3.a() - dp3 * p1.a());

            // last coefficient of the dual plane equation
            RT d = (alpha * p1.a() + beta * p1.b() + gamma * p1.c()) / dp1;

            std::cout << "p1.a : " << p1.a() << std::endl;
            std::cout << "alpha : " << alpha << std::endl;
            std::cout << "p1.b : " << p1.b() << std::endl;
            std::cout << "beta : " << beta << std::endl;
            std::cout << "p1.c : " << p1.c() << std::endl;
            std::cout << "gamma : " << gamma << std::endl;
            std::cout << "d : " << d << std::endl;

            // Primal vertex associated to the current dual plane
            // TODO: add origin
            Point_3 p(origin.x() - alpha / d, origin.y() - beta / d, origin.z() - gamma / d);

            std::cout << p << std::endl;

            B.add_vertex(p);
            primal_vertices[it] = n;
        }

        // Then, add facets to the primal polyhedron
        // To do this, for each dual vertex, we circulate around this vertex
        // and we add an edge between each facet we encounter
        for (Vertex_const_iterator it = _dual.vertices_begin();
             it != _dual.vertices_end(); ++it, ++n) {
            typename Polyhedron_dual::Halfedge_around_vertex_const_circulator
                h0 = it->vertex_begin(), hf = h0;
            B.begin_facet();
            do {
                B.add_vertex_to_facet(primal_vertices[hf->facet()]);
            } while (++hf != h0);
            B.end_facet();
        }

        B.end_surface();
    }
    };

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
    typedef Build_primal_polyhedron<K, Polyhedron_dual_3, Polyhedron_3>
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

