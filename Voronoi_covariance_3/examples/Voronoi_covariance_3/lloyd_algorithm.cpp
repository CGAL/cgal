#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_covariance_3/halfspaces_intersection.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Convex_hull_3.h>

#include <vector>
#include <sstream>

#include "include/off.h"

// Chrono
#include <chrono>
#include <ctime>

#define SQUARE(x) ( (x) * (x) )

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;

typedef CGAL::Convex_hull_traits_3<K> Traits;
typedef Traits::Polyhedron_3 Polyhedron;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Plane_iterator Plane_iterator;

typedef CGAL::Delaunay_triangulation_3<K> DT;
typedef DT::Vertex_handle Vertex_handle;

#include <CGAL/point_generators_3.h>

template <typename K>
typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
    typename K::Vector_3 v(p.x(), p.y(), p.z());
    v = v / sqrt(v.squared_length());
    typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);

    return plane;
}

// Centroid of a polyhedron
Point compute_centroid (Polyhedron &P) {
    typedef Polyhedron::Halfedge_around_facet_circulator Hafc;

    float volume = 0.f;
    float cx = 0.f, cy = 0.f, cz = 0.f;
    Vector ex(1, 0, 0),
           ey(0, 1, 0),
           ez(0, 0, 1);

    for (Facet_iterator fit = P.facets_begin();
         fit != P.facets_end();
         fit++) {
        Hafc h0 = fit->facet_begin(), hf = h0--, hs = hf;
        hs ++;

        while (1) {
            Point a = h0->vertex()->point(),
                  b = hf->vertex()->point(),
                  c = hs->vertex()->point();

            Vector va = CGAL::ORIGIN - a;
            Vector vb = CGAL::ORIGIN - b;
            Vector vc = CGAL::ORIGIN - c;

            // Volume
            Vector nhat = CGAL::cross_product(b - a, c - a);
            volume += nhat * va;

            // Centroid
            // X
            cx += (nhat * ex) * (
                SQUARE(ex * (va + vb)) +
                SQUARE(ex * (vb + vc)) +
                SQUARE(ex * (vc + va)) );

            // Y
            cy += (nhat * ey) * (
                SQUARE(ey * (va + vb)) +
                SQUARE(ey * (vb + vc)) +
                SQUARE(ey * (vc + va)) );

            // Z
            cz += (nhat * ez) * (
                SQUARE(ez * (va + vb)) +
                SQUARE(ez * (vb + vc)) +
                SQUARE(ez * (vc + va)) );

            if (hs == h0)
                break;

            hs++; hf++;
        }
    }

    volume /= 6;
    volume = fabs(volume);
    /* std::cout << "volume = " << volume << std::endl; */
    Vector centroid(cx, cy, cz);
    centroid = centroid / (48 * volume);
    /* std::cout << "centroid = " << centroid << std::endl; */

    return (CGAL::ORIGIN + centroid);
}

template <class PolyIterator>
void lloyd_algorithm (PolyIterator poly_begin,
                      PolyIterator poly_end,
                      std::vector<Point> & points) {
    std::list<Plane> planes;
    std::list<Point> centroids;

    // Compute Delaunay triangulation
    DT dt(points.begin(), points.end());

    for (std::vector<Point>::iterator pit = points.begin();
         pit != points.end();
         pit++) {
        planes.clear();
        // Voronoi cells
        std::list<Vertex_handle> vertices;
        Vertex_handle v = dt.nearest_vertex(*pit);
        dt.incident_vertices(v, std::back_inserter(vertices));
        for (std::list<Vertex_handle>::iterator it = vertices.begin();
             it != vertices.end();
             it++) {
            Vector p = ((*it)->point() - v->point()) / 2;
            planes.push_back (Plane(CGAL::ORIGIN + p, p));
        }

        // Add planes of the polyhedron faces
        for (PolyIterator it = poly_begin;
             it != poly_end;
             it++) {
            planes.push_back(*it);
        }

        // Intersection
        Polyhedron P;
        CGAL::Convex_hull_3::halfspaces_intersection(planes.begin(),
                                                     planes.end(),
                                                     P,
                                                     K());
        // Centroid
        Point centroid = compute_centroid(P);
        centroids.push_back(centroid);
    }

    // Replace the initial points by the computed centroids
    points.clear();
    for (std::list<Point>::iterator it = centroids.begin();
         it != centroids.end();
         it++) {
        points.push_back(*it);
    }
}

int main (int argc, char *argv[]) {
    // Cube
    std::list<Plane> planes;
    /* planes.push_back(Plane(1, 0, 0, -1)); */
    /* planes.push_back(Plane(-1, 0, 0, -1)); */
    /* planes.push_back(Plane(0, 1, 0, -1)); */
    /* planes.push_back(Plane(0, -1, 0, -1)); */
    /* planes.push_back(Plane(0, 0, 1, -1)); */
    /* planes.push_back(Plane(0, 0, -1, -1)); */

    std::vector<Point> points;
    int N, steps;
    // Number of points
    if (argc > 1) {
        N = atoi(argv[1]);
    } else {
        N = 50;
    }

    // Number of steps
    if (argc > 2) {
        steps = atoi(argv[2]);
    } else {
        steps = 10;
    }

    /* CGAL::Random_points_in_cube_3<Point> g(1); */
    CGAL::Random_points_on_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        Point p = *g++;
        points.push_back(p);
        planes.push_back(tangent_plane<K>(p));
        std::cout << p << std::endl;
    }
    Polyhedron P_before, P_after;
    CGAL::convex_hull_3(points.begin(), points.end(), P_before);
    convertToOFF<K, Polyhedron>("before.off", P_before);

    // Apply Lloyd algorithm
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_time;
    for (int i = 0; i < steps; i++) {
        std::cout << "iter " << i << std::endl;
        start = std::chrono::system_clock::now();
        lloyd_algorithm(planes.begin(),
                        planes.end(),
                        points);
        end = std::chrono::system_clock::now();

        elapsed_time = end - start;
        std::cout << "Execution time : " << elapsed_time.count() << "s\n";
    }

    for (std::vector<Point>::iterator it = points.begin();
         it != points.end();
         it++) {
        std::cout << *it << std::endl;
    }

    CGAL::convex_hull_3(points.begin(), points.end(), P_after);
    convertToOFF<K, Polyhedron>("after.off", P_after);

    return 0;
}

