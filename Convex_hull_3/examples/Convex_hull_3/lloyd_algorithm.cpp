#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;

typedef CGAL::Convex_hull_traits_3<K> Traits;
typedef Traits::Polygon_mesh Polyhedron;

typedef CGAL::Delaunay_triangulation_3<K> DT;
typedef DT::Vertex_handle Vertex_handle;

// Function object that computes the volume and the centroid of a polyhedron.
template <class K>
class Centroid_volume_accumulator {
    public:
        typedef typename K::Point_3 Point;
        typedef typename K::Vector_3 Vector;

        Centroid_volume_accumulator() : vol(0),
                                        cx(0), cy(0), cz(0) {}

        void operator () (const Point &a,
                          const Point &b,
                          const Point &c) {
            Vector ex(1, 0, 0),
                   ey(0, 1, 0),
                   ez(0, 0, 1);

            Vector va = CGAL::ORIGIN - a;
            Vector vb = CGAL::ORIGIN - b;
            Vector vc = CGAL::ORIGIN - c;

            // Updating the volume...
            Vector nhat = CGAL::cross_product(b - a, c - a);
            vol += nhat * va;

            // ... and the centroid
            // X
            cx += (nhat * ex) * (
                CGAL::square(ex * (va + vb)) +
                CGAL::square(ex * (vb + vc)) +
                CGAL::square(ex * (vc + va)) );

            // Y
            cy += (nhat * ey) * (
                CGAL::square(ey * (va + vb)) +
                CGAL::square(ey * (vb + vc)) +
                CGAL::square(ey * (vc + va)) );

            // Z
            cz += (nhat * ez) * (
                CGAL::square(ez * (va + vb)) +
                CGAL::square(ez * (vb + vc)) +
                CGAL::square(ez * (vc + va)) );
        }

        void end () {
            vol /= 6;
            cx /= (48 * vol);
            cy /= (48 * vol);
            cz /= (48 * vol);
        }

        Point centroid () const {
            return Point(cx, cy, cz);
        }

        typename K::FT volume () const {
            return vol;
        }

        void reset () {
            vol = 0;
            cx = 0;
            cy = 0;
            cz = 0;
        }

    private:
        // Volume
        typename K::FT vol;

        // Centroid
        typename K::FT cx, cy, cz;
};

// Apply a function object to all the triangles composing the faces of a polyhedron.
template <typename Polyhedron, class F>
F& apply_function_object_polyhedron (Polyhedron &P,
                                     F &f) {
    typedef typename Polyhedron::Halfedge_around_facet_circulator Hafc;
    typedef typename Polyhedron::Facet_iterator Facet_iterator;

    f.reset();

    for (Facet_iterator fit = P.facets_begin();
         fit != P.facets_end();
         fit++) {
        Hafc h0 = fit->facet_begin(), hf = h0--, hs = hf;
        hs ++;

        while (1) {
            // Apply 'f' on each triangle of the polyhedron's facet
            f ( h0->vertex()->point(),
                hf->vertex()->point(),
                hs->vertex()->point() );

            if (hs == h0)
                break;

            hs++; hf++;
        }
    }

    f.end();

    return f;
}

// Lloyd algorithm
// Generate points uniformly sampled inside a polyhedron.
// An initial set of points needs to be given.
template <class PolyIterator>
void lloyd_algorithm (PolyIterator poly_begin,
                      PolyIterator poly_end,
                      std::vector<Point>& points) {
    std::list<Plane> planes;
    std::list<Point> centroids;
    Centroid_volume_accumulator<K> centroid_acc;

    // Compute Delaunay triangulation
    DT dt(points.begin(), points.end());

    for (DT::Finite_vertices_iterator vit = dt.finite_vertices_begin();
         vit != dt.finite_vertices_end();
         ++vit) {
        planes.clear();
        // Voronoi cells
        std::list<Vertex_handle> vertices;
        dt.incident_vertices(vit, std::back_inserter(vertices));
        for (std::list<Vertex_handle>::iterator it = vertices.begin();
             it != vertices.end();
             it++)
        {
            if (dt.is_infinite(*it)) continue;
            Vector p(vit->point(),(*it)->point());
            planes.push_back (Plane(CGAL::midpoint((*it)->point(), vit->point()), p));
        }

        // Add planes of the polyhedron faces
        for (PolyIterator it = poly_begin;
             it != poly_end;
             it++) {
            planes.push_back(*it);
        }

        // Intersection
        Polyhedron P;
        CGAL::halfspace_intersection_3(planes.begin(),
                                       planes.end(),
                                       P,
                                       boost::make_optional(vit->point()));

        // Centroid
        apply_function_object_polyhedron(P, centroid_acc);
        centroids.push_back(centroid_acc.centroid());
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
    planes.push_back(Plane(1, 0, 0, -1));
    planes.push_back(Plane(-1, 0, 0, -1));
    planes.push_back(Plane(0, 1, 0, -1));
    planes.push_back(Plane(0, -1, 0, -1));
    planes.push_back(Plane(0, 0, 1, -1));
    planes.push_back(Plane(0, 0, -1, -1));

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

    CGAL::Random_points_in_sphere_3<Point> g;
    for (int i = 0; i < N; i++) {
        Point p = *g++;
        points.push_back(p);
    }

    std::ofstream bos("before_lloyd.xyz");
    std::copy(points.begin(), points.end(),
              std::ostream_iterator<Point>(bos, "\n"));

    // Apply Lloyd algorithm: will generate points
    // uniformly sampled inside a cube.
    for (int i = 0; i < steps; i++) {
        std::cout << "iteration " << i + 1 << std::endl;

        CGAL::Timer timer;
        timer.start();
        lloyd_algorithm(planes.begin(),
                        planes.end(),
                        points);
        timer.stop();

        std::cout << "Execution time : " << timer.time() << "s\n";
    }

    std::ofstream aos("after_lloyd.xyz");
    std::copy(points.begin(), points.end(),
              std::ostream_iterator<Point>(aos, "\n"));

    return 0;
}

