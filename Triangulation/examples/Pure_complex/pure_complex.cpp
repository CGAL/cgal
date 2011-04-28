#include <CGAL/point_generators_d.h>
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Pure_complex.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <iterator>
#include <iostream>
#include <vector>

int main()
{
    const int D = 5;   // we work in euclidean 5-space
    const int N = 100; // we will insert 100 points

    // |Cartesian_d| is a model of the concept PureComplexTraits
    typedef CGAL::Cartesian_d<double> K;

    // |Filtered_kernel_d|  provides exact geometric predicates
    typedef CGAL::Filtered_kernel_d<K> FK;

    // Here is our Pure_complex type:
    typedef CGAL::Pure_complex<FK> PC;

    typedef PC::Point_d Point;

    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 1

    // Instanciate a random point generator
    CGAL::Random rng;
    typedef CGAL::Random_points_in_iso_box_d<Point> Random_points_iterator;
    Random_points_iterator rand_it(D, 1.0, rng);

    // Generate N random points
    std::vector<Point> points;
    CGAL::copy_n(rand_it, N, std::back_inserter(points));

    PC pc(D);
    assert(pc.empty());

    // insert the points in the pure complex
    pc.insert(points.begin(), points.end());
    assert( pc.is_valid() );

    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 2

    typedef PC::Face Face;
    typedef std::vector<Face> Faces;
    Faces edges;
    std::back_insert_iterator<Faces> out(edges);
    pc.gather_incident_upper_faces(pc.infinite_vertex(), 1, out);
    // Count the number of points on the convex hull
    std::cout << "\nThere are " << edges.size() << " vertices on the convex hull.";
    edges.clear();

    // cleanup
    pc.clear();
    assert(pc.empty());

    std::cout << std::endl;
    return 0;
}
