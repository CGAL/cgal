#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <iterator>
#include <iostream>
#include <vector>

int main()
{
    const int D = 5;   // we work in euclidean 5-space
    const int N = 100; // we will insert 100 points

    // |Cartesian_d| is a model of the concept TriangulationTraits
    typedef CGAL::Cartesian_d<double> K;

    // |Filtered_kernel_d|  provides exact geometric predicates
    typedef CGAL::Filtered_kernel_d<K> FK;

    // Here is our Triangulation type:
    typedef CGAL::Triangulation<FK> T;

    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 1

    // Instanciate a random point generator
    CGAL::Random rng;
    typedef CGAL::Random_points_in_cube_d<T::Point> Random_points_iterator;
    Random_points_iterator rand_it(D, 1.0, rng);

    // Generate N random points
    std::vector<T::Point> points;
    CGAL::copy_n(rand_it, N, std::back_inserter(points));

    T t(D);
    assert(t.empty());

    // insert the points in the triangulation
    t.insert(points.begin(), points.end());
    assert( t.is_valid() );

    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 2

    typedef T::Face Face;
    typedef std::vector<Face> Faces;
    Faces edges;
    std::back_insert_iterator<Faces> out(edges);
    t.incident_faces(t.infinite_vertex(), 1, out);
    // Count the number of points on the convex hull
    std::cout << "There are " << edges.size() << " vertices on the (triangulated) convex hull.";
    edges.clear();

    // cleanup
    t.clear();
    assert(t.empty());

    std::cout << std::endl;
    return 0;
}
