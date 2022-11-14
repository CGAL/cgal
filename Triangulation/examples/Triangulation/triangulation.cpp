#include <CGAL/config.h>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <cassert>

typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag >  K;
typedef CGAL::Triangulation<K>                        Triangulation;

int main()
{
    const int D = 5;   // we work in Euclidean 5-space
    const int N = 100; // we will insert 100 points
    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 1
    CGAL::Random_points_in_cube_d<Triangulation::Point> rand_it(D, 1.0);
    std::vector<Triangulation::Point> points;
    std::copy_n(rand_it, N, std::back_inserter(points));

    Triangulation t(D);                      // create triangulation
    assert(t.empty());
    t.insert(points.begin(), points.end());  // compute triangulation
    assert( t.is_valid() );
    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 2
    typedef Triangulation::Face Face;
    typedef std::vector<Face> Faces;
    Faces edges;
    std::back_insert_iterator<Faces> out(edges);
    t.tds().incident_faces(t.infinite_vertex(), 1, out);
    // collect faces of dimension 1 (edges) incident to the infinite vertex
    std::cout << "There are " << edges.size()
              << " vertices on the convex hull." << std::endl;

#include "triangulation1.cpp" // See below
#include "triangulation2.cpp"

    return 0;
}
