#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
typedef K::Point_2                                                                                                         Point_2;
typedef K::Triangle_2                                                                                                 Triangle_2;
typedef std::vector<Point_2>                                                                                Container;
typedef CGAL::Random_points_in_triangle_2<Point_2>                                         Point_generator;

int main() {
        std::cout << "Creating 100 random points in a triangle in 2D." << std::endl;

        // The input triangle is as follows
        Triangle_2 tri(Point_2(0,0),Point_2(1,0),Point_2(0,1));

        // generated points are in that container
        Container points;

        // creating the generator, input is the Triangle_2 tri
        Point_generator g(tri);

        // get 100 random points in tri
        std::copy_n(g, 100, std::back_inserter(points));

        // Check that we have really created 100 points.
        assert( points.size() == 100);

        // print the first point that was generated
        std::cout << points[0] << std::endl;

        return 0;
}

