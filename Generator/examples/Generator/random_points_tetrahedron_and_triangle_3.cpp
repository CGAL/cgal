#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
typedef K::Point_3                                                                                                         Point_3;
typedef K::Triangle_3                                                                                                 Triangle_3;
typedef K::Tetrahedron_3                                                                                         Tetrahedron_3;
typedef CGAL::Random_points_in_triangle_3<Point_3>                                         Point_generator_i;
typedef CGAL::Random_points_in_tetrahedron_3<Point_3>                                 Point_generator_ii;

int main() {
        std::cout << "This example does two things:" << std::endl;
        std::cout << "  (i) it creates 100 random points in a triangle in 3D; and" << std::endl;
        std::cout << "  (ii) it creates 100 random points in a tetrahedron in 3D." << std::endl;

        // The input triangle is as follows
        Triangle_3 tri(Point_3(0,0,0),Point_3(1,0,0),Point_3(0,1,0));
        Tetrahedron_3 tet(Point_3(0,0,0),Point_3(1,0,0),Point_3(0,1,0),Point_3(0,0,1));

        // we get output points in these containers
        std::vector<Point_3> points_in_tri, points_in_tet;

        // creating the first generator, input is the Triangle_3 tri
        Point_generator_i g_i(tri);

        // creating the second generator, input is the Tetrahedron_3 tet
        Point_generator_ii g_ii(tet);

        // get 100 random points in tri
        std::copy_n(g_i, 100, std::back_inserter(points_in_tri));

        // get 100 random points in tet
        std::copy_n(g_ii, 100, std::back_inserter(points_in_tet));

        // Check that we have really created 100 points.
        assert( points_in_tri.size() == 100);

        // Check that we have really created 100 points.
        assert( points_in_tet.size() == 100);

        // print the first points
        std::cout << "In triangle: " << points_in_tri[0] << std::endl;
        std::cout << "In tetrahedron: " << points_in_tet[0] << std::endl;

        return 0;
}

