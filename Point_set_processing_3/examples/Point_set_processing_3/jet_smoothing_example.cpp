#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/jet_smooth_point_set.h>
#include <deque>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
    // generate point set
    std::deque<Point> points;
    points.push_back(Point( 0.0, 0.0, 0.01));
    points.push_back(Point(-0.1,-0.1, 0.02));
    points.push_back(Point(-0.1,-0.2, 0.01));
    points.push_back(Point(-0.1, 0.1, 0.02));
    points.push_back(Point( 0.1,-0.1, 0.00));
    points.push_back(Point( 0.1, 0.2, 0.01));
    points.push_back(Point( 0.2, 0.0, 0.02));
    points.push_back(Point( 0.2, 0.1, 0.00));
    points.push_back(Point( 0.0,-0.1, 0.01));

    // smoothing
    const unsigned int nb_neighbors = 8;
    CGAL::jet_smooth_point_set(points.begin(),points.end(),nb_neighbors);

    return EXIT_SUCCESS;
}

