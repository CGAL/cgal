#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/jet_smoothing_3.h>

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

    // smoothing
    std::deque<Point> output;
    const unsigned int nb_neighbors = 5;
    CGAL::jet_smoothing_3(points.begin(),points.end(),std::back_inserter(output),nb_neighbors);

    return EXIT_SUCCESS;
}

