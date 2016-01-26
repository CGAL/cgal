#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/vcm_estimate_normals.h>
#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_3<Kernel> Point_3;
typedef CGAL::Vector_3<Kernel> Vector_3;
typedef std::pair<Point_3, Vector_3> PointWithNormal;

int main (void) {
    // Generate points on a plane
    int k = 100;
    double r = 10;

    std::vector<PointWithNormal> points;
    points.push_back(std::make_pair(Point_3(0, 0, 0), Vector_3(0, 0, 0)));
    for (int i = 0; i < k; ++i) {
        double theta = 2 * i * CGAL_PI / k;
        points.push_back(std::make_pair(Point_3(r * cos(theta), r * sin(theta), 0),
                                        Vector_3(0, 0, 0)));
    }

    // Estimate the normals using VCM
    double R = 20;
    vcm_estimate_normals(points.begin(), points.end(),
                         CGAL::First_of_pair_property_map<PointWithNormal>(),
                         CGAL::Second_of_pair_property_map<PointWithNormal>(),
                         R, 0.0);

    std::cout << "Normal is " << points[0].second << std::endl;

    // The normal at the origin should be (0, 0, 1)
    double epsilon=2e-5;
    assert(points[0].second.x() < epsilon && points[0].second.x() > -epsilon);
    assert(points[0].second.y() < epsilon && points[0].second.y() > -epsilon);
    assert(points[0].second.z() < 0 || (points[0].second.z() < 1+epsilon && points[0].second.z() > 1-epsilon));
    assert(points[0].second.z() > 0 || (points[0].second.z() < -1+epsilon && points[0].second.z() > -1-epsilon));

    return 0;
}

