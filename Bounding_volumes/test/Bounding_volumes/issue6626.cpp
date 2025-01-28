#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>
#include "stdio.h"

typedef CGAL::Cartesian_d<double> Kernel;
typedef CGAL::MP_Float ET;
typedef CGAL::Approximate_min_ellipsoid_d_traits_d<Kernel, ET> Traits;
typedef Traits::Point Point;
typedef CGAL::Approximate_min_ellipsoid_d<Traits> Ellipsoid;

// Cube of non-zero volume (8 corner points).
std::vector<Point> cube() {
        // 8 corner points.
        double l[3] = { 390,  971, 8957 }; // lower (0)
        double u[3] = { 490, 1071, 9057 }; // upper (1)
        double p[8][3] = {
                { l[0], l[1], l[2] }, // (0, 0, 0)
                { u[0], l[1], l[2] }, // (1, 0, 0)
                { u[0], u[1], l[2] }, // (1, 1, 0)
                { l[0], u[1], l[2] }, // (0, 1, 0)
                { l[0], l[1], u[2] }, // (0, 0, 1)
                { u[0], l[1], u[2] }, // (1, 0, 1)
                { u[0], u[1], u[2] }, // (1, 1, 1)
                { l[0], u[1], u[2] }, // (0, 1, 1)
        };
        std::vector<Point> points;
        for (size_t i = 0; i < 8; ++i)
                points.push_back(Point(3, std::begin(p[i]), std::end(p[i])));
        // Center of cube.
        double center[3] = { 0, 0, 0 };
        for (size_t i = 0; i < points.size(); ++i) {
                center[0] += points[i][0];
                center[1] += points[i][1];
                center[2] += points[i][2];
        }
        center[0] /= points.size();
        center[1] /= points.size();
        center[2] /= points.size();
        // Point-center distances, their average and standard deviation.
        double distances[8];
        double distances_avg = 0;
        for (size_t i = 0; i < points.size(); ++i) {
                distances[i] = sqrt(
                        pow(points[i][0] - center[0], 2) +
                        pow(points[i][1] - center[1], 2) +
                        pow(points[i][2] - center[2], 2)
                );
                distances_avg += distances[i];
        }
        distances_avg /= points.size();
        // Assert that cube has non-zero volume.
        assert( (distances_avg > 86.60) && (distances_avg < 86.61) );
        double distances_err = 0;
        for (size_t i = 0; i < points.size(); ++i) {
                distances_err = pow(distances[i] - distances_avg, 2);
        }
        distances_err /= points.size();
        distances_err = sqrt(distances_err);
        // Assert that points describe cube.
        assert(distances_err <= 0.00001);
        return points;
}

// Assertion fails.
int main() {
        // Cube of non-zero volume (8 corner points).
        auto points = cube();
        // Ellipsoid enclosing cube.
        Ellipsoid ellipsoid(0.01, points.begin(), points.end(), Traits());
        // Assert that ellipsoid enclosing non-zero volume cube has non-zero volume.
        assert(ellipsoid.is_full_dimensional());
        // Access center.
        for(auto it = ellipsoid.center_cartesian_begin(); it != ellipsoid.center_cartesian_end(); ++it){
          std::cout << *it << " ";
        }
        std::cout << std::endl;

        // Access axes.
        ellipsoid.axes_lengths_begin();
        for(auto it = ellipsoid.axes_lengths_begin(); it != ellipsoid.axes_lengths_end(); ++it){
          std::cout << *it << " ";
        }
        std::cout << std::endl;

        return 0;
}
