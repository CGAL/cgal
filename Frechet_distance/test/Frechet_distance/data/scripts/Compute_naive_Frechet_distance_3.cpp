#include <CGAL/Frechet_naive.h>
#include <CGAL/Frechet_distance_traits_3.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/number_utils.h>
#include <iostream>
#include <sstream>
#include <fstream>

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Frechet_distance_traits_3<Kernel>;
using Point = Traits::Point;
using Points = std::vector<Point>;
using Curve = CGAL::Frechet_distance_::internal::Curve<Traits>;

void readCurve(std::ifstream& curve_file, Points& points)
{
    // Read everything into a stringstream.
    std::stringstream ss;
    ss << curve_file.rdbuf();
    CGAL::set_ascii_mode(ss);

    Point p;
    auto ignore_count = std::numeric_limits<std::streamsize>::max();
    while (ss >> p) {
        ss.ignore(ignore_count, '\n');

        if ((!points.empty()) && (p == points.back())) {
            continue;
        }
        points.push_back(p);
    }
}

int main(int argc, char* argv[])
{
    // TODO AF: why do we need this here?
    CGAL::Protect_FPU_rounding<true> p;

    assert(argc == 3);
    double epsilon = 10e-10;

    std::ifstream file1(argv[1]);
    std::ifstream file2(argv[2]);

    Points points1, points2;
    readCurve(file1, points1);
    readCurve(file2, points2);
    auto curve1 = Curve(points1);
    auto curve2 = Curve(points2);

    CGAL::Frechet_distance_::internal::FrechetNaive<Curve> frechet_naive;
    auto dist = frechet_naive.calcDistance(curve1, curve2, epsilon);
    std::cout.precision(17);
    std::cout << (dist.second - dist.first)/2. << std::endl;
}
