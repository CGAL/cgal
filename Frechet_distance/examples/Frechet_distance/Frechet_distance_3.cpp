#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_3.h>
#include <CGAL/Simple_cartesian.h>

#include <ostream>

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Frechet_distance_traits_3<Kernel>;
using Point = Traits::Point_d;

int main(int argc, char* argv[])
{
    std::vector<Point> A, B;
    std::pair<double, double> res = CGAL::approximate_Frechet_distance<Traits>(A, B, 0.000001);
    std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    return 0;
}
