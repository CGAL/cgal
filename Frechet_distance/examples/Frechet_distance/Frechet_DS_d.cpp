#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_d.h>
#include <CGAL/Frechet_distance/Neighbor_search.h>
#include <CGAL/Epick_d.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>
#include <filesystem>

using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>;
using Traits = CGAL::Frechet_distance_traits_d<Kernel>;
using Point = Traits::Point_d;
using Curve = std::vector<Point>;
using Curves = std::vector<Curve>;

int main()
{
    Curves curves;
    Curve polylineA = { Point(0,0,0,0), Point(1,0,0,0), Point(1,1,0,1),Point(1,1,1,0)};
    Curve polylineB = { Point(0,0,0,0), Point(1,0,0,0), Point(1,1,0,0),Point(1,1,1,0)};
    Curve polylineC = { Point(0,0,0,1), Point(1,0,0,0), Point(1,1,1,0),Point(1,1,0,0)};
    curves.push_back(polylineA);
    curves.push_back(polylineB);
    curves.push_back(polylineC);

    // last curve is the query curve
    Curve query = curves.back();
    curves.pop_back();

    CGAL::Frechet_distance::Neighbor_search<Curve,Traits> ds(curves);

    for(const Curve& c : curves){
        std::pair<double, double> res = CGAL::bounded_error_Frechet_distance(c, query, 0.000001);
        std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    }
    double distance = 1.1;
    std::vector<std::size_t> result = ds.get_close_curves(query, distance);

    std::cout << result.size() << " curves at Frechet distance closer than " << distance << std::endl;
}
