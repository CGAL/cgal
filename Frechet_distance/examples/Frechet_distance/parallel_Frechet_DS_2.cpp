#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_2.h>
#include <CGAL/Frechet_distance/Neighbor_search.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>
#include <filesystem>

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
using Point = Traits::Point_d;
using Curve = std::vector<Point>;
using Curves = std::vector<Curve>;

int main()
{
    Curves curves;

    const std::filesystem::path data{"./data_2d"};
    std::vector<std::string> filenames;
    for (auto const& dir_entry : std::filesystem::directory_iterator{data}) {
        filenames.push_back(dir_entry.path().string());
    }
    std::sort(filenames.begin(), filenames.end());

    for (auto const& filename : filenames) {
        std::cout << filename << std::endl;
        std::ifstream in(filename);
        curves.push_back(Curve());
        CGAL::IO::read_linestring_WKT(in, curves.back());
    }

    // last curve is the query curve
    Curve query = curves.back();
    curves.pop_back();

    using Neighbor_search = CGAL::Frechet_distance::Neighbor_search<Curve, Traits>;
    Neighbor_search ds(curves, CGAL::Parallel_if_available_tag());

    for(const Curve& c : curves){
        std::pair<double, double> res = CGAL::bounded_error_Frechet_distance(c, query, 0.000001);
        std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    }
    double distance = 16;
    std::vector<std::size_t> result = ds.get_close_curves<CGAL::Parallel_if_available_tag>(query, distance);

    std::cout << result.size() << " curves at Frechet distance closer than " << distance << std::endl;
}
