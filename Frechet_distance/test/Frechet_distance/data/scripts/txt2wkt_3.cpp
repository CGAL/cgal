#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/WKT.h>

#include <filesystem>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Curve = std::vector<Point>;

int main(int , char* argv[])
{
  const std::filesystem::path data{std::string(argv[1])};
  for (auto const& dir_entry : std::filesystem::directory_iterator{data}){
    std::filesystem::path path = dir_entry.path();
    if(path.stem() == std::filesystem::path("dataset")){
      continue;
    }
    std::ifstream in(path.string());
    Point p;
    Curve curve;
    while(in >> p){
      curve.push_back(p);
    }
    path.replace_extension(".wkt");
    std::ofstream out(path.string());
    out.precision(17);
    CGAL::IO::write_linestring_WKT(out, curve);
    out << std::endl;
  }
  return 0;
}
