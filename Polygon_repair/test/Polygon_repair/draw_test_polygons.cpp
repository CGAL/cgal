#include <iostream>
#include <fstream>
#include <sstream>

// work around for old compilers (Apple clang < 11 for example)
#define HAS_FILESYSTEM 1
#if defined(__has_include)
#if !__has_include(<filesystem>)
#undef HAS_FILESYSTEM
#define HAS_FILESYSTEM 0
#endif
#endif


#if HAS_FILESYSTEM

#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main() {

  for (const auto& file: std::filesystem::directory_iterator("data/in")) {
    if (file.path().filename().extension() != ".wkt") continue;
    std::cout << "Reading " << file.path().filename() << "..." << std::endl;

    if (file.path().filename() != "nesting-spike.wkt") continue;

    std::string in;
    std::getline(std::ifstream(file.path()), in);
    std::istringstream iss(in);
    Multipolygon_with_holes_2 rmp;

    if (in.find("POLYGON") == 0) {
      Polygon_with_holes_2 p;
      if (in != "POLYGON()") { // maybe should be checked in WKT reader
        CGAL::IO::read_polygon_WKT(iss, p);
      } CGAL::draw(p);
      rmp = CGAL::Polygon_repair::repair(p, CGAL::Polygon_repair::Even_odd_rule());
    } else if (in.find("MULTIPOLYGON") == 0) {
      Multipolygon_with_holes_2 mp;
      CGAL::IO::read_multi_polygon_WKT(iss, mp);
      CGAL::draw(mp);
      rmp = CGAL::Polygon_repair::repair(mp, CGAL::Polygon_repair::Even_odd_rule());
    } std::ostringstream oss;
    CGAL::IO::write_multi_polygon_WKT(oss, rmp);
    std::string out = oss.str();
    std::cout << "\tin:  " << in << std::endl;
    std::cout << "\tout: " << out;
    CGAL::draw(rmp);

  }

  return 0;
}

#else

int main()
{
  std::cout << "Warning: filesystem feature is not present on the system, nothing will be tested\n";
  return 0;
}

#endif
