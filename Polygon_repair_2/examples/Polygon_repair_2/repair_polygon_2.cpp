#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair_2 = CGAL::Polygon_repair_2::Polygon_repair_2<Kernel>;

int main(int argc, char* argv[]) {

  for (const auto& file: std::__fs::filesystem::directory_iterator("../../test/Polygon_repair_2/data/in")) {
    std::cout << "Reading " << file.path().filename() << "..." << std::endl;

    // if (file.path().filename() != "hole-outside.wkt") continue;

    std::string in;
    std::getline(std::ifstream(file.path()), in);
    std::istringstream iss(in);
    std::size_t is_polygon = in.find("POLYGON");
    std::size_t is_multipolygon = in.find("MULTIPOLYGON");
    Multipolygon_with_holes_2 rmp;
    if (is_polygon == 0) {
      Polygon_with_holes_2 p;
      CGAL::IO::read_polygon_WKT(iss, p);
      CGAL::draw(p);
      rmp = CGAL::Polygon_repair_2::repair(p);
    } else if (is_multipolygon == 0) {
      Multipolygon_with_holes_2 mp;
      CGAL::IO::read_multi_polygon_WKT(iss, mp);
      CGAL::draw(mp);
      rmp = CGAL::Polygon_repair_2::repair(mp);
    } std::ostringstream oss;
    CGAL::IO::write_multi_polygon_WKT(oss, rmp);
    std::string out = oss.str();
    std::cout << "\tin:  " << in << std::endl;
    std::cout << "\tout: " << out;
    CGAL::draw(rmp);

  }

  return 0;
}
