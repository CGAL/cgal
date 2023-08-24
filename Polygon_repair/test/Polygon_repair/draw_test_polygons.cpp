#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/Polygon_repair.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/Polygon_repair/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main(int argc, char* argv[]) {

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
      rmp = CGAL::Polygon_repair::repair_odd_even(p);
    } else if (in.find("MULTIPOLYGON") == 0) {
      Multipolygon_with_holes_2 mp;
      CGAL::IO::read_multi_polygon_WKT(iss, mp);
      CGAL::draw(mp);
      rmp = CGAL::Polygon_repair::repair_odd_even(mp);
    } std::ostringstream oss;
    CGAL::IO::write_multi_polygon_WKT(oss, rmp);
    std::string out = oss.str();
    std::cout << "\tin:  " << in << std::endl;
    std::cout << "\tout: " << out;
    CGAL::draw(rmp);

  }

  return 0;
}
