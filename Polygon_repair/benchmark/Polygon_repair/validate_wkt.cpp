#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main() {

//  std::string folder = "/Users/ken/Downloads/big polygons/";
  std::string folder = "data/in";

  for (const auto& file: std::filesystem::directory_iterator(folder)) {
    if (file.path().filename().extension() != ".wkt") continue;
    std::cout << "Testing " << file.path().filename() << "... ";

    // Read test file
    std::string in;
    std::getline(std::ifstream(file.path()), in);

    // Load test file
    std::istringstream iss(in);
    bool valid = true;
    if (in.find("POLYGON") == 0) {
      Polygon_with_holes_2 p;
      if (in != "POLYGON()") { // maybe should be checked in WKT reader
        CGAL::IO::read_polygon_WKT(iss, p);
        valid = CGAL::Polygon_repair::is_valid(p);
      }
    } else if (in.find("MULTIPOLYGON") == 0) {
      Multipolygon_with_holes_2 mp;
      CGAL::IO::read_multi_polygon_WKT(iss, mp);
      valid = CGAL::Polygon_repair::is_valid(mp);
    } if (valid) std::cout << "Valid" << std::endl;
  }

  return 0;
}
