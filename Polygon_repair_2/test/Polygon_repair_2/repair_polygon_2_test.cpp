#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair_2 = CGAL::Polygon_repair_2::Polygon_repair_2<Kernel>;

int main(int argc, char* argv[]) {

  for (const auto& file: std::filesystem::directory_iterator("data/in")) {
    if (file.path().filename().extension() != ".wkt") continue;
    std::cout << "Testing " << file.path().filename() << "... ";

    // Read test file
    std::string in;
    std::getline(std::ifstream(file.path()), in);

    // Load test file and repair to create output
    std::istringstream iss(in);
    Multipolygon_with_holes_2 rmp;
    if (in.find("POLYGON") == 0) {
      Polygon_with_holes_2 p;
      if (in != "POLYGON()") { // maybe should be checked in WKT reader
        CGAL::IO::read_polygon_WKT(iss, p);
      } rmp = CGAL::Polygon_repair_2::repair_odd_even(p);
    } else if (in.find("MULTIPOLYGON") == 0) {
      Multipolygon_with_holes_2 mp;
      CGAL::IO::read_multi_polygon_WKT(iss, mp);
      rmp = CGAL::Polygon_repair_2::repair_odd_even(mp);
    } std::ostringstream oss;
    CGAL::IO::write_multi_polygon_WKT(oss, rmp);
    std::string out = oss.str();

    // Read reference file
    std::string ref_path = "data/ref/";
    ref_path += file.path().filename();
    std::ifstream ref_ifs(ref_path);
    if (ref_ifs.fail()) {
      std::cout << std::endl << "\tin:  " << in << std::endl;
      std::cout << "\tout: " << out;
      std::cout << "\tno reference output -> skipped" << std::endl;
      continue;
    } std::string ref;
    std::getline(ref_ifs, ref);
    ref += "\n";

    // Compare output with reference file
    if (ref == out) {
      std::cout << "ok" << std::endl;
    } else {
      std::cout << "fail" << std::endl;
      std::cout << "\tin:  " << in << std::endl;
      std::cout << "\tout: " << out;
      std::cout << "\tref: " << ref;
    } CGAL_assertion(ref == out);

    // Test orientations
    for (auto const& polygon: rmp.polygons()) {
      CGAL_assertion(polygon.outer_boundary().orientation() == CGAL::COUNTERCLOCKWISE);
      for (auto const &hole: polygon.holes()) {
        CGAL_assertion(hole.orientation() == CGAL::CLOCKWISE);
      }
    }
  }

  return 0;
}
