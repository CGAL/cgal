#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>

// work around for old compilers (Apple clang < 11 for example)
#define HAS_FILESYSTEM 1
#if defined(__has_include)
#if !__has_include(<filesystem>)
#undef HAS_FILESYSTEM
#define HAS_FILESYSTEM 0
#endif
#endif


#if HAS_FILESYSTEM

#include <fstream>
#include <sstream>
#include <filesystem>
#include <cassert>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main() {

  for (const auto& file: std::filesystem::directory_iterator("data/in")) {
    if (file.path().filename().extension() != ".wkt") continue;
    std::cout << "Testing " << file.path().filename() << "... ";

    // Read test file
    std::string in;
    std::getline(std::ifstream(file.path()), in);

    // Load test file and repair to create output
    std::istringstream iss(in);
    Multipolygon_with_holes_2 rmp, refmp;
    if (in.find("POLYGON") == 0) {
      Polygon_with_holes_2 p;
      if (in != "POLYGON()") { // maybe should be checked in WKT reader
        CGAL::IO::read_polygon_WKT(iss, p);
      } rmp = CGAL::Polygon_repair::repair(p, CGAL::Polygon_repair::Even_odd_rule());
    } else if (in.find("MULTIPOLYGON") == 0) {
      Multipolygon_with_holes_2 mp;
      CGAL::IO::read_multi_polygon_WKT(iss, mp);
      rmp = CGAL::Polygon_repair::repair(mp, CGAL::Polygon_repair::Even_odd_rule());
    } std::stringstream oss;
    CGAL::IO::write_multi_polygon_WKT(oss, rmp);
    std::string out = oss.str();
    rmp.clear();
    CGAL::IO::read_multi_polygon_WKT(oss, rmp);

    // Read reference file
    std::string ref_path = "data/ref/";
    ref_path += file.path().filename().string();
    std::ifstream ref_ifs(ref_path);
    if (ref_ifs.fail()) {
      std::cout << std::endl << "\tin:  " << in << std::endl;
      std::cout << "\tout: " << out;
      std::cout << "\tno reference output -> skipped" << std::endl;
      continue;
    } std::string ref;
    std::getline(ref_ifs, ref);
    ref += "\n";
    std::stringstream refss(ref);
    CGAL::IO::read_multi_polygon_WKT(refss, refmp);

    // Compare output with reference file
    if (rmp == refmp) {
      std::cout << "ok" << std::endl;
    } else {
      std::cout << "fail" << std::endl;
      std::cout << "\tin:  " << in << std::endl;
      std::cout << "\tout: " << out << std::flush;
      std::cout << "\tref: " << ref << std::flush;
    }
    assert(rmp == refmp);

    // Test orientations
    for (auto const& polygon: rmp.polygons_with_holes()) {
      assert(polygon.outer_boundary().orientation() == CGAL::COUNTERCLOCKWISE);
      for (auto const &hole: polygon.holes()) {
        assert(hole.orientation() == CGAL::CLOCKWISE);
      }
    }
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
