#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

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

  // std::ifstream ifs("/Users/ken/Downloads/180927.wkt");
  // std::ofstream ofs("/Users/ken/Downloads/1.geojson");

  std::ifstream ifs("/Users/ken/Downloads/2018418.wkt");
  std::ofstream ofs("/Users/ken/Downloads/2.geojson");

  std::string in;
  std::getline(ifs, in);
  std::istringstream iss(in);
  Polygon_with_holes_2 p;
  CGAL::IO::read_polygon_WKT(iss, p);
  Polygon_repair pr;
  pr.add_to_triangulation_even_odd(p);
  pr.label_triangulation_even_odd();
  pr.reconstruct_multipolygon();
  Multipolygon_with_holes_2 rmp = pr.multipolygon();
  std::ostringstream oss;
  CGAL::IO::write_multi_polygon_WKT(oss, rmp);
  std::string out = oss.str();

  ofs << std::fixed;
  ofs << std::setprecision(15);

  ofs << "{" << std::endl;
  ofs << "\t\"type\": \"MultiPolygon\"," << std::endl;
  ofs << "\t\"coordinates\": [" << std::endl;
  for (int i = 0; i < rmp.polygons_with_holes().size(); ++i) {
    ofs << "\t\t[" << std::endl;

    ofs << "\t\t\t[" << std::endl;
    for (int j = 0; j < rmp.polygons_with_holes()[i].outer_boundary().size(); ++j) {
      ofs << "\t\t\t\t[" << rmp.polygons_with_holes()[i].outer_boundary()[j].x() <<
      ", " << rmp.polygons_with_holes()[i].outer_boundary()[j].y() << "]";
      if (j < rmp.polygons_with_holes()[i].outer_boundary().size()-1) ofs << ",";
      ofs << std::endl;
    } ofs << "\t\t\t]";
    if (rmp.polygons_with_holes()[i].number_of_holes() > 0) ofs << ",";
    ofs << std::endl;

    for (int j = 0; j < rmp.polygons_with_holes()[i].holes().size(); ++j) {
      ofs << "\t\t\t[" << std::endl;
      for (int k = 0; k < rmp.polygons_with_holes()[i].holes()[j].size(); ++k) {
        ofs << "\t\t\t\t[" << rmp.polygons_with_holes()[i].holes()[j][k].x() <<
        ", " << rmp.polygons_with_holes()[i].holes()[j][k].y() << "]";
        if (k < rmp.polygons_with_holes()[i].holes()[j].size()-1) ofs << ",";
        ofs << std::endl;
      }
      ofs << "\t\t\t]";
      if (j < rmp.polygons_with_holes()[i].holes().size()-1) ofs << ",";
      ofs << std::endl;
    }

    ofs << "\t\t]";
    if (i < rmp.polygons_with_holes().size()-1) ofs << ",";
    ofs << std::endl;
  } ofs << "\t]" << std::endl;
  ofs << "}" << std::endl;

  return 0;
}
