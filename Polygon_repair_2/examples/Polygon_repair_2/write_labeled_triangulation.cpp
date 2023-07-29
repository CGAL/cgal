#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

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

  // std::ifstream ifs("/Users/ken/Downloads/180927.wkt");
  // std::ofstream ofs("/Users/ken/Downloads/2.geojson");

  // std::ifstream ifs("/Users/ken/Downloads/2018418.wkt");
  // std::ofstream ofs("/Users/ken/Downloads/1.geojson");

  std::ifstream ifs("../../test/Polygon_repair_2/data/in/nesting-spike.wkt");
  std::ofstream ofs("/Users/ken/Downloads/triangulation.geojson");

  std::string in;
  std::getline(ifs, in);
  std::istringstream iss(in);
  Polygon_repair_2 pr;

  if (in.find("POLYGON") == 0) {
    Polygon_with_holes_2 p;
    if (in != "POLYGON()") { // maybe should be checked in WKT reader
      CGAL::IO::read_polygon_WKT(iss, p);
    } pr.add_to_triangulation(p);
  } else if (in.find("MULTIPOLYGON") == 0) {
    Multipolygon_with_holes_2 mp;
    CGAL::IO::read_multi_polygon_WKT(iss, mp);
    pr.add_to_triangulation(mp);
  } pr.label_triangulation();

  // ofs << std::fixed;
  // ofs << std::setprecision(15);

  ofs << "{" << std::endl;
  ofs << "\t\"type\": \"FeatureCollection\"," << std::endl;
  ofs << "\t\"features\": [" << std::endl;

  for (Polygon_repair_2::Triangulation::Finite_faces_iterator face = pr.triangulation().finite_faces_begin();
       face != pr.triangulation().finite_faces_end(); ++face) {
    ofs << "\t\t{" << std::endl;
    ofs << "\t\t\t\"type\": \"Feature\"," << std::endl;
    ofs << "\t\t\t\"properties\": {" << std::endl;
    ofs << "\t\t\t\t\"label\": " << face->label() << std::endl;
    ofs << "\t\t\t}," << std::endl;
    ofs << "\t\t\t\"geometry\": {" << std::endl;
    ofs << "\t\t\t\t\"type\": \"Polygon\"," << std::endl;
    ofs << "\t\t\t\t\"coordinates\": [" << std::endl;
    ofs << "\t\t\t\t\t[" << std::endl;
    ofs << "\t\t\t\t\t\t[" << face->vertex(0)->point().x() << ", " << face->vertex(0)->point().y() << "]," << std::endl;
    ofs << "\t\t\t\t\t\t[" << face->vertex(1)->point().x() << ", " << face->vertex(1)->point().y() << "]," << std::endl;
    ofs << "\t\t\t\t\t\t[" << face->vertex(2)->point().x() << ", " << face->vertex(2)->point().y() << "]" << std::endl;
    ofs << "\t\t\t\t\t]" << std::endl;
    ofs << "\t\t\t\t]" << std::endl;
    ofs << "\t\t\t}" << std::endl;
    ofs << "\t\t}";
    Polygon_repair_2::Triangulation::Finite_faces_iterator next_face = face;
    ++next_face;
    if (next_face != pr.triangulation().finite_faces_end()) ofs << ",";
    ofs << std::endl;
  }

  ofs << "\t]" << std::endl;
  ofs << "}" << std::endl;

  pr.reconstruct_multipolygon();
  Multipolygon_with_holes_2 rmp = pr.multipolygon();
  CGAL::IO::write_multi_polygon_WKT(std::cout, rmp);

  return 0;
}
