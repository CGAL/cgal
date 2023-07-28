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

  std::ifstream ifs("../../test/Polygon_repair_2/data/in/nesting.wkt");
  std::ofstream ofs("/Users/ken/Downloads/triangulation.geojson");

  std::string in;
  std::getline(ifs, in);
  std::istringstream iss(in);
  Multipolygon_with_holes_2 mp;
  CGAL::IO::read_multi_polygon_WKT(iss, mp);
  Polygon_repair_2 pr;
  pr.add_to_triangulation(mp);
  pr.label_triangulation();

  // ofs << std::fixed;
  // ofs << std::setprecision(15);

  // ofs << "{" << std::endl;
  // ofs << "\t\"type\": \"FeatureCollection\"," << std::endl;
  // ofs << "\t\"features\": [" << std::endl;

  // for (Polygon_repair_2::Triangulation::Finite_faces_iterator face = pr.triangulation().finite_faces_begin();
  //      face != pr.triangulation().finite_faces_end(); ++face) {
  //   ofs << "\t\t{" << std::endl;
  //   ofs << "\t\t\t\"type\": \"Feature\"," << std::endl;
  //   ofs << "\t\t\t\"properties\": {" << std::endl;
  //   ofs << "\t\t\t\t\"label\": " << face->label() << std::endl;
  //   ofs << "\t\t\t}," << std::endl;
  //   ofs << "\t\t\t\"geometry\": {" << std::endl;
  //   ofs << "\t\t\t\t\"type\": \"Polygon\"," << std::endl;
  //   ofs << "\t\t\t\t\"coordinates\": [" << std::endl;
  //   ofs << "\t\t\t\t\t[" << std::endl;
  //   ofs << "\t\t\t\t\t\t[" << face->vertex(0)->point().x() << ", " << face->vertex(0)->point().y() << "]," << std::endl;
  //   ofs << "\t\t\t\t\t\t[" << face->vertex(1)->point().x() << ", " << face->vertex(1)->point().y() << "]," << std::endl;
  //   ofs << "\t\t\t\t\t\t[" << face->vertex(2)->point().x() << ", " << face->vertex(2)->point().y() << "]" << std::endl;
  //   ofs << "\t\t\t\t\t]" << std::endl;
  //   ofs << "\t\t\t\t]" << std::endl;
  //   ofs << "\t\t\t}" << std::endl;
  //   ofs << "\t\t}";
  //   Polygon_repair_2::Triangulation::Finite_faces_iterator next_face = face;
  //   ++next_face;
  //   if (next_face != pr.triangulation().finite_faces_end()) ofs << ",";
  //   ofs << std::endl;
  // }

  // ofs << "\t]" << std::endl;
  // ofs << "}" << std::endl;

  pr.compute_nesting();
  pr.reconstruct_multipolygon();
  Multipolygon_with_holes_2 rmp = pr.multipolygon();
  // std::ostringstream oss;
  // CGAL::IO::write_multi_polygon_WKT(oss, rmp);
  // std::string out = oss.str();

  // ofs << std::fixed;
  // ofs << std::setprecision(15);

  
  // ofs << "\t\"coordinates\": [" << std::endl;
  // for (int i = 0; i < rmp.polygons().size(); ++i) {
  //   ofs << "\t\t[" << std::endl;

  //   ofs << "\t\t\t[" << std::endl;
  //   for (int j = 0; j < rmp.polygons()[i].outer_boundary().size(); ++j) {
  //     ofs << "\t\t\t\t[" << rmp.polygons()[i].outer_boundary()[j].x() <<
  //     ", " << rmp.polygons()[i].outer_boundary()[j].y() << "]";
  //     if (j < rmp.polygons()[i].outer_boundary().size()-1) ofs << ",";
  //     ofs << std::endl;
  //   } ofs << "\t\t\t]";
  //   if (rmp.polygons()[i].number_of_holes() > 0) ofs << ",";
  //   ofs << std::endl;

  //   for (int j = 0; j < rmp.polygons()[i].holes().size(); ++j) {
  //     ofs << "\t\t\t[" << std::endl;
  //     for (int k = 0; k < rmp.polygons()[i].holes()[j].size(); ++k) {
  //       ofs << "\t\t\t\t[" << rmp.polygons()[i].holes()[j][k].x() <<
  //       ", " << rmp.polygons()[i].holes()[j][k].y() << "]";
  //       if (k < rmp.polygons()[i].holes()[j].size()-1) ofs << ",";
  //       ofs << std::endl;
  //     }
  //     ofs << "\t\t\t]";
  //     if (j < rmp.polygons()[i].holes().size()-1) ofs << ",";
  //     ofs << std::endl;
  //   }

  //   ofs << "\t\t]";
  //   if (i < rmp.polygons().size()-1) ofs << ",";
  //   ofs << std::endl;
  // } ofs << "\t]" << std::endl;
  

  return 0;
}
