#include <iostream>
#include <fstream>

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

  std::ifstream ifs("data/in/nesting-spike.wkt");

  Multipolygon_with_holes_2 mp;
  CGAL::IO::read_multi_polygon_WKT(ifs, mp);

  Polygon_repair pr;
  pr.add_to_triangulation_even_odd(mp);
  pr.label_triangulation_even_odd();

  std::cout << "{" << std::endl;
  std::cout << "\t\"type\": \"FeatureCollection\"," << std::endl;
  std::cout << "\t\"features\": [" << std::endl;

  for (Polygon_repair::Triangulation::Finite_faces_iterator face = pr.triangulation().finite_faces_begin();
       face != pr.triangulation().finite_faces_end(); ++face) {
    std::cout << "\t\t{" << std::endl;
    std::cout << "\t\t\t\"type\": \"Feature\"," << std::endl;
    std::cout << "\t\t\t\"properties\": {" << std::endl;
    std::cout << "\t\t\t\t\"label\": " << face->label() << std::endl;
    std::cout << "\t\t\t}," << std::endl;
    std::cout << "\t\t\t\"geometry\": {" << std::endl;
    std::cout << "\t\t\t\t\"type\": \"Polygon\"," << std::endl;
    std::cout << "\t\t\t\t\"coordinates\": [" << std::endl;
    std::cout << "\t\t\t\t\t[" << std::endl;
    std::cout << "\t\t\t\t\t\t[" << face->vertex(0)->point().x() << ", " << face->vertex(0)->point().y() << "]," << std::endl;
    std::cout << "\t\t\t\t\t\t[" << face->vertex(1)->point().x() << ", " << face->vertex(1)->point().y() << "]," << std::endl;
    std::cout << "\t\t\t\t\t\t[" << face->vertex(2)->point().x() << ", " << face->vertex(2)->point().y() << "]" << std::endl;
    std::cout << "\t\t\t\t\t]" << std::endl;
    std::cout << "\t\t\t\t]" << std::endl;
    std::cout << "\t\t\t}" << std::endl;
    std::cout << "\t\t}";
    Polygon_repair::Triangulation::Finite_faces_iterator next_face = face;
    ++next_face;
    if (next_face != pr.triangulation().finite_faces_end()) std::cout << ",";
    std::cout << std::endl;
  }

  std::cout << "\t]" << std::endl;
  std::cout << "}" << std::endl;

  pr.reconstruct_multipolygon();
  Multipolygon_with_holes_2 repaired = pr.multipolygon();

  return 0;
}
