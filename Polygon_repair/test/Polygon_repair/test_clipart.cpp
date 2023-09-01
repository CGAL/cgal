#include <iostream>
#include <sstream>
#include <pugixml.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/Polygon_repair.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main(int argc, char* argv[]) {

  std::string folder = "/Users/ken/Downloads/repaired";

  for (const auto& file: std::filesystem::directory_iterator(folder)) {
    if (file.path().filename().extension() != ".svg") continue;
    if (file.path().filename().stem() != "182377") continue;
    std::cout << "Testing " << file.path().filename() << "... ";

    // Read test file and create multipolygon with outer boundaries
    pugi::xml_document doc;
    if (doc.load_file(file.path().string().c_str())) {
      Multipolygon_with_holes_2 mp;
      std::list<Polygon_2> holes;
      for (auto child: doc.child("svg").children()) {
        if (strcmp(child.name(), "polygon") == 0) {
          Polygon_2 p;
          std::string points(child.attribute("points").value());
          std::string color(child.attribute("fill").value());

          std::istringstream polygonss(points);
          std::string point;
          while (polygonss >> point) {
            std::istringstream pointss(point);
            std::string coordinate;
            double x, y;
            getline(pointss, coordinate, ',');
            x = std::stod(coordinate);
            getline(pointss, coordinate);
            y = std::stod(coordinate);
//            std::cout << "(" << x << ", " << y << ")" << std::endl;
            p.push_back(Point_2(x, y));
          }

          if (color == "black") {
            mp.add_polygon(p);
          } else {
            holes.push_back(p);
          }
        }
      }

      // Put holes in correct polygon
      for (auto const& hole: holes) {
        std::set<std::size_t> matches;
        for (auto const& vertex: hole.vertices()) {
          for (std::size_t pn = 0; pn < mp.number_of_polygons(); ++pn) {
            if (mp.polygons()[pn].outer_boundary().bounded_side(vertex) == CGAL::ON_BOUNDED_SIDE) {
              matches.insert(pn);
            }
          }
        } if (matches.size() == 1) {
          std::cout << "Found match" << std::endl;
          mp.polygons()[*matches.begin()].add_hole(hole);
        } else {
          std::cout << "Error: couldn't find polygon for hole" << std::endl;
        }
      }

      // Check validity
      std::cout << CGAL::Polygon_repair::is_valid(mp) << std::endl;
    }
  }

  return 0;
}
