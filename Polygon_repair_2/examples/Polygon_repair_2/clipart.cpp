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

  std::string folder_in = "/Volumes/Toshiba/out_fix";
  std::string folder_out = "/Volumes/Toshiba/repaired";
  double desired_width = 500.0;
  int current = 0, how_many = 10;

  for (const auto& file: std::filesystem::directory_iterator(folder_in)) {
    if (file.path().filename().extension() != ".obj") continue;
    std::cout << "Reading " << file.path().filename() << "..." << std::endl;

    Polygon_repair_2 pr;
    std::vector<Point_2> vertices;
    std::vector<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>> edges;

    std::ifstream ifs(file.path());
    std::string line;

    while (std::getline(ifs, line)) {
      std::istringstream iss(line);

      char c;
      iss >> c;

      if (c == 'v') {
        double x, y;
        iss >> x >> y;
        vertices.emplace_back(x, y);
      } else if (c == 'l') {
        unsigned int a, b;
        iss >> a >> b;
        edges.push_back(std::make_pair(vertices[a-1], vertices[b-1]));
      }
    }

    ifs.close();

    std::unordered_set<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>,
                       boost::hash<std::pair<typename Kernel::Point_2, typename Kernel::Point_2>>> edges_to_insert;

    for (auto const& edge: edges) {
      if (edge.first == edge.second) continue;
      std::pair<typename Kernel::Point_2, typename Kernel::Point_2> pair = (edge.first < edge.second)?
      std::make_pair(edge.first, edge.second) : std::make_pair(edge.second, edge.first);
      auto inserted = edges_to_insert.insert(pair);
      if (!inserted.second) edges_to_insert.erase(inserted.first);
    }

    for (auto const& edge: edges_to_insert) {
      pr.triangulation().odd_even_insert_constraint(edge.first, edge.second);
    }

    // std::cout << pr.triangulation().number_of_faces() << " faces in the triangulation" << std::endl;

    if (pr.triangulation().number_of_faces() > 0) {
      pr.label_triangulation();
      pr.reconstruct_multipolygon();
    } Multipolygon_with_holes_2 mp = pr.multipolygon();

    // std::cout << mp << std::endl;
    // std::cout << mp.number_of_polygons() << " polygons" << std::endl;

    if (mp.number_of_polygons() > 0) {
      CGAL::Bbox_2 bbox = mp.polygons().front().bbox();
      for (auto const& polygon: mp.polygons()) {
        bbox += polygon.outer_boundary().bbox();
      } // std::cout << bbox.xmin() << " " << bbox.xmax() << " " << bbox.ymin() << " " << bbox.ymax() << std::endl;
      Kernel::Vector_2 translate(-bbox.xmin(), -bbox.ymin());
      double scale = desired_width/(bbox.xmax()-bbox.xmin());


      std::ofstream ofs(folder_out + "/" + std::string(file.path().stem()) + ".svg");
      ofs << "<svg viewBox=\"0 0 " << desired_width << " " << scale*(bbox.ymax()-bbox.ymin()) << "\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;

      for (auto const& polygon: mp.polygons()) {
        // std::cout << polygon << std::endl;
        ofs << "\t<polygon points=\"";
        for (auto const& vertex: polygon.outer_boundary()) {
          ofs << scale*(vertex.x()+translate.x()) << "," << scale*(vertex.y()+translate.y()) << " ";
        } ofs << "\" fill=\"black\"/>" << std::endl;
      }

      for (auto const& polygon: mp.polygons()) {
        for (auto const& hole: polygon.holes()) {
          // std::cout << hole << std::endl;
          ofs << "\t<polygon points=\"";
          for (auto const& vertex: hole) {
            ofs << scale*(vertex.x()+translate.x()) << "," << scale*(vertex.y()+translate.y()) << " ";
          } ofs << "\" fill=\"white\"/>" << std::endl;
        }
      }

      ofs << "</svg>";
      ofs.close();
    } ++current;
    
//    if (current >= how_many) break;
  }

  return 0;
}
