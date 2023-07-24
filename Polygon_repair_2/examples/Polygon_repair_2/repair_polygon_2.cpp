#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>
#include <CGAL/Polygon_repair_2/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair_2 = CGAL::Polygon_repair_2::Polygon_repair_2<Kernel>;

int main(int argc, char* argv[]) {
  // std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);

  // Square
  // Point_2 ps[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  // Polygon_2 p(ps, ps+4);
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p);

  // Bowtie
  Point_2 ps[] = {Point_2(0,0), Point_2(1,1), Point_2(1,0), Point_2(0,1)};
  Polygon_2 p(ps, ps+4);
  Multipolygon_with_holes_2 mp;
  mp.add_polygon(p);

  // Overlapping edge
  // Point_2 ps1[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  // Polygon_2 p1(ps1, ps1+4);
  // Point_2 ps2[] = {Point_2(1,0), Point_2(2,0), Point_2(2,1), Point_2(1,1)};
  // Polygon_2 p2(ps2, ps2+4);
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p1);
  // mp.add_polygon(p2);

  // Edge partly overlapping (start)
  // Point_2 ps1[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  // Polygon_2 p1(ps1, ps1+4);
  // Point_2 ps2[] = {Point_2(1,0), Point_2(2,0), Point_2(2,0.5), Point_2(1,0.5)};
  // Polygon_2 p2(ps2, ps2+4);
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p1);
  // mp.add_polygon(p2);

  // Edge partly overlapping (middle)
  // Point_2 ps1[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  // Polygon_2 p1(ps1, ps1+4);
  // Point_2 ps2[] = {Point_2(1,0.25), Point_2(2,0.25), Point_2(2,0.75), Point_2(1,0.75)};
  // Polygon_2 p2(ps2, ps2+4);
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p1);
  // mp.add_polygon(p2);

  // Square with hole
  // Point_2 ps1[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  // Polygon_with_holes_2 p(Polygon_2(ps1, ps1+4));
  // Point_2 ps2[] = {Point_2(0.25,0.25), Point_2(0.75,0.25), Point_2(0.75,0.75), Point_2(0.25,0.75)};
  // Polygon_2 h(ps2, ps2+4);
  // p.add_hole(h);
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p);

  // Square with hole touching boundary
  // Point_2 ps1[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1)};
  // Polygon_with_holes_2 p(Polygon_2(ps1, ps1+4));
  // Point_2 ps2[] = {Point_2(0.25,0.25), Point_2(0.75,0.25), Point_2(0.75,0.75), Point_2(0,1)};
  // Polygon_2 h(ps2, ps2+4);
  // p.add_hole(h);
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p);

  // Square with hole touching boundary (self-intersecting loop)
  // Point_2 ps[] = {Point_2(0,0), Point_2(1,0), Point_2(1,1), Point_2(0,1),
  //               Point_2(0.25,0.25), Point_2(0.75,0.25), Point_2(0.75,0.75), Point_2(0,1)};
  // Polygon_2 p(ps, ps+8);
  // std::cout << p << std::endl;
  // Multipolygon_with_holes_2 mp;
  // mp.add_polygon(p);

  // Polygon_repair_2 pr;
  // pr.add_to_triangulation(mp);
  // pr.label_triangulation();

  // for (auto f = pr.triangulation().interior_faces_begin(); f != pr.triangulation().interior_faces_end(); ++f) {
  //   std::cout << f->label() << " ";
  // } std::cout << std::endl;

  // pr.reconstruct_multipolygon();
  // std::cout << pr.multipolygon() << std::endl;

  // CGAL::IO::write_polygon_WKT(std::cout, p);
  // CGAL::IO::write_multi_polygon_WKT(std::cout, mp);

  Multipolygon_with_holes_2 rmp = CGAL::Polygon_repair_2::repair(p);
  CGAL::IO::write_multi_polygon_WKT(std::cout, rmp);
  CGAL::draw(rmp);

  // std::cout << "Orientations good?" << std::endl;
  // for (auto const& polygon: rmp.polygons()) {
  //   std::cout << (polygon.outer_boundary().orientation() == CGAL::COUNTERCLOCKWISE) << std::endl;
  //   for (auto const &hole: polygon.holes()) {
  //     std::cout << (hole.orientation() == CGAL::CLOCKWISE) << std::endl;
  //   }
  // }
  

  return 0;
}
