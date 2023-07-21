#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>
#include <CGAL/Polygon_repair_2/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes;
typedef CGAL::Multipolygon_with_holes_2<Kernel> Multipolygon;
typedef CGAL::Polygon_repair_2::Polygon_repair_2<Kernel> Polygon_repair;

int main(int argc, char* argv[]) {
  // std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);

  // Square
  // Point ps[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1)};
  // Polygon_2 p(ps, ps+4);
  // Multipolygon mp;
  // mp.add_polygon(p);

  // Bowtie
  Point ps[] = {Point(0,0), Point(1,1), Point(1,0), Point(0,1)};
  Polygon_2 p(ps, ps+4);
  Multipolygon mp;
  mp.add_polygon(p);

  // Overlapping edge
  // Point ps1[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1)};
  // Polygon_2 p1(ps1, ps1+4);
  // Point ps2[] = {Point(1,0), Point(2,0), Point(2,1), Point(1,1)};
  // Polygon_2 p2(ps2, ps2+4);
  // Multipolygon mp;
  // mp.add_polygon(p1);
  // mp.add_polygon(p2);

  // Edge partly overlapping (start)
  // Point ps1[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1)};
  // Polygon_2 p1(ps1, ps1+4);
  // Point ps2[] = {Point(1,0), Point(2,0), Point(2,0.5), Point(1,0.5)};
  // Polygon_2 p2(ps2, ps2+4);
  // Multipolygon mp;
  // mp.add_polygon(p1);
  // mp.add_polygon(p2);

  // Edge partly overlapping (middle)
  // Point ps1[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1)};
  // Polygon_2 p1(ps1, ps1+4);
  // Point ps2[] = {Point(1,0.25), Point(2,0.25), Point(2,0.75), Point(1,0.75)};
  // Polygon_2 p2(ps2, ps2+4);
  // Multipolygon mp;
  // mp.add_polygon(p1);
  // mp.add_polygon(p2);

  // Square with hole
  // Point ps1[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1)};
  // Polygon_with_holes p(Polygon_2(ps1, ps1+4));
  // Point ps2[] = {Point(0.25,0.25), Point(0.75,0.25), Point(0.75,0.75), Point(0.25,0.75)};
  // Polygon_2 h(ps2, ps2+4);
  // p.add_hole(h);
  // Multipolygon mp;
  // mp.add_polygon(p);

  // Square with hole touching boundary
  // Point ps1[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1)};
  // Polygon_with_holes p(Polygon_2(ps1, ps1+4));
  // Point ps2[] = {Point(0.25,0.25), Point(0.75,0.25), Point(0.75,0.75), Point(0,1)};
  // Polygon_2 h(ps2, ps2+4);
  // p.add_hole(h);
  // Multipolygon mp;
  // mp.add_polygon(p);

  // Square with hole touching boundary (self-intersecting loop)
  // Point ps[] = {Point(0,0), Point(1,0), Point(1,1), Point(0,1),
  //               Point(0.25,0.25), Point(0.75,0.25), Point(0.75,0.75), Point(0,1)};
  // Polygon_2 p(ps, ps+8);
  // std::cout << p << std::endl;
  // Multipolygon mp;
  // mp.add_polygon(p);

  // Polygon_repair pr;
  // pr.add_to_triangulation(mp);
  // pr.label_triangulation();

  // for (auto f = pr.triangulation().interior_faces_begin(); f != pr.triangulation().interior_faces_end(); ++f) {
  //   std::cout << f->label() << " ";
  // } std::cout << std::endl;

  // pr.reconstruct_multipolygon();
  // std::cout << pr.multipolygon() << std::endl;

  // CGAL::IO::write_polygon_WKT(std::cout, p);
  // CGAL::IO::write_multi_polygon_WKT(std::cout, mp);

  Multipolygon rmp = CGAL::Polygon_repair_2::repair(p);
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
