#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

#include <CGAL/enum.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_3                                                  Point;

typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<K>          Gt;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Gt>                Tr;

int main(int, char**)
{
  Gt traits(Point(1, 2, 3), 5);

  std::vector<Point> points;
  points.emplace_back( 2, 1, 1);
  points.emplace_back(-2, 1, 1); // not on the sphere
  points.emplace_back( 0, 1, 1);
  points.emplace_back( 1, 2, 1);
  points.emplace_back( 0, 1, 1); // duplicate
  points.emplace_back( 1, 0, 1);
  points.emplace_back( 1, 1, 2);

  // test constructors
  Tr tr0;
  assert(tr0.geom_traits().center() == Point(0, 0, 0));
  assert(tr0.geom_traits().radius() == 1);

  Tr tr(traits);
  assert(tr.geom_traits().center() == Point(1, 2, 3));
  assert(tr.geom_traits().radius() == 5);

  Tr trbis(tr);
  assert(trbis.geom_traits().center() == tr.geom_traits().center());
  assert(trbis.geom_traits().radius() == tr.geom_traits().radius());
  assert(trbis.number_of_vertices() == tr.number_of_vertices());
  assert(trbis.number_of_edges() == tr.number_of_edges());
  assert(trbis.number_of_solid_faces() == tr.number_of_solid_faces());
  assert(trbis.number_of_ghost_faces() == tr.number_of_ghost_faces());

  Tr tr2(tr);
  assert(tr2.geom_traits().center() == tr.geom_traits().center());
  assert(tr2.geom_traits().radius() == tr.geom_traits().radius());

  Tr tr3 = tr2;
  assert(tr3.geom_traits().center() == tr2.geom_traits().center());
  assert(tr3.geom_traits().radius() == tr2.geom_traits().radius());

  Tr tr4(points.begin(), points.end(), Point(1, 1, 1), 1);
  assert(tr4.number_of_vertices() == 5);

  tr4.set_radius(0.123);
  assert(tr4.number_of_vertices() == 0);

  Tr tr5(points.begin(), points.end());
  assert(tr5.geom_traits().center() == Point(0, 0, 0));
  assert(tr5.geom_traits().radius() == 1);
  assert(tr5.number_of_vertices() == 0);

  Gt traits2(Point(1, 1, 1), 1);
  Tr tr6(points.begin(), points.end(), traits2);
  assert(tr6.number_of_vertices() == 5);

  tr6.set_center(Point(1,2,3));
  assert(tr6.number_of_vertices() == 0);

  // center / radius setting
  tr.insert(Point(6, 2, 3));
  assert(tr.number_of_vertices() == 1);

  tr.set_center_and_radius(Point(1, 1, 1), 1);
  assert(tr3.number_of_vertices() == 0);
  assert(tr.geom_traits().center() == Point(1, 1, 1));
  assert(tr.geom_traits().radius() == 1);

  tr.insert(points.begin(), points.end());
  CGAL::IO::write_OFF("test_dtos.off", tr);
  assert(tr.is_valid());

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // test ranges
  Tr::Vertex_handles vs = tr.vertex_handles();
  assert(vs.size() == tr.number_of_vertices() && vs.size() == 5);

  Tr::All_edges es = tr.all_edges();
  Tr::Solid_edges ses = tr.solid_edges();
  assert(es.size() == 9 && ses.size() == 8);

  Tr::All_face_handles afs = tr.all_face_handles();
  assert(afs.size() == tr.number_of_faces() && afs.size() == 6);

  Tr::Solid_face_handles sfs = tr.solid_faces();
  assert(sfs.size() == 4);
  assert(tr.number_of_solid_faces() + tr.number_of_ghost_faces() == tr.number_of_faces());

  Tr::Face_handle fh = *tr.solid_faces().begin();
  for(auto f: tr.solid_faces()) {
    assert(f == fh);
    break;
  }

  Tr::Points pts = tr.points();
  assert(pts.size() == tr.number_of_vertices() && pts.size() == 5);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // test iterators

  //   All vertices iterator
  Tr::Vertices_iterator vit = tr.vertices_begin(), vend = tr.vertices_end();
  assert(std::distance(vit, vend) == 5);

  const Tr::Point& p0 = tr.point(vit);
  ++vit;
  const Tr::Point& p1 = tr.point(vit++);
  assert(p0 != p1);
  std::advance(vit, -2);
  assert(p0 == vit->point());

  //   All edges iterator
  Tr::All_edges_iterator eit = tr.all_edges_begin(), eend = tr.all_edges_end();
  assert(std::distance(eit, eend) == 9);

  Tr::Edge e0 = *eit;
  ++eit;
  Tr::Edge e1 = *eit++;
  assert(e0 != e1);
  std::advance(eit, -2);
  assert(e0 == *eit);

  //   Solid edges iterator
  Tr::Solid_edges_iterator seit = tr.solid_edges_begin(), seend = tr.solid_edges_end();
  assert(std::distance(seit, seend) == 8);

  Tr::Edge se0 = *seit;
  ++seit;
  Tr::Edge se1 = *seit++;
  assert(se0 != se1);
  std::advance(seit, -2);
  assert(se0 == *seit);

  //   Contour edges iterator
  Tr::Contour_edges_iterator ceit=tr.contour_edges_begin(), ceend=tr.contour_edges_end();
  assert(std::distance(ceit, ceend) == 4);

  std::set<Tr::Edge> unique_edges;
  while(ceit != ceend)
  {
    const Tr::Edge& e = *ceit;
    unique_edges.insert(e);
    assert(tr.is_ghost(e.first) != tr.is_ghost(tr.mirror_edge(e).first)); // xor
    ++ceit;
  }

  assert(unique_edges.size() == 4);
  std::advance(ceit, -4);
  assert(ceit == tr.contour_edges_begin());

  //   All faces iterator
  Tr::All_faces_iterator fit = tr.all_faces_begin(), fend = tr.all_faces_end();
  assert(std::distance(fit, fend) == 6);

  Tr::Face_handle f0 = fit;
  ++fit;
  Tr::Face_handle f1 = fit++;
  assert(f0 != f1);
  std::advance(fit, -2);
  assert(f0 == fit);

  //   Solid faces iterator
  Tr::Solid_faces_iterator sfit = tr.solid_faces_begin(), sfend = tr.solid_faces_end();
  assert(std::distance(sfit, sfend) == 4);

  Tr::Face_handle sf0 = sfit;
  ++sfit;
  Tr::Face_handle sf1 = sfit++;
  assert(sf0 != sf1);
  std::advance(sfit, -2);
  Tr::Face_handle sf2 = sfit;
  assert(sf0 == sf2);

  //   Face Circulator
  vit = tr.vertices_begin();
  std::advance(vit, 3); // top of the pyramid
  Tr::Face_circulator fc = tr.incident_faces(vit), end = fc;
  assert(std::distance(++fc, end) + 1 == 4);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
