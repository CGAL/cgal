// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/make_skin_surface_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K1;
typedef CGAL::Exact_predicates_exact_constructions_kernel   K2;
typedef K1::FT                                              FT;
typedef K1::Point_3                                         Point_3;
typedef CGAL::Weighted_point<Point_3,FT>                    Weighted_point;
typedef CGAL::Polyhedron_3<K1>                              Polyhedron;
typedef CGAL::Skin_surface_3<K2>                            Skin_surface_3;

int main(int argc, char *argv[]) {
  std::list<Weighted_point> l;
  FT                        shrinkfactor = 0.5;
  Skin_surface_3            skin_surface;
  Polyhedron                p;
  
  l.push_front(Weighted_point(Point_3(0,0,0), 1));
  l.push_front(Weighted_point(Point_3(0,1,0), 2));
  l.push_front(Weighted_point(Point_3(0,0,2), 1));

  CGAL::make_skin_surface_3(l.begin(), l.end(), shrinkfactor, skin_surface);
  CGAL::mesh_skin_surface_3(skin_surface, p);

  return 0;
}
