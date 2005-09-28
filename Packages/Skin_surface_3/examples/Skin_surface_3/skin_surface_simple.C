// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>

#include <list>

typedef CGAL::Skin_surface_traits_3<>                    Skin_surface_traits;
typedef Skin_surface_traits::Regular_traits              Regular_traits;
typedef Regular_traits::Bare_point                       Reg_point;
typedef Regular_traits::Weighted_point                   Reg_weighted_point;
typedef CGAL::Polyhedron_3<Skin_surface_traits::Polyhedron_kernel> Polyhedron;

int main(int argc, char *argv[]) {
  std::list<Reg_weighted_point> l;
  
  l.push_front(Reg_weighted_point(Reg_point(0,0,0), 1));
  l.push_front(Reg_weighted_point(Reg_point(0,1,0), 2));
  l.push_front(Reg_weighted_point(Reg_point(0,0,2), 1));

  Polyhedron p;
  Skin_surface_traits skin_surface_traits(.5);
  CGAL::skin_surface_3(l.begin(), l.end(), p, skin_surface_traits);

  return 0;
}
