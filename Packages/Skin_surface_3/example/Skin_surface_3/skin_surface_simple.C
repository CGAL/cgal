// examples/Skin_surface_3/skin_surface_simple.C

#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

typedef CGAL::Skin_surface_traits_3<>                 Skin_traits;
typedef CGAL::Regular_triangulation_euclidean_traits_3<Skin_traits::Regular_K>
                                                      Regular_traits;
typedef CGAL::Regular_triangulation_3<Regular_traits> Regular;
typedef Regular_traits::Weighted_point                Reg_weighted_point;

int main(int argc, char *argv[]) {
  Regular regular;

  std::ifstream is("./data/caffeine.cin");
  Reg_weighted_point wp;
  while (is >> wp) regular.insert(wp);
}
