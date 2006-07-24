// test/Skin_surface_3/subdivision_test.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <list>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::RT                                  RT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

typedef CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> Poly_items_skin;
typedef CGAL::Polyhedron_3<K,Poly_items_skin>                 Polyhedron_skin;

#include <fstream>

Skin_surface_3 create_skin_surface(char * filename, double shrink) {
  std::list<Weighted_point> l;
  std::ifstream in(filename);
  CGAL_assertion(in.is_open());
  Weighted_point wp;
  while (in >> wp) l.push_front(wp);

  return Skin_surface_3(l.begin(), l.end(), shrink);
}

template < class Skin_surface_3, class Polyhedron>
bool construct_and_subdivide_mesh(Skin_surface_3 &skin_surface,
				  Polyhedron &polyhedron) 
{
  CGAL::mesh_skin_surface_3(skin_surface, polyhedron);
  CGAL_assertion(polyhedron.is_valid() && polyhedron.is_closed());

  CGAL::subdivide_skin_surface_mesh_3(polyhedron, skin_surface);
  return (polyhedron.is_valid() && polyhedron.is_closed());
}

int main(int argc, char *argv[]) {
  bool result; 

  std::vector<char *> filenames;
  filenames.push_back("data/caffeine.cin");
  filenames.push_back("data/ball.cin");
  filenames.push_back("data/degenerate.cin");
  filenames.push_back("data/test1.cin");
  filenames.push_back("data/test2.cin");
  filenames.push_back("data/test3.cin");
  filenames.push_back("data/test4.cin");
  filenames.push_back("data/test5.cin");
  filenames.push_back("data/test6.cin");
  filenames.push_back("data/test7.cin");
  filenames.push_back("data/test8.cin");
  filenames.push_back("data/test9.cin");
  filenames.push_back("data/test10.cin");
  filenames.push_back("data/test11.cin");

  Polyhedron p;
  Polyhedron_skin p_skin;
  for (std::vector<char *>::size_type i=0; i<filenames.size(); i++) {
    Skin_surface_3 skin_surface = create_skin_surface(filenames[i], .5);
    Polyhedron p;

    result = construct_and_subdivide_mesh(skin_surface, p);
    CGAL_assertion(result);
    p.clear();

    result = construct_and_subdivide_mesh(skin_surface, p_skin);
    CGAL_assertion(result);
    p_skin.clear();
  }

  return 0;
}
