#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#include <list>
#include <string>
#include <fstream>
#include <algorithm>

typedef CGAL::Exact_predicates_exact_constructions_kernel   K;
typedef CGAL::Skin_surface_traits_3<K>                     Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;

typedef Skin_surface_3::Bare_point                          Bare_point;
typedef Skin_surface_3::Weighted_point                      Weighted_point;

typedef CGAL::Exact_predicates_inexact_constructions_kernel IK;
typedef CGAL::Polyhedron_3<IK>                              Polyhedron;

class Test_file
{
public:
  Test_file(double shrink) : s(shrink) {
  }

  void operator()(std::string filename) {
    std::cout << filename << std::endl;

    std::list<Weighted_point> l;
    std::ifstream in(filename.c_str());
    assert(in.is_open());
    Weighted_point wp;
    double x, y, z, w;
    while (in >> x >> y >> z >> w ) {
      l.push_front(Weighted_point(Bare_point(x,y,z),w));
    }

    Skin_surface_3 skin_surface(l.begin(), l.end(), s);

    Polyhedron p;
    CGAL::mesh_skin_surface_3(skin_surface, p);

    assert(p.is_valid() && p.is_closed());

    //std::cout << p << std::endl;
  }

private:
  double s;
};

int main(int, char **)
{
  std::vector<std::string> filenames;
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

  std::for_each(filenames.begin(), filenames.end(), Test_file(.5));
  std::for_each(filenames.begin(), filenames.end(), Test_file(.85));

  return 0;
}
