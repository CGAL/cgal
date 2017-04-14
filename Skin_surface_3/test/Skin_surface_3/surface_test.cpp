#define  CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/triangulate_mixed_complex_3.h>

#include <list>
#include <string>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Inexact_K;
typedef CGAL::Exact_predicates_exact_constructions_kernel   Exact_K;
typedef CGAL::Skin_surface_traits_3<Inexact_K>              Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface;
typedef Skin_surface::Simplex                               Simplex;
typedef Skin_surface::FT                                    FT;

typedef Skin_surface::Bare_point                            Bare_point;
typedef Skin_surface::Weighted_point                        Weighted_point;

typedef CGAL::Polyhedron_3<Traits>                          Polyhedron;

CGAL::Cartesian_converter<Exact_K,Inexact_K> e2i_converter;
CGAL::Cartesian_converter<Inexact_K,Exact_K> i2e_converter;

// Triangulated_mixed_complex:
typedef Skin_surface::TMC                                 TMC;

typedef TMC::Vertex_iterator TMC_Vertex_iterator;
typedef TMC::Cell_iterator   TMC_Cell_iterator;
typedef TMC::Vertex_handle   TMC_Vertex_handle;
typedef TMC::Cell_handle     TMC_Cell_handle;
typedef TMC::Finite_vertices_iterator TMC_Finite_vertices_iterator;
typedef TMC::Finite_cells_iterator TMC_Finite_cells_iterator;

#include <fstream>

class Test_file {
public:
  Test_file(double shrink) : s(shrink) {
  }

  void operator()(std::string & filename)
  {
    std::cout << s << " " << filename << std::endl;

    std::ifstream in(filename.c_str());

    std::list<Weighted_point> l;
    double x,y,z,w;
    while (in >> x >> y >> z >> w)
      l.push_front(Weighted_point(Bare_point(x,y,z),w));

    Skin_surface skin_surface(l.begin(), l.end(), s);

    TMC &tmc = skin_surface.triangulated_mixed_complex();
    //     CGAL::Triangulated_mixed_complex_observer_3<TMC, Skin_surface>
    //       observer(skin_surface.shrink_factor());
    //     triangulate_mixed_complex_3(skin_surface.get_regular_triangulation(),
    //                                 skin_surface.shrink_factor(),
    //                                 tmc,
    //                                 observer,
    //                                 false);

    for (TMC_Finite_vertices_iterator vit = tmc.finite_vertices_begin();
         vit != tmc.finite_vertices_end(); vit++) {

      if (tmc.is_infinite(vit->cell())) {
        std::cerr << "ERROR: is_infinite (main)" << std::endl;
      }
      Exact_K::FT val = vit->cell()->info().second->value(vit->point());
      std::list<TMC_Cell_handle> cells;
      tmc.incident_cells(vit, std::back_inserter(cells));
      for (std::list<TMC_Cell_handle>::iterator cell =
           cells.begin();
           cell != cells.end(); cell++) {
        if (!tmc.is_infinite(*cell)) {
          Exact_K::FT val2 = (*cell)->info().second->value(vit->point());
          // NGHK: Make exact:
          //assert(val == val2);
          assert(std::abs(CGAL::to_double(val-val2)) < 1e-8);
        }
      }
    }
  }

private:
  double s;
};

int main(int , char **)
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
  filenames.push_back("data/degenerate.cin");

  std::for_each(filenames.begin(), filenames.end(), Test_file(.5));
  std::for_each(filenames.begin(), filenames.end(), Test_file(.85));

  return 0;
}
