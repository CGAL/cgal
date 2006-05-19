#define  CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#include <list>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Skin_surface_3::Bare_point                          Bare_point;
typedef CGAL::Polyhedron_3<Traits>                          Polyhedron;
typedef Skin_surface_3::Quadratic_surface                   Quadratic_surface;


typedef Skin_surface_3::Triangulated_mixed_complex  Triangulated_mixed_complex;
typedef Triangulated_mixed_complex::Vertex_handle   Tmc_Vertex_handle;
typedef Triangulated_mixed_complex::Finite_vertices_iterator 
                                                    Tmc_Finite_vertices_iterator;

#include <fstream>

int main(int argc, char *argv[]) {
  double shrink = .5;
  std::list<Weighted_point> l;

  while (argc>1) {
    argc--; argv++;
    std::ifstream in(argv[0]);

    std::list<Weighted_point> l;
    Weighted_point wp;
    while (in >> wp) l.push_front(wp);
    
    Skin_surface_3 skin_surface(l.begin(), l.end(), shrink);

    const Skin_surface_3::Triangulated_mixed_complex &tmc = 
      skin_surface.triangulated_mixed_complex();

    for (Tmc_Finite_vertices_iterator vit = tmc.finite_vertices_begin();
	 vit != tmc.finite_vertices_end(); vit++) {
      if (tmc.is_infinite(vit->cell())) {
	std::cerr << "ERROR: is_infinite (main)" << std::endl;
      }
      Quadratic_surface::RT val = vit->cell()->surf->value(vit->point());
      std::list<Triangulated_mixed_complex::Cell_handle> cells;
      tmc.incident_cells(vit, std::back_inserter(cells));
      for (std::list<Triangulated_mixed_complex::Cell_handle>::iterator cell = 
	     cells.begin();
	   cell != cells.end(); cell++) {
	if (!tmc.is_infinite(*cell)) {
	  Quadratic_surface::RT val2 = (*cell)->surf->value(vit->point());
 	  CGAL_assertion(std::abs(val - val2) < 10e-6);
	}
      }
    }
  }

  return 0;
}
