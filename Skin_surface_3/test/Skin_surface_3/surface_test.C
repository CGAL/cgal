// test/Skin_surface_3/surface_test.C
// #include <CGAL/Skin_surface_traits_3.h>
// #include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// #include <CGAL/Regular_triangulation_3.h>
// #include <CGAL/Triangulated_mixed_complex_3.h>
// #include <CGAL/triangulate_mixed_complex_3.h>
// #include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
// #include <CGAL/marching_tetrahedra_3.h>
// #include <CGAL/Polyhedron_3.h>

// typedef CGAL::Skin_surface_traits_3<>                  Skin_traits;
// typedef Skin_traits::Regular_traits                    Regular_traits;

// typedef CGAL::Regular_triangulation_3<Regular_traits> Regular;
// typedef Regular_traits::Weighted_point                Reg_weighted_point;
// typedef Regular_traits::Bare_point                    Reg_point;

// typedef CGAL::Triangulated_mixed_complex_3<Skin_traits> Tr2;
// typedef Tr2::Cell_handle                       Tr2_cell_handle;
// typedef Tr2::Finite_cells_iterator             Tr2_Fin_cells_it;
// typedef Tr2::Finite_vertices_iterator          Tr2_Fin_vertices_it;

// typedef Skin_traits::Polyhedron_traits         Polyhedron_kernel;
// typedef CGAL::Polyhedron_3<Polyhedron_kernel>  Polyhedron;
// typedef Polyhedron_kernel::RT                  Polyhedron_rt;

// typedef CGAL::Marching_tetrahedra_traits_skin_surface_3<
//   Tr2, Polyhedron, Skin_traits::T2P_converter> Marching_tetrahedra_traits;

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
  double shrink = .85;
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
	  CGAL_assertion(val == val2);
	}
      }
    }
  }

  return 0;
}
