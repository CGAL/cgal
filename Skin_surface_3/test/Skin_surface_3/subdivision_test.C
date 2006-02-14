// test/Skin_surface_3/subdivision_test.C
#include <CGAL/skin_surface_3.h>
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulated_mixed_complex_3.h>
#include <CGAL/triangulate_mixed_complex_3.h>
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>
#include <CGAL/Marching_tetrahedra_observer_skin_surface_3.h>
#include <CGAL/marching_tetrahedra_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Skin_surface_refinement_traits_3.h>
#include <CGAL/Skin_surface_refinement_traits_with_face_info_3.h>

#include <CGAL/Skin_surface_polyhedral_items_3.h>

#include <CGAL/skin_surface_sqrt3_3.h>

typedef CGAL::Skin_surface_traits_3<>                  Skin_traits;
typedef Skin_traits::Regular_traits                    Regular_traits;

typedef CGAL::Regular_triangulation_3<Regular_traits> Regular;
typedef Regular_traits::Weighted_point                Reg_weighted_point;
typedef Regular_traits::Bare_point                    Reg_point;

typedef CGAL::Triangulated_mixed_complex_3<Skin_traits> Tr2;
typedef Tr2::Cell_handle                       Tr2_cell_handle;
typedef Tr2::Finite_cells_iterator             Tr2_Fin_cells_it;
typedef Tr2::Finite_vertices_iterator          Tr2_Fin_vertices_it;

typedef Skin_traits::Polyhedron_traits         Polyhedron_kernel;
typedef CGAL::Polyhedron_3<Polyhedron_kernel>  Polyhedron;

typedef CGAL::Skin_surface_polyhedral_items_3<Tr2> Polyhedral_items;
typedef CGAL::Polyhedron_3<Polyhedron_kernel, Polyhedral_items>
                                                     Polyhedron_plus;

typedef Polyhedron_kernel::RT                  Polyhedron_rt;

typedef CGAL::Marching_tetrahedra_traits_skin_surface_3<
  Tr2, Polyhedron, Skin_traits::T2P_converter> Marching_tetrahedra_traits;

typedef CGAL::Marching_tetrahedra_observer_default_3<
  Tr2, Polyhedron>     Marching_observer_default;
typedef CGAL::Marching_tetrahedra_observer_skin_surface_3<
  Tr2, Polyhedron_plus>     Marching_observer_skin_surface;

typedef CGAL::Skin_surface_refinement_traits_3<
          Tr2, 
          Polyhedron, 
          Skin_traits::T2P_converter,
          Skin_traits::P2T_converter>    Skin_refinement_traits;

typedef CGAL::Skin_surface_refinement_traits_with_face_info_3<
          Tr2, 
          Polyhedron_plus, 
          Skin_traits::T2P_converter,
          Skin_traits::P2T_converter>    Skin_refinement_traits_plus;


#include <fstream>


template < class Triangulated_mixed_complex,
	   class Polyhedron,
	   class Marching_traits,
	   class Marching_observer,
	   class Skin_surface_refinement_traits>
void construct_and_subdivide_mesh(
  Triangulated_mixed_complex &triangulated_mixed_complex,
  Polyhedron &polyhedron,
  Marching_traits &marching_traits,
  Marching_observer &marching_observer,
  Skin_surface_refinement_traits &refinement_traits) {

  CGAL::marching_tetrahedra_3(
    triangulated_mixed_complex, polyhedron, marching_traits, marching_observer);
  
  CGAL::skin_surface_sqrt3(polyhedron, refinement_traits, 4);
}

int main(int argc, char *argv[]) {
  double shrink = .85;
  Skin_traits skin_traits(shrink);
  Marching_tetrahedra_traits marching_traits;
  while (argc>1) {
    argc--; argv++;
    std::ifstream is(argv[0]);
  
    Regular regular;
    Reg_weighted_point wp;

    while (is >> wp) regular.insert(wp);
    CGAL::skin_surface_construct_bounding_box_3(regular, skin_traits);

    // Triangulate mixed complex:
    Tr2 triangulated_mixed_complex;
    CGAL::triangulate_mixed_complex_3(
      regular, triangulated_mixed_complex, skin_traits);

    
    Polyhedron                  polyhedron;
    Marching_tetrahedra_traits  marching_traits;
    Marching_observer_default   marching_observer;
    Skin_refinement_traits      refinement_traits(triangulated_mixed_complex);

    // Without additional face information in the polyhedron
    construct_and_subdivide_mesh(
      triangulated_mixed_complex, polyhedron, 
      marching_traits, marching_observer,
      refinement_traits);
    
    Polyhedron_plus polyhedron_plus;
    Marching_observer_skin_surface marching_observer_plus;
    Skin_refinement_traits_plus
      refinement_traits_plus(triangulated_mixed_complex);

    // With additional face information in the polyhedron
    construct_and_subdivide_mesh(
      triangulated_mixed_complex, polyhedron_plus, 
      marching_traits, marching_observer_plus,
      refinement_traits_plus);
    
  }

  return 0;
}
