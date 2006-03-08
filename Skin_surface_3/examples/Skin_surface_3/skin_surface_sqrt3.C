// examples/Skin_surface_3/skin_surface_sqrt3.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/skin_surface_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/skin_surface_sqrt3_3.h>
#include <CGAL/Skin_surface_refinement_traits_with_face_info_3.h>

#include <CGAL/Marching_tetrahedra_observer_skin_surface_3.h>

#include <list>
#include <fstream>

typedef CGAL::Skin_surface_traits_3<>                    Skin_surface_traits;
typedef Skin_surface_traits::Regular_traits              Regular_traits;
typedef CGAL::Regular_triangulation_3<Regular_traits>    Regular;
typedef Regular_traits::Bare_point                       Reg_point;
typedef Regular_traits::Weighted_point                   Reg_weighted_point;
typedef CGAL::Triangulated_mixed_complex_3<Skin_surface_traits>
                                                     Triangulated_mixed_complex;
typedef Skin_surface_traits::Polyhedron_traits       Polyhedron_kernel;
typedef CGAL::Skin_surface_polyhedral_items_3<
          Triangulated_mixed_complex>                Polyhedral_items;

typedef CGAL::Polyhedron_3<Polyhedron_kernel, Polyhedral_items> Polyhedron;

typedef CGAL::Marching_tetrahedra_traits_skin_surface_3<
  Triangulated_mixed_complex,
  Polyhedron,
  Skin_surface_traits::T2P_converter>  Marching_tetrahedra_traits;
typedef CGAL::Marching_tetrahedra_observer_skin_surface_3<
  Triangulated_mixed_complex, Polyhedron>     Marching_tetrahedra_observer;

typedef CGAL::Skin_surface_refinement_traits_with_face_info_3<
          Triangulated_mixed_complex, 
          Polyhedron, 
          Skin_surface_traits::T2P_converter,
          Skin_surface_traits::P2T_converter>    Skin_surface_refinement_traits;

int main(int argc, char *argv[]) {
  std::list<Reg_weighted_point> l;
  Skin_surface_traits           skin_surface_traits(.5);
  Regular                       regular;
  Triangulated_mixed_complex    triangulated_mixed_complex;
  Polyhedron                    polyhedron;
  
  l.push_front(Reg_weighted_point(Reg_point(0,0,0), 1.1));
  l.push_front(Reg_weighted_point(Reg_point(0,1,0), 2));
  l.push_front(Reg_weighted_point(Reg_point(0,0,2), 1));

  // Construct regular triangulation and bounding box
  regular.insert(l.begin(), l.end());
  CGAL::skin_surface_construct_bounding_box_3(regular, skin_surface_traits);

  // Construct the triangulated mixed complex
  CGAL::triangulate_mixed_complex_3(
    regular, triangulated_mixed_complex, skin_surface_traits);

  // Extract the coarse mesh using marching_tetrahedra
  Marching_tetrahedra_traits marching_traits;
  Marching_tetrahedra_observer marching_observer;
  CGAL::marching_tetrahedra_3(
    triangulated_mixed_complex, polyhedron, marching_traits, marching_observer);

  // Subdivide mesh:
  Skin_surface_refinement_traits refinement_traits(triangulated_mixed_complex);
  CGAL::skin_surface_sqrt3(polyhedron, refinement_traits, 3);
  
  return 0;
}
