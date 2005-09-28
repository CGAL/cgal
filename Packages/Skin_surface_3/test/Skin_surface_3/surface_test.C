// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulated_mixed_complex_3.h>
#include <CGAL/triangulate_mixed_complex_3.h>
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/marching_tetrahedra_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Skin_surface_traits_3<>                  Skin_traits;
typedef Skin_traits::Regular_traits                    Regular_traits;

typedef CGAL::Regular_triangulation_3<Regular_traits> Regular;
typedef Regular_traits::Weighted_point                Reg_weighted_point;
typedef Regular_traits::Bare_point                    Reg_point;

typedef CGAL::Triangulated_mixed_complex_3<Skin_traits> Tr2;
typedef Tr2::Cell_handle                       Tr2_cell_handle;
typedef Tr2::Finite_cells_iterator             Tr2_Fin_cells_it;
typedef Tr2::Finite_vertices_iterator          Tr2_Fin_vertices_it;

typedef Skin_traits::Polyhedron_kernel         Polyhedron_kernel;
typedef CGAL::Polyhedron_3<Polyhedron_kernel>  Polyhedron;
typedef Polyhedron_kernel::RT                  Polyhedron_rt;

typedef CGAL::Marching_tetrahedra_traits_skin_surface_3<
  Tr2, Polyhedron, Skin_traits::T2P_converter> Marching_tetrahedra_traits;

#include <fstream>

int main(int argc, char *argv[]) {
  double shrink = .85;
  Skin_traits skin_traits(shrink);
  Marching_tetrahedra_traits marching_traits;
  while (argc>1) {
    argc--; argv++;
    std::ifstream is(argv[0]);
  
    Regular regular;
    Reg_weighted_point wp;
    double max_weight=1;
    CGAL::Bbox_3 box;
    while (is >> wp) {
      max_weight = std::max(max_weight, wp.weight());
      box = box + wp.point().bbox();
      regular.insert(wp);
    }

    // add a bounding octahedron:
    Reg_point mid((box.xmin() + box.xmax())/2,
      (box.ymin() + box.ymax())/2,
      (box.zmin() + box.zmax())/2);
    double size = 1.5*((box.xmax() - box.xmin() +
			 box.ymax() - box.ymin() +
			 box.zmax() - box.zmin())/2 + max_weight);
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x()+size,mid.y(),mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x()-size,mid.y(),mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y()+size,mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y()-size,mid.z()),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()+size),-1));
    regular.insert(
      Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()-size),-1));

    // Triangulate mixed complex:
    Tr2 triangulated_mixed_complex;
    triangulate_mixed_complex_3(regular, triangulated_mixed_complex,skin_traits);
    for (Tr2_Fin_vertices_it vit = triangulated_mixed_complex.finite_vertices_begin();
	 vit != triangulated_mixed_complex.finite_vertices_end(); vit++) {
      if (triangulated_mixed_complex.is_infinite(vit->cell())) {
	std::cerr << "ERROR: is_infinite (main)" << std::endl;
      }
      Polyhedron_rt val = marching_traits.value(vit->cell(), vit->point());
      std::cout << vit->cell()->surf->dimension() << " - "
		<< val
		<< std::endl;
      std::list<Tr2_cell_handle> cells;
      triangulated_mixed_complex.incident_cells(vit, std::back_inserter(cells));
      for (std::list<Tr2_cell_handle>::iterator cell = cells.begin();
	   cell != cells.end(); cell++) {
	if (!triangulated_mixed_complex.is_infinite(*cell)) {
	  std::cout << vit->cell()->surf->dimension() << " "
		    << (*cell)->surf->dimension() << " - "
		    << (marching_traits.value(*cell, vit->point())/val) << std::endl;
	}
      }
    }
  }

  return 0;
}
