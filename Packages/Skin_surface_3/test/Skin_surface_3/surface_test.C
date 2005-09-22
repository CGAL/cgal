// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Mixed_complex_builder_3.h>
#include <CGAL/Marching_tetrahedra.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Skin_surface_traits_3<>                 Skin_traits;
typedef Skin_traits::Regular_traits                   Regular_traits;
typedef Skin_traits::Regular                          Regular;
typedef Regular_traits::Weighted_point                Reg_weighted_point;
typedef Regular_traits::Bare_point                    Reg_point;
typedef Skin_traits::Simplicial                       Simplicial;
typedef Simplicial::Cell_handle                       Simpl_cell_handle;
typedef Simplicial::Finite_cells_iterator             Simpl_Fin_cells_it;
typedef Simplicial::Finite_vertices_iterator          Simpl_Fin_vertices_it;

typedef CGAL::Mixed_complex_builder_3<Skin_traits>    Mixed_complex_builder;

typedef Skin_traits::Mesh                             Mesh;
typedef Skin_traits::Mesh_K::RT                       Mesh_rt;
typedef CGAL::Marching_tetrahedra_3<Simplicial, Mesh> Marching_tetrahedra;

typedef CGAL::Skin_surface_extractor<
  Simplicial, Mesh, Skin_traits::S2M_converter>       Extractor;
#include <fstream>

int main(int argc, char *argv[]) {
  double shrink = .85;
  Extractor extr;
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
    Simplicial simplicial;
    Mixed_complex_builder(regular, simplicial, shrink);
    for (Simpl_Fin_vertices_it vit = simplicial.finite_vertices_begin();
	 vit != simplicial.finite_vertices_end(); vit++) {
      if (simplicial.is_infinite(vit->cell())) {
	std::cerr << "ERROR: is_infinite (main)" << std::endl;
      }
      Mesh_rt val = extr.value(vit->cell(), vit->point());
      std::cout << vit->cell()->surf->dimension() << " - "
		<< val
		<< std::endl;
      std::list<Simpl_cell_handle> cells;
      simplicial.incident_cells(vit, std::back_inserter(cells));
      for (std::list<Simpl_cell_handle>::iterator cell = cells.begin();
	   cell != cells.end(); cell++) {
	if (!simplicial.is_infinite(*cell)) {
	  std::cout << vit->cell()->surf->dimension() << " "
		    << (*cell)->surf->dimension() << " - "
		    << (extr.value(*cell, vit->point())/val) << std::endl;
	}
      }
    }
  }

  return 0;
}
