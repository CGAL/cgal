// geometric traits class
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

// triangulation
#include <CGAL/Regular_triangulation_3.h>

// vertex
#include <CGAL/PSC_mesh_vertex_base_3.h>

// cell 
#include <CGAL/Volume_mesh_cell_base_3.h>

// c2t3
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>

// indexed 3D image
// PSC
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Volume_bounded_by_polyhedron_mesh_traits_3.h>

// meshing criteria
#include <CGAL/PSC_mesh_default_edges_criteria_3.h>
#include <CGAL/PSC_mesh_default_facets_criteria_3.h>
#include <CGAL/PSC_mesh_default_cells_criteria_3.h>

// meshing function
#include <CGAL/make_mesh_3_from_piecewise_smooth_surfaces.h>

// input/output
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_medit.h>
#include <fstream>


struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Regular_triangulation_filtered_traits_3<K> Regular_traits;
typedef CGAL::Robust_circumcenter_traits_3<Regular_traits> Traits;

typedef CGAL::PSC_mesh_vertex_base_3<Traits> Vb;
typedef CGAL::Volume_mesh_cell_base_3<Traits> Cb;

typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds> Tr;

typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;

typedef CGAL::Polyhedron_3<Traits> Polyhedron;
typedef CGAL::Volume_bounded_by_polyhedron_mesh_traits_3<Polyhedron> Mesh_traits;

typedef CGAL::PSC_mesh_default_edges_criteria_3<Tr> Edges_criteria;
typedef CGAL::PSC_mesh_default_facets_criteria_3<Tr> Facets_criteria;
typedef CGAL::PSC_mesh_default_tetrahedra_criteria_3<Tr> Tets_criteria;

int main(int, char **) {
  Tr tr;           // 3D Delaunay triangulation
  C2t3 c2t3(tr);   // 2D complex in 3D-Delaunay triangulation

  Polyhedron polyhedron;
  std::ifstream input("input.off");

  
  Edges_criteria edges_criteria(5,  // uniform radius upper bound,
                                1); // distance bound

  Facets_criteria facets_criteria(30, // angle upper bound,
                                  5,  // uniform radius upper bound,
                                  1); // distance bound

  Tets_criteria tets_criteria(4,  // radius-edge ratio upper bound
                              5); // uniform radius bound

  CGAL::make_mesh_3_from_piecewise_smooth_surfaces(c2t3,
                                                   mesh_traits,
                                                   edges_criteria,
                                                   facets_criteria,
                                                   tets_criteria,
                                                   0.5); 
  // 0.5 is radius-radius ratio upper bound for the sliver exudation
  // process
                       
  std::ofstream file_medit("out.mesh");
  CGAL::output_to_medit_file(file_medit, c2t3);
}
