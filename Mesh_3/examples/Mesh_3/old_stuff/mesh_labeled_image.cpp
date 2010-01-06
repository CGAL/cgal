// geometric traits class
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>

// triangulation
#include <CGAL/Regular_triangulation_3.h>

// vertex
#include <CGAL/Surface_mesh_vertex_base_3.h>

// cell 
#include <CGAL/Multilabel_mesh_cell_base_3.h>

// c2t3
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>

// indexed 3D image
#include <CGAL/Labeled_image_3.h>
#include <CGAL/Labeled_image_volume_mesh_traits_3.h>

// meshing criteria
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Volume_mesh_default_criteria_3.h>

// meshing function
#include <CGAL/make_mesh_3_for_multivolumes.h>

// output
#include <CGAL/IO/File_medit.h>
#include <fstream>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Regular_triangulation_filtered_traits_3<K> Traits;

typedef CGAL::Surface_mesh_vertex_base_3<Traits> Vb;
typedef CGAL::Multilabel_mesh_cell_base_3<Traits, unsigned char> Cb;

typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds> Tr;

typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;

typedef CGAL::Labeled_image_3<unsigned char> Image_3;
typedef CGAL::Labeled_image_volume_mesh_traits_3<Image_3,Traits> Mesh_traits;

typedef CGAL::Surface_mesh_default_criteria_3<Tr> Facets_criteria;
typedef CGAL::Volume_mesh_default_criteria_3<Tr> Tets_criteria;

int main(int, char **) {
  Tr tr;           // 3D Delaunay triangulation
  C2t3 c2t3(tr);   // 2D complex in 3D-Delaunay triangulation

  // Loads image
  Image_3 image("segmented_head.inr.gz");
  Mesh_traits mesh_traits(image);

  Facets_criteria facets_criteria(30, // angle upper bound,
                                  5,  // uniform radius upper bound,
                                  1); // distance bound

  Tets_criteria tets_criteria(4,  // radius-edge ratio upper bound
                              5); // uniform radius bound

  CGAL::make_mesh_3_for_multivolumes(c2t3,
                                     mesh_traits,
                                     facets_criteria,
                                     tets_criteria,
                                     0.5); // radius-radius ratio upper
                                           // bound for the sliver
                                           // exudation process
   
  std::ofstream medit_file("out.mesh");
  CGAL::output_to_medit_file(medit_file, c2t3);
}
