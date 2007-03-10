#define CGAL_SURFACE_MESHER_VERBOSE 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Skin_surface_3.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel;
typedef CGAL::Regular_triangulation_euclidean_traits_3<Kernel> Traits;
typedef CGAL::Skin_surface_3<Traits>                           Skin_surface_3;
typedef Skin_surface_3::RT                                     RT;
typedef Skin_surface_3::Weighted_point                         Weighted_point;
typedef Weighted_point::Point                                  Bare_point;

typedef CGAL::Surface_mesh_vertex_base_3<Kernel>            Vb;
typedef CGAL::Surface_mesh_cell_base_3<Kernel>              Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>        Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>         Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;
typedef Kernel::Point_3                                     Point_3;
typedef Kernel::FT                                          FT;


int main() {
  // Construct the Skin_surface_3 object:
  std::list<Weighted_point> l;
  RT                        shrinkfactor = 0.5;

  l.push_front(Weighted_point(Bare_point(0,0,0), 1));
  l.push_front(Weighted_point(Bare_point(0,1,0), 2));
  l.push_front(Weighted_point(Bare_point(0,0,2), 1));

  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  // Construct the types needed for the surface mesher:
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.1,  // radius bound
                                                     0.1); // distance bound
  // meshing surface
  make_surface_mesh(c2t3, skin_surface, criteria, CGAL::Non_manifold_tag());

  std::ofstream out("delaunay_mesh.off");
  CGAL::output_surface_facets_to_off(out, c2t3);

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
