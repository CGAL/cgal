#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_surface_mesh.h>

#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

struct Kernel : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Surface_mesh_vertex_base_3<Kernel> Vb;
typedef CGAL::Surface_mesh_cell_base_3<Kernel> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;
typedef Kernel::Sphere_3 Sphere_3;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Gray_level_image_3<Kernel::FT> Function;
typedef CGAL::Implicit_surface_3<Kernel, Function> Surface_3;

int main(int, char **) {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  Function function("ImageIO/data/skull_2.9.inr.gz");
  Surface_3 surface(function, 
                    Sphere_3(Point_3(250., 250., 250.), 500.),
                    1e-03);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                     0.1,
                                                     0.1);
  // meshing surface
  make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
