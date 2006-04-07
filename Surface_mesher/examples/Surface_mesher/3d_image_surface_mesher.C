// file examples/Surface_mesher/3d_image_surface_mesher.C
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

typedef CGAL::Gray_level_image_3<Kernel::FT, Kernel::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<Kernel, Gray_level_image> Surface_3;

int main(int, char **) {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // the 'function' is a 3D gray level image
  Gray_level_image image("ImageIO/data/skull_2.9.inr.gz", 2.9);

  // Carefully choosen bounding sphere: the center must be inside the
  // surface defined by 'image' and the radius must be high enough so that
  // the sphere actually bounds the whole image.
  Kernel::Point_3 bounding_sphere_center(122., 102., 117.);
  Kernel::FT bounding_sphere_squared_radius = 200.*200.*2.;
  Kernel::Sphere_3 bounding_sphere(bounding_sphere_center,
                                   bounding_sphere_squared_radius);

  // definition of the surface, with 10^-2 as relative precision
  Surface_3 surface(image, bounding_sphere, 1e-2);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                     5.,
                                                     5.);

  // meshing surface, with the "manifold without boundary" algorithm
  make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
