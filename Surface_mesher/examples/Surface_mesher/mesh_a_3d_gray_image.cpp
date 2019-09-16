#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

int main() {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // the 'function' is a 3D gray level image
  Gray_level_image image("data/skull_2.9.inr", 2.9f);

  // Carefully chosen bounding sphere: the center must be inside the
  // surface defined by 'image' and the radius must be high enough so that
  // the sphere actually bounds the whole image.
  GT::Point_3 bounding_sphere_center(122., 102., 117.);
  GT::FT bounding_sphere_squared_radius = 200.*200.*2.;
  GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                   bounding_sphere_squared_radius);

  // definition of the surface, with 10^-5 as relative precision
  Surface_3 surface(image, bounding_sphere, 1e-5);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                     5.,
                                                     5.);

  // meshing surface, with the "manifold without boundary" algorithm
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  std::ofstream out("out.off");
  CGAL::output_surface_facets_to_off (out, c2t3);
  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
