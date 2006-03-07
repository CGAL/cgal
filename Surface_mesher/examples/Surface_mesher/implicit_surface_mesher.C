#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

struct Kernel : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Surface_mesh_vertex_base_3<Kernel> Vb;
typedef CGAL::Surface_mesh_cell_base_3<Kernel> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds> Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Kernel::Sphere_3 Sphere_3;

typedef double (*Function)(double, double, double);
typedef CGAL::Implicit_surface_3<Kernel, Function> Surface_3;

double sphere_function (double x, double y, double z) {
  const double x2=x*x, y2=y*y, z2=z*z;
  return x2+y2+z2-1;
}

int main(int, char **) {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  Surface_3 surface(sphere_function, 
                    Sphere_3(CGAL::ORIGIN, 2.),
                    1e-06);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                     0.05,
                                                     0.01);
  // meshing surface
  make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
