#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>

#define CGAL_SURFACE_MESHER_TEST 1
#include <CGAL/Surface_mesh_triangulation_generator_3.h> // undocumented

#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Implicit_surface_3.h>

#include <iostream>

double sphere(double x, double y, double z)
{
  return x*x+y*y+z*z-1.;
}

struct Kernel : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Surface_mesh_triangulation_generator_3<Kernel>::Type Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2T3;

typedef Kernel::Sphere_3 Sphere_3;

using CGAL::make_surface_mesh;

int main(int argc, char **argv) {
  Tr tr;

  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 0.05, 0.01);

  std::cout << "Initial number of points: " << tr.number_of_vertices() 
            << std::endl;

  // 2D-complex in 3D-Delaunay triangulation
  C2T3 c2t3 (tr);

  // Surface meshing
  make_surface_mesh(c2t3,
		    CGAL::make_implicit_surface_3(Kernel(),
						  sphere,
						  Sphere_3(CGAL::ORIGIN, 2.),
						  1e-06),
		    criteria,
		    CGAL::Non_manifold_tag());
  
  std::cout << "Final number of points: " << tr.number_of_vertices() 
            << std::endl;
}
