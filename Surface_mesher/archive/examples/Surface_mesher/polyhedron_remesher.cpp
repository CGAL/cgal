#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_piecewise_smooth_surface_mesh.h>

#include <CGAL/Polyhedral_surface_3.h>
#include <fstream>

#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Piecewise_smooth_surface_mesh_default_edges_criteria_3.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef CGAL::Polyhedral_surface_3<GT> Polyhedral_surface;


int main() {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  std::ifstream file_input("data/triceratops.off");
  Polyhedral_surface surface(file_input);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr>
    facets_criteria(30.,  // angular bound
                    0.5,  // radius bound
                    0.5); // distance bound
  CGAL::Surface_mesh_default_edges_criteria_3<Tr>
    edges_criteria(0.5,  // radius bound
                   0.5); // distance bound

  // meshing surface
  CGAL::make_piecewise_smooth_surface_mesh(c2t3, surface,
                                           facets_criteria, edges_criteria,
                                           CGAL::Manifold_tag());

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
