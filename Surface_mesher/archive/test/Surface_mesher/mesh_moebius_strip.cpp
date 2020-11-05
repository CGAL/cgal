#define CGAL_SURFACE_MESHER_VERBOSE

// traits class
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Delaunay
#include <CGAL/Delaunay_triangulation_3.h>

// vertex and cell bases
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>

#include <cmath>

#include <cassert>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// vertex and cell types
typedef CGAL::Surface_mesh_vertex_base_3<K> Vb;
typedef CGAL::Surface_mesh_cell_base_3<K> Cb;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

FT moebius_function (Point_3 p) {
  const double x = CGAL::to_double(p.x());
  const double y = CGAL::to_double(p.y());
  const double z = CGAL::to_double(p.z());

  const double r = CGAL::sqrt(x*x + y*y);
  const double halftheta = std::atan2(y,x) / 2.;
  return std::log(r)*std::sin(halftheta) - z * cos(halftheta);
//   const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
//   return x2+y2+z2-1;
}

int main() {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  Surface_3 surface(moebius_function,             // pointer to function
                    Sphere_3(Point_3(0.0001, -0.0003, 0.), 2.)); // bounding sphere

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     0.05,  // radius bound
                                                     0.05); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  std::ofstream out("out.off");

#ifndef NDEBUG
  const bool result =
#endif
    CGAL::output_surface_facets_to_off(out, c2t3,
                                       CGAL::Surface_mesher::IO_VERBOSE |
                                       CGAL::Surface_mesher::IO_ORIENT_SURFACE);

  assert(result == false);

  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}
