#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>
#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

// default triangulation for Surface_mesher
using Tr =  CGAL::Surface_mesh_default_triangulation_3;

using C2t3 = CGAL::Complex_2_in_triangulation_3<Tr>;
using Sphere_3 = Kernel::Sphere_3;
using Point_3 =  Kernel::Point_3;
using FT =  Kernel::FT;

typedef FT (*Function)(Point_3);
using  Surface_3 = CGAL::Implicit_surface_3<Kernel, Function>;
using Surface_mesh =  CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

FT sphere_function (Point_3 p){
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}

int main() {

  Tr tr;
  C2t3 c2t3(tr);

  Surface_3 surface(sphere_function, Sphere_3(CGAL::ORIGIN, 2.));
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 0.1, 0.1);
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  Surface_mesh sm;
  Surface_mesh deformed;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
  deformed = sm;

  Surface_mesh quad_cage;

  const Point_3 p0(2, -2, -2), p0_new(5, -5, -5);
  const Point_3 p1(2, 2, -2), p1_new(3, 3, -3);
  const Point_3 p2(-2, 2, -2), p2_new(-2, 2, -2);
  const Point_3 p3(-2, -2, -2), p3_new(-3, -3, -3);

  const Point_3 p4(-2, -2, 2), p4_new(-3, -3, 3);
  const Point_3 p5(2, -2, 2), p5_new(4, -4, 4);
  const Point_3 p6(2, 2, 2), p6_new(2, 2, 3);
  const Point_3 p7(-2, 2, 2), p7_new(-3, 3, 3);

  CGAL::make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, quad_cage);
  PMP::triangulate_faces(faces(quad_cage), quad_cage);

  CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Surface_mesh, Kernel> mv(quad_cage);
  auto vertex_to_point_map = get_property_map(CGAL::vertex_point, deformed);

  std::vector<FT> coords;
  std::vector<Point_3> target_cube{p0_new, p1_new, p2_new, p3_new,
                                   p4_new, p5_new, p6_new, p7_new};

  for(auto& v : vertices(deformed)){

    const Point_3 vertex_val = get(vertex_to_point_map, v);
    coords.clear();
    mv(vertex_val, std::back_inserter(coords));

    FT x = FT(0), y = FT(0), z = FT(0);
    for(std::size_t i = 0; i < 8; i++){

      x += target_cube[i].x() * coords[i];
      y += target_cube[i].y() * coords[i];
      z += target_cube[i].z() * coords[i];
    }

    put(vertex_to_point_map, v, Point_3(x, y, z));
  }

  std::ofstream out_original("sphere.off");
  out_original << sm << std::endl;

  std::ofstream out_deformed("deformed_sphere.off");
  out_deformed << deformed << std::endl;

  return EXIT_SUCCESS;
}
