#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using EK = CGAL::Exact_predicates_exact_constructions_kernel;
using SM = CGAL::Surface_mesh<K::Point_3>;
using ESM = CGAL::Surface_mesh<EK::Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;

template<class K>
typename K::Aff_transformation_3
rotation(double a, double b, double c)
{
  double ca = cos(a), cb = cos(b), cc = cos(c);
  double sa = sin(a), sb = sin(b), sc = sin(c);

  typename K::Aff_transformation_3 aff(cb * cc, cc* sa* sb - ca * sc, ca* cc* sb + sa * sc,
                                       cb* sc, ca* cc + sa * sb * sc, ca* sb* sc - cc * sa,
                                       -sb, cb* sa, ca* cb);
  return aff;
}

template<class Mesh>
void test_kernel_on_mesh(const Mesh &input, std::size_t expected_nb_vertices, std::size_t expected_nb_edges, std::size_t expected_nb_faces, double expected_volume = 0){
  // CGAL_assertion(PMP::is_kernel_empty(input) == (expected_nb_vertices == 0));
  // CGAL_assertion((PMP::kernel_point(input, CGAL::parameters::allow_non_manifold_non_watertight_input(true)) != std::nullopt) == (expected_nb_vertices < 3 && expected_nb_faces != 1));
  // CGAL_assertion((PMP::kernel_point(input, CGAL::parameters::allow_non_manifold_non_watertight_input(true).require_strictly_inside(false)) != std::nullopt) == (expected_nb_vertices != 0));

  Mesh kernel = kernel = PMP::kernel(input, CGAL::parameters::allow_non_manifold_non_watertight_input(true).use_bounding_box_filtering(false).shuffle_planes(false).starting_cube(PMP::bbox(input)));
  CGAL_assertion(vertices(kernel).size() == expected_nb_vertices);
  CGAL_assertion(edges(kernel).size()    == expected_nb_edges);
  CGAL_assertion(faces(kernel).size()    == expected_nb_faces);
  if(expected_volume != 0){
    PMP::triangulate_faces(kernel);
    CGAL_assertion(PMP::volume(kernel) > expected_volume * 0.99 && PMP::volume(kernel) < expected_volume * 1.01);
  }
  clear(kernel);

  kernel = kernel = PMP::kernel(input, CGAL::parameters::allow_non_manifold_non_watertight_input(true));
  CGAL_assertion(vertices(kernel).size() == expected_nb_vertices);
  CGAL_assertion(edges(kernel).size()    == expected_nb_edges);
  CGAL_assertion(faces(kernel).size()    == expected_nb_faces);
  if(expected_volume != 0){
    PMP::triangulate_faces(kernel);
    CGAL_assertion(PMP::volume(kernel) > expected_volume * 0.99 && PMP::volume(kernel) < expected_volume * 1.01);
  }
}

template<class Mesh, class K>
void tests(){
  using P = typename K::Point_3;
  Mesh m;

  // A simple cube
  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  test_kernel_on_mesh(m, 8, 12, 8, 1);
  clear(m);

  // rotated cubes
  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  PMP::transform(rotation<K>(0, 0, 1), m);
  test_kernel_on_mesh(m, 8, 12, 8, 1);
  clear(m);

  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  PMP::transform(rotation<K>(0, 1, 1), m);
  test_kernel_on_mesh(m, 8, 12, 8, 1);
  clear(m);

  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  PMP::transform(rotation<K>(1, 1, 1), m);
  test_kernel_on_mesh(m, 8, 12, 8, 1);
  clear(m);

  // Degenerate to a segment
  // make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  // kernel = PMP::internal::kernel(m, faces(m), CGAL::parameters::allow_non_manifold_non_watertight_input(true).use_bounding_box_filtering(false));
  // CGAL_assertion(kernel.vertices().size()==8);

  // Degenerate to a face
  // make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  // make_hexahedron(P(0,0,0).bbox()+P(-1,1,1).bbox(), m);
  // kernel = PMP::internal::kernel(m, faces(m), CGAL::parameters::allow_non_manifold_non_watertight_input(true));
  // CGAL_assertion(kernel.vertices().size()==4);

  // Degenerate to a segment
  // m.clear();
  // kernel.clear();
  // make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  // make_hexahedron(P(0,0,-0.5).bbox()+P(-1,-1,0.5).bbox(), m);
  // kernel = PMP::internal::kernel(m, faces(m), CGAL::parameters::allow_non_manifold_non_watertight_input(true).use_bounding_box_filtering(false));
  // CGAL_assertion(kernel.vertices().size()==2);
  // P a = kernel.point(*kernel.vertices().begin());
  // P b = kernel.point(*(++kernel.vertices().begin()));
  // CGAL_assertion((a==P(0,0,0) && b==P(0,0,0.5)) || (b==P(0,0,0) && a==P(0,0,0.5)));

  // Degenerate to a point
  // m.clear();
  // kernel.clear();
  // make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  // make_hexahedron(P(0,0,0).bbox()+P(-1,-1,-1).bbox(), m);
  // kernel = PMP::internal::kernel(m, faces(m), CGAL::parameters::allow_non_manifold_non_watertight_input(true).use_bounding_box_filtering(false));
  // CGAL_assertion(kernel.vertices().size()==1);
  // a = kernel.point(*kernel.vertices().begin());
  // CGAL_assertion((a==P(0,0,0)));

}

int main(int argc, char** argv)
{
  tests<SM, K>();
  tests<ESM, EK>();
  tests<CGAL::Polyhedron_3<K>,K>();
  tests<CGAL::Polyhedron_3<EK>, EK>();
}
