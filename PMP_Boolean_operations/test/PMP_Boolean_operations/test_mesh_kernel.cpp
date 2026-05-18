#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/boost/graph/Euler_operations.h>
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


std::size_t i=0;
template<class Mesh, class Plane_3>
void test_clip_convex_on_mesh(const Mesh &m, const Plane_3 &pl, std::size_t expected_nb_vertices, std::size_t expected_nb_edges, std::size_t expected_nb_faces){
  using K =typename CGAL::Kernel_traits<typename Plane_3::value_type>::Kernel;
  using Pl = typename K::Plane_3;
  using Pl_3P = typename PMP::internal::Three_point_cut_plane_traits<K>::Plane_3;

  auto m_copy = m;
  PMP::internal::clip_convex(m_copy, Pl_3P(pl[0], pl[1], pl[2]), CGAL::parameters::geom_traits(PMP::internal::Three_point_cut_plane_traits<K>()));

  assert(vertices(m_copy).size() == expected_nb_vertices);
  assert(edges(m_copy).size()    == expected_nb_edges);
  assert(faces(m_copy).size()    == expected_nb_faces);

  m_copy = m;
  PMP::internal::clip_convex(m_copy, Pl(pl[0], pl[1], pl[2]));

  assert(vertices(m_copy).size() == expected_nb_vertices);
  assert(edges(m_copy).size()    == expected_nb_edges);
  assert(faces(m_copy).size()    == expected_nb_faces);
}

template<class Mesh>
void test_kernel_on_mesh(const Mesh &input, std::size_t expected_nb_vertices, std::size_t expected_nb_edges, std::size_t expected_nb_faces, double expected_volume = 0){
  assert(PMP::has_empty_kernel(input, CGAL::parameters::allow_open_input(true)) == (expected_nb_vertices == 0));
  assert((PMP::kernel_point(input, CGAL::parameters::allow_open_input(true).strictly_inside(false)) != std::nullopt) == (expected_nb_vertices != 0));

  Mesh kernel;
  PMP::kernel(input, kernel, CGAL::parameters::allow_open_input(true));
#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "nb of vertices: " << vertices(kernel).size() << " ( " << expected_nb_vertices << " expected)" << std::endl
            << "nb of edges: " << edges(kernel).size() << " ( " << expected_nb_edges << " expected)" << std::endl
            << "nb of faces: " << faces(kernel).size() << " ( " << expected_nb_faces << " expected)" << std::endl;
#endif
  assert(vertices(kernel).size() == expected_nb_vertices);
  assert(edges(kernel).size()    == expected_nb_edges);
  assert(faces(kernel).size()    == expected_nb_faces);
  if(expected_volume != 0){
    PMP::triangulate_faces(kernel);
#ifdef TEST_MESH_KERNEL_VERBOSE
    std::cout << "volume: " << PMP::volume(kernel) << " ( " << expected_volume << " expected)" << std::endl;
#endif
    assert(PMP::volume(kernel) > expected_volume * 0.99 && PMP::volume(kernel) < expected_volume * 1.01);
  }

  clear(kernel);
  using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;
  std::map<face_descriptor, face_descriptor> storage;
  boost::associative_property_map< std::map<face_descriptor, face_descriptor> > f2f_map(storage);
  auto bb = CGAL::Polygon_mesh_processing::bbox(input);
  bb = CGAL::Bbox_3(bb.xmin()-3,bb.ymin()-3,bb.zmin()-3,bb.xmax()+3,bb.ymax()+3,bb.zmax()+3);
  make_hexahedron(bb, kernel);
  PMP::kernel(faces(input), input, kernel, CGAL::parameters::allow_open_input(true).use_bounding_box_filtering(false).shuffle_planes(false), CGAL::parameters::face_to_face_map(f2f_map));
#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "nb of vertices: " << vertices(kernel).size() << " ( " << expected_nb_vertices << " expected)" << std::endl
            << "nb of edges: " << edges(kernel).size() << " ( " << expected_nb_edges << " expected)" << std::endl
            << "nb of faces: " << faces(kernel).size() << " ( " << expected_nb_faces << " expected)" << std::endl;
#endif
  assert(vertices(kernel).size() == expected_nb_vertices);
  assert(edges(kernel).size()    == expected_nb_edges);
  assert(faces(kernel).size()    == expected_nb_faces);
  for(face_descriptor f: faces(kernel))
    assert(get(f2f_map, f) != boost::graph_traits<Mesh>::null_face());
  if(expected_volume != 0){
    PMP::triangulate_faces(kernel);
#ifdef TEST_MESH_KERNEL_VERBOSE
    std::cout << "volume: " << PMP::volume(kernel) << " ( " << expected_volume << " expected)" << std::endl;
#endif
    assert(PMP::volume(kernel) > expected_volume * 0.99 && PMP::volume(kernel) < expected_volume * 1.01);
  }
}

template<class Mesh, class K>
void tests(){
  using P = typename K::Point_3;
  using Pl = std::array<P, 3>;
  auto opposite = [](const Pl &a){
    return Pl({a[0], a[2], a[1]});
  };

  Mesh m;

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "\n Test " << typeid(Mesh).name() << "\n" << std::endl;
  std::cout << "Test cube" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  test_kernel_on_mesh(m, 8, 12, 6, 1);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test rotated cube" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  PMP::transform(rotation<K>(0, 0, 1), m);
  test_kernel_on_mesh(m, 8, 12, 6, 1);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test rotated tetrahedron 1" << std::endl;
#endif
  make_tetrahedron(P(0,0,0),
                   P(0,0,1),
                   P(1,0,0),
                   P(0,1,0),
                   m);
  PMP::transform(rotation<K>(0, 1, 1), m);
  test_kernel_on_mesh(m, 4, 6, 4, 1./6.);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test rotated tetrahedron 2" << std::endl;
#endif
  make_tetrahedron(P(0,0,0),
                   P(0,0,1),
                   P(1,0,0),
                   P(0,1,0),
                   m);
  PMP::transform(rotation<K>(1, 1, 1), m);
  test_kernel_on_mesh(m, 4, 6, 4, 1./6.);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test star" << std::endl;
#endif
  make_tetrahedron(P(0,0,2),
                   P(0,2,0),
                   P(2,0,0),
                   P(0,0,0),
                   m);
  auto f1 = *(faces(m).begin());
  auto f2 = *(++faces(m).begin());
  auto f3 = *(++(++faces(m).begin()));
  auto f4 = *(++(++(++(faces(m).begin()))));
  auto v1 = CGAL::Euler::add_center_vertex(halfedge(f1, m), m);
  auto v2 = CGAL::Euler::add_center_vertex(halfedge(f2, m), m);
  auto v3 = CGAL::Euler::add_center_vertex(halfedge(f3, m), m);
  auto v4 = CGAL::Euler::add_center_vertex(halfedge(f4, m), m);
  put(CGAL::get_property_map(CGAL::vertex_point, m), target(v1, m), P( 4, 4, 4));
  put(CGAL::get_property_map(CGAL::vertex_point, m), target(v2, m), P( 1,-6, 1));
  put(CGAL::get_property_map(CGAL::vertex_point, m), target(v3, m), P( 1, 1,-6));
  put(CGAL::get_property_map(CGAL::vertex_point, m), target(v4, m), P(-6, 1, 1));
  test_kernel_on_mesh(m, 8, 18, 12, 2.2963);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test rotated star" << std::endl;
#endif
  PMP::transform(rotation<K>(1, 1, 1), m);
  test_kernel_on_mesh(m, 8, 18, 12, 2.2963);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test degenerated to a face 1" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  make_hexahedron(P(0,0,0).bbox()+P(-1,1,1).bbox(), m);
  test_kernel_on_mesh(m, 4, 4, 1, 0);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test degenerated to a face 2" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(2,2,2).bbox(), m);
  make_hexahedron(P(0,1,1).bbox()+P(-2,3,3).bbox(), m);
  test_kernel_on_mesh(m, 4, 4, 1, 0);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test degenerated to a segment 1" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  make_hexahedron(P(0,0,0).bbox()+P(-1,-1,1).bbox(), m);
  test_kernel_on_mesh(m, 2, 0, 0, 0);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test degenerated to a segment 2" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(2,2,2).bbox(), m);
  make_hexahedron(P(0,0,1).bbox()+P(-2,-2,3).bbox(), m);
  test_kernel_on_mesh(m, 2, 0, 0, 0);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test degenerated to a vertex" << std::endl;
#endif
  make_hexahedron(P(0,0,0).bbox()+P(2,2,2).bbox(), m);
  make_hexahedron(P(0,0,0).bbox()+P(-2,-2,-2).bbox(), m);
  test_kernel_on_mesh(m, 1, 0, 0, 0);
  clear(m);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test empty" << std::endl;
#endif
  make_hexahedron(P(1,1,1).bbox()+P(3,3,3).bbox(), m);
  make_hexahedron(P(0,0,0).bbox()+P(-2,-2,-2).bbox(), m);
  test_kernel_on_mesh(m, 0, 0, 0, 0);
  clear(m);

  make_hexahedron(P(0,0,0).bbox()+P(1,1,1).bbox(), m);
  Pl pl;
#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip general" << std::endl;
#endif
  pl = {P(1,1,0.5),P(1,0,0.5),P(0,1,0.5)};
  test_clip_convex_on_mesh(m, pl, 8, 12, 6);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip one vertex on" << std::endl;
#endif
  pl = {P(1,1,1),P(1,0,0.8),P(0,1,0.8)};
  test_clip_convex_on_mesh(m, pl, 7, 11, 6);
  test_clip_convex_on_mesh(m, opposite(pl), 8, 12, 6);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip two vertices on" << std::endl;
#endif
  pl = {P(1,1,1),P(0,0,1),P(0,1,0.8)};
  test_clip_convex_on_mesh(m, pl, 4, 6, 4);
  test_clip_convex_on_mesh(m, opposite(pl), 8, 13, 7);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip three vertices on" << std::endl;
#endif
  pl = {P(1,1,1),P(0,0,1),P(0,1,0)};
  test_clip_convex_on_mesh(m, pl, 4, 6, 4);
  test_clip_convex_on_mesh(m, opposite(pl), 7, 12, 7);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip full outside except one edge" << std::endl;
#endif
  pl = {P(1,1,1),P(1,1,0),P(2,0,0)};
  test_clip_convex_on_mesh(m, pl, 2, 0, 0);
#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip full inside except one edge" << std::endl;
#endif
  test_clip_convex_on_mesh(m, opposite(pl), 8, 12, 6);

#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip full outside except one vertex" << std::endl;
#endif
  pl = {P(1,1,1),P(0,3,0),P(3,0,0)};
  test_clip_convex_on_mesh(m, pl, 1, 0, 0);
#ifdef TEST_MESH_KERNEL_VERBOSE
  std::cout << "Test clip full inside except one vertex" << std::endl;
#endif
  test_clip_convex_on_mesh(m, opposite(pl), 8, 12, 6);
}

int main(/*int argc, char** argv*/)
{
  tests<SM, K>();
  tests<ESM, EK>();
  tests<CGAL::Polyhedron_3<K>,K>();
  tests<CGAL::Polyhedron_3<EK>, EK>();
}
