#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template <class TriangleMesh>
void test()
{
  // test with a clipper mesh
  TriangleMesh tm1, tm2;
  std::ifstream input("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;
  input.close();

  PMP::clip(tm1, tm2,
            params::clip_volume(false)
              .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
            params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2))
  );
  assert(!CGAL::is_closed(tm1));

  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;
  input.close();

  PMP::clip(tm1, tm2, params::clip_volume(true)
              .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
            params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(CGAL::is_closed(tm1));

  // test with a plane
  CGAL::clear(tm1);
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();

  K::Plane_3 plane(0, 0, 1, -1);

  PMP::clip(tm1, plane, params::clip_volume(true));
  assert(CGAL::is_closed(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, plane, params::clip_volume(false)
              .use_compact_clipper(false));
  assert(!CGAL::is_closed(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, plane, params::clip_volume(false)
                        .use_compact_clipper(true));
  assert(CGAL::is_closed(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, K::Plane_3(-0.236474, 0.437732, 0.867451, -0.838791), params::clip_volume(true));
  assert(CGAL::is_closed(tm1));
  assert(!CGAL::is_empty(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, K::Plane_3(0, 0, 1, 2));
  assert(CGAL::is_empty(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, K::Plane_3(0, 0, 1, -2));
  assert(!CGAL::is_empty(tm1));
  CGAL::clear(tm1);
}

int main()
{
  test<Surface_mesh>();
  test<Polyhedron>();

  return 0;
}
