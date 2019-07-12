#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <iostream>
#include <fstream>
#include <sstream>

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
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // test with a plane
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

  // clipping with identity
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(true)
                      .use_compact_clipper(true)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(num_vertices(tm1)==8);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(CGAL::is_empty(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(true)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(num_vertices(tm1)==8);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(CGAL::is_empty(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(true)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==4);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // test orientation + patch without input vertex
  CGAL::make_tetrahedron(
    K::Point_3(0.53, -1.3, 0.2),
    K::Point_3(0.53, 1.1, 0.2),
    K::Point_3(0.53, -1.3, 0.4),
    K::Point_3(0.73, -1.3, 0.2),
    tm2);
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==6);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  CGAL::make_tetrahedron(
    K::Point_3(0.53, -1.3, 0.2),
    K::Point_3(0.53, 1.1, 0.2),
    K::Point_3(0.53, -1.3, 0.4),
    K::Point_3(0.73, -1.3, 0.2),
    tm2);
  PMP::reverse_face_orientations(tm2);
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==6+8);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // clip meshes with intersection polyline opened
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(1, 0, 0, -2));
  assert(vertices(tm1).size()==4);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(-1, 0, 0, 2));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);

  // test with clipper on border edge
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, 1, 0 , 0));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, -1, 0 , 0));
  assert(vertices(tm1).size()==4);
  CGAL::clear(tm1);

  // test with clipper on border edge: full triangle
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(true));
  assert(vertices(tm1).size()!=0);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(false));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);

  // test tangencies
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(1, 0, 0, -1));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(-1, 0, 0, 1));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2, params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                      params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::reverse_face_orientations(tm2);
  PMP::clip(tm1, tm2, params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                      params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // test special case
  input.open("data-clip/tm_1.off");
  input >> tm1;
  input.close();
  input.open("data-clip/clipper_1.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2, params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                      params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(is_valid_polygon_mesh(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // non-manifold border vertices
  std::stringstream ss;
  ss << "OFF\n 5 2 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 0 1 4\n3 1 2 3\n";
  ss >> tm1;
  PMP::clip(tm1, K::Plane_3(-1,0,0,2));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);

  ss = std::stringstream();
  ss << "OFF\n 7 4 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 3 5\n";
  ss >> tm1;
  CGAL::Euler::remove_face(halfedge(*CGAL::cpp11::prev(faces(tm1).end()),tm1),tm1);
  PMP::clip(tm1, K::Plane_3(-1,0,0,2));
  assert(vertices(tm1).size()==6);
  CGAL::clear(tm1);

  ss = std::stringstream();
  ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
  ss >> tm1;
  for (int i=0;i<3;++i)
    CGAL::Euler::remove_face(halfedge(*CGAL::cpp11::prev(faces(tm1).end()),tm1),tm1);
  PMP::clip(tm1, K::Plane_3(-1,0,0,2));
  assert(vertices(tm1).size()==7);
  CGAL::clear(tm1);

  ss = std::stringstream();
  ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
  ss >> tm1;
  for (int i=0;i<3;++i)
    CGAL::Euler::remove_face(halfedge(*CGAL::cpp11::prev(faces(tm1).end()),tm1),tm1);
  PMP::clip(tm1, K::Plane_3(0,1,0,0));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);

  ss = std::stringstream();
  ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
  ss >> tm1;
  for (int i=0;i<3;++i)
    CGAL::Euler::remove_face(halfedge(*CGAL::cpp11::prev(faces(tm1).end()),tm1),tm1);
  PMP::clip(tm1, K::Plane_3(0,-1,0,0));
  assert(vertices(tm1).size()==7);
  CGAL::clear(tm1);
}

int main()
{
  test<Surface_mesh>();
  test<Polyhedron>();

  return 0;
}
