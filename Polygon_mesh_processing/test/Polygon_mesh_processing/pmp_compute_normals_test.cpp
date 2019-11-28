#define CGAL_PMP_COMPUTE_NORMAL_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel         EPECK;

typedef CGAL::Surface_mesh<EPICK::Point_3>                        EPICK_SM;
typedef CGAL::Surface_mesh<EPECK::Point_3>                        EPECK_SM;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename K, typename Mesh, typename VertexNormalPmap, typename FaceNormalPmap>
void test(const Mesh& mesh,
          const VertexNormalPmap& vnormals,
          const FaceNormalPmap& fnormals)
{
  typedef typename K::Vector_3                                            Vector;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor             face_descriptor;

  typedef typename CGAL::GetVertexPointMap<Mesh>::const_type              VPMap;
  VPMap vpmap = get_const_property_map(CGAL::vertex_point, mesh);

  const vertex_descriptor first_vertex = *(vertices(mesh).begin());
  const face_descriptor first_face = *(faces(mesh).begin());

  PMP::compute_face_normals(mesh, fnormals);
  PMP::compute_face_normals(mesh, fnormals, PMP::parameters::vertex_point_map(vpmap));
  PMP::compute_face_normals(mesh, fnormals, CGAL::parameters::vertex_point_map(vpmap)
                                                             .geom_traits(K()));

  Vector f0n = PMP::compute_face_normal(first_face, mesh);
  assert(f0n == get(fnormals, first_face));
  PMP::compute_face_normal(first_face, mesh, PMP::parameters::vertex_point_map(vpmap));

  PMP::compute_vertex_normals(mesh, vnormals);
  PMP::compute_vertex_normals(mesh, vnormals, PMP::parameters::vertex_point_map(vpmap));
  PMP::compute_vertex_normals(mesh, vnormals, CGAL::parameters::vertex_point_map(vpmap)
                                                               .geom_traits(K()));

  Vector v0n = PMP::compute_vertex_normal(first_vertex, mesh);
  assert(v0n == get(vnormals, first_vertex));
  v0n = PMP::compute_vertex_normal(first_vertex, mesh, PMP::parameters::vertex_point_map(vpmap)
                                                                       .face_normal_map(fnormals));
  std::cout.precision(17);
  std::cout << v0n << " versus " << get(vnormals, first_vertex) << std::endl;
  assert(v0n == get(vnormals, first_vertex));

  PMP::compute_normals(mesh, vnormals, fnormals);
  PMP::compute_normals(mesh, vnormals, fnormals, PMP::parameters::vertex_point_map(vpmap));
  PMP::compute_normals(mesh, vnormals, fnormals, CGAL::parameters::vertex_point_map(vpmap)
                                                                  .geom_traits(K()));

  for(vertex_descriptor v : vertices(mesh)) {
    assert(get(vnormals, v) != CGAL::NULL_VECTOR);
  }

  for(face_descriptor f : faces(mesh)) {
    assert(get(fnormals, f) != CGAL::NULL_VECTOR);
  }
}

template<typename K>
void test_SM(const char* file_name)
{
  typedef CGAL::Surface_mesh<typename K::Point_3>                         SM;
  typedef typename boost::graph_traits<SM>::vertex_descriptor             vertex_descriptor;
  typedef typename boost::graph_traits<SM>::face_descriptor               face_descriptor;

  typedef typename K::Vector_3                                            Vector;

  SM mesh;
  std::ifstream input(file_name);
  if(!(input >> mesh))
  {
    std::cerr << "Error: cannot read " << file_name << " as a Surface_mesh\n";
    assert(false);
  }

  typename SM::template Property_map<vertex_descriptor, Vector> vnormals;
  vnormals = mesh.template add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
  typename SM::template Property_map<face_descriptor, Vector> fnormals;
  fnormals = mesh.template add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;

  test<K>(mesh, vnormals, fnormals);
}

template<typename K>
void test_Polyhedron(const char* file_name)
{
  typedef CGAL::Polyhedron_3<K>                                           Polyhedron;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor       face_descriptor;

  typedef typename K::Vector_3                                            Vector;

  Polyhedron mesh;
  std::ifstream input(file_name);
  if(!(input >> mesh))
  {
    std::cerr << "Error: cannot read " << file_name << " as a Polyhedron\n";
    assert(false);
  }

  typedef std::map<vertex_descriptor, Vector>                             VertexNormalMap;
  typedef boost::associative_property_map<VertexNormalMap>                VertexNormalPMap;
  typedef std::map<face_descriptor, Vector>                               FaceNormalMap;
  typedef boost::associative_property_map<FaceNormalMap>                  FaceNormalPMap;

  VertexNormalMap vn_map;
  FaceNormalMap fn_map;

  for(vertex_descriptor v : vertices(mesh))
    vn_map[v] = CGAL::NULL_VECTOR;
  for(face_descriptor f : faces(mesh))
    fn_map[f] = CGAL::NULL_VECTOR;

  VertexNormalPMap vnormals(vn_map);
  FaceNormalPMap fnormals(fn_map);

  test<K>(mesh, vnormals, fnormals);
}

void test(const char* filename)
{
  std::cout << "test " << filename << "..." << std::endl;

  test_SM<EPICK>(filename);
  test_SM<EPECK>(filename);
  test_Polyhedron<EPICK>(filename);
  test_Polyhedron<EPECK>(filename);
}

int main()
{
  CGAL::Set_ieee_double_precision pfr;

  test("data/elephant.off");
  test("data/folded_star.off");
  test("data/joint_refined.off");
  test("data/mannequin-devil.off");
  test("data/U.off");

  std::cerr << "All done." << std::endl;
  return EXIT_SUCCESS;
}
