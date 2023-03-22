// #define CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/centroid.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK;

typedef CGAL::Surface_mesh<EPICK::Point_3>                           EPICK_SM;
typedef CGAL::Surface_mesh<EPECK::Point_3>                           EPECK_SM;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename K, typename Mesh, typename VertexNormalPmap, typename FaceNormalPmap>
void test(const Mesh& mesh,
          const VertexNormalPmap& vnormals,
          const FaceNormalPmap& fnormals)
{
  typedef typename K::Vector_3                                            Vector;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor             face_descriptor;

  typedef typename CGAL::GetVertexPointMap<Mesh>::const_type              VPMap;
  VPMap vpmap = get_const_property_map(CGAL::vertex_point, mesh);

  assert(!is_empty(mesh));

  const vertex_descriptor first_vertex = *(vertices(mesh).begin());
  const face_descriptor first_face = *(faces(mesh).begin());

  PMP::compute_face_normals(mesh, fnormals);
  PMP::compute_face_normals(mesh, fnormals, CGAL::parameters::vertex_point_map(vpmap));
  PMP::compute_face_normals(mesh, fnormals, CGAL::parameters::vertex_point_map(vpmap)
                                                             .geom_traits(K()));

  Vector f0n = PMP::compute_face_normal(first_face, mesh);
  assert(f0n == get(fnormals, first_face));
  PMP::compute_face_normal(first_face, mesh, CGAL::parameters::vertex_point_map(vpmap));

  PMP::compute_vertex_normals(mesh, vnormals);
  PMP::compute_vertex_normals(mesh, vnormals, CGAL::parameters::vertex_point_map(vpmap));
  PMP::compute_vertex_normals(mesh, vnormals, CGAL::parameters::vertex_point_map(vpmap)
                                                               .geom_traits(K()));

  Vector v0n = PMP::compute_vertex_normal(first_vertex, mesh);
  assert(v0n == get(vnormals, first_vertex));
  v0n = PMP::compute_vertex_normal(first_vertex, mesh, CGAL::parameters::vertex_point_map(vpmap)
                                                                        .face_normal_map(fnormals));
  std::cout.precision(17);
  assert(v0n == get(vnormals, first_vertex));

  PMP::compute_normals(mesh, vnormals, fnormals);
  PMP::compute_normals(mesh, vnormals, fnormals, CGAL::parameters::vertex_point_map(vpmap));
  PMP::compute_normals(mesh, vnormals, fnormals, CGAL::parameters::vertex_point_map(vpmap)
                                                                  .geom_traits(K()));

#if 1//def CGAL_PMP_COMPUTE_NORMAL_DEBUG_PP
  std::ofstream vn_out("vertex_normals.cgal.polylines.txt");
  std::ofstream fn_out("face_normals.cgal.polylines.txt");

  const CGAL::Bbox_3 bb = PMP::bbox(mesh);
  const auto bbox_diagonal = CGAL::sqrt(CGAL::square(bb.xmax() - bb.xmin()) +
                                        CGAL::square(bb.ymax() - bb.ymin()) +
                                        CGAL::square(bb.zmax() - bb.zmin()));

  for(vertex_descriptor v : vertices(mesh))
    vn_out << "2 " << get(vpmap, v) << " " << get(vpmap, v) + 0.1 * bbox_diagonal * get(vnormals, v) << "\n";
  vn_out.close();

  for(face_descriptor f : faces(mesh))
  {
    std::list<typename K::Point_3> vertices;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh))
      vertices.push_back(get(vpmap, target(h, mesh)));

    const typename K::Point_3& c = CGAL::centroid(vertices.begin(), vertices.end());
    fn_out << "2 " << c << " " << c + 0.1 * bbox_diagonal * get(fnormals, f) << "\n";
  }
  fn_out.close();
#endif

  // Check sanity of output
  for(face_descriptor f : faces(mesh))
  {
    // tests on non triangular meshes are @todo
    if(CGAL::is_triangle(halfedge(f, mesh), mesh))
    {
      if (PMP::is_degenerate_triangle_face(f, mesh))
      {
        if (std::is_same<K, EPECK>())
          assert(get(fnormals, f) == CGAL::NULL_VECTOR);
      }
      else
        assert(get(fnormals, f) != CGAL::NULL_VECTOR);
    }
  }

  for(vertex_descriptor v : vertices(mesh))
  {
    if(get(vnormals, v) == CGAL::NULL_VECTOR)
    {
      for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(v, mesh), mesh))
      {
        if(!is_border(h, mesh))
        {
          // There are other cases where a vertex normal can be null without the face normals being null,
          // (only two incident faces with opposite normals, for example), but they are not tested for now.
          assert(get(fnormals, face(h, mesh)) == CGAL::NULL_VECTOR);
        }
      }
    }
  }
}

template<typename K>
void test_SM(const std::string file_name)
{
  typedef CGAL::Surface_mesh<typename K::Point_3>                         SM;
  typedef typename boost::graph_traits<SM>::vertex_descriptor             vertex_descriptor;
  typedef typename boost::graph_traits<SM>::face_descriptor               face_descriptor;

  typedef typename K::Vector_3                                            Vector;

  std::cout << "Test with Surface_mesh, and kernel: " << typeid(K).name() << std::endl;

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
void test_Polyhedron(const std::string file_name)
{
  typedef CGAL::Polyhedron_3<K>                                           Polyhedron;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor       face_descriptor;

  typedef typename K::Vector_3                                            Vector;

  std::cout << "Test with Polyhedron, and kernel: " << typeid(K).name() << std::endl;

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

void test(const std::string filename)
{
  std::cout << "test " << filename << "..." << std::endl;

  // EPECK disabled because it takes too long for the testsuite due to sq roots comparisons,
  // but it passed.
  test_SM<EPICK>(filename);
//  test_SM<EPECK>(filename);
  test_Polyhedron<EPICK>(filename);
//  test_Polyhedron<EPECK>(filename);
}

int main()
{
  std::cout.precision(17);

  CGAL::Set_ieee_double_precision pfr;

  test(CGAL::data_file_path("meshes/elephant.off"));
  test("data/folded_star.off");
  test("data/joint_refined.off");
  test(CGAL::data_file_path("meshes/mannequin-devil.off"));
  test("data/U.off");

  test("data_degeneracies/deg_on_border.off");
  test("data_degeneracies/degtri_edge.off");
  test("data_degeneracies/degtri_three.off");
  test("data_degeneracies/degtri_four.off");
  test("data_degeneracies/degtri_nullface.off");
  test("data_degeneracies/degtri_single.off");
  test("data_degeneracies/existing_flip.off");
  test("data_degeneracies/fused_vertices.off");
  test("data_degeneracies/small_ccs.off");
  test("data_degeneracies/trihole.off");

  std::cerr << "All done." << std::endl;
  return EXIT_SUCCESS;
}
