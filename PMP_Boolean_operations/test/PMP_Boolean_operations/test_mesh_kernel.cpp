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

using Surface_mesh = SM;
struct Test_visitor: public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<Surface_mesh>
{
  using halfedge_descriptor = typename boost::graph_traits<Surface_mesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<Surface_mesh>::face_descriptor;
  using vertex_descriptor = typename boost::graph_traits<Surface_mesh>::vertex_descriptor;

  Surface_mesh& sm;
  Surface_mesh::Property_map<Surface_mesh::Face_index, int> fid_map;
  Surface_mesh::Property_map<Surface_mesh::Vertex_index, int> vid_map;
  Surface_mesh::Property_map<Surface_mesh::Halfedge_index, int> hid_map;

  int fid=-1;
  int hid=-1;
  std::size_t nbf=0;
  std::size_t nbe=0;
  std::size_t nb_input_v=sm.number_of_vertices();
  std::size_t nb_input_v_on=0;
  std::size_t nb_new_v_on=0;

  Test_visitor(Surface_mesh& sm)
    : sm(sm)
  {
    bool is_new=false;
    std::tie(fid_map, is_new)=sm.add_property_map<Surface_mesh::Face_index,int>("f:id", -1);
    assert(is_new);
    int i=0; for (auto f : faces(sm)) put(fid_map, f, i++);
    std::tie(vid_map, is_new)=sm.add_property_map<Surface_mesh::Vertex_index,int>("v:id", -1);
    assert(is_new);
    i=0; for (auto v : vertices(sm)) put(vid_map, v, i++);
    std::tie(hid_map, is_new)=sm.add_property_map<Surface_mesh::Halfedge_index,int>("h:id", -1);
    assert(is_new);
    i=0;
    for (auto e : edges(sm))
    {
      auto h = halfedge(e, sm);
      put(hid_map, opposite(h, sm), i);
      put(hid_map, h, i++);
    }
  }

  void before_subface_creations(face_descriptor f_split, const Surface_mesh& pm)
  {
    fid=get(fid_map, f_split);
    assert(fid!=-1);
    assert(&pm==&sm);
  }
  void after_subface_creations(const Surface_mesh& pm)
  {
    assert(&pm==&sm);
  }
  void before_subface_created(const Surface_mesh& pm)
  {
    nbf=sm.number_of_faces();
    assert(&pm==&sm);
  }
  void after_subface_created(face_descriptor f_new, const Surface_mesh& pm)
  {
    assert(&pm==&sm);
    assert(get(fid_map, f_new)==-1);
    put(fid_map, f_new, fid);
    assert(nbf+1==sm.number_of_faces());
  }

  void before_edge_split(halfedge_descriptor h, Surface_mesh& pm)
  {
    hid=get(hid_map, h);
    assert(hid!=-1);
    assert(&pm==&sm);
    nbe=sm.number_of_edges();
  }
  void edge_split(halfedge_descriptor hnew, Surface_mesh& pm)
  {
    assert(&pm==&sm);
    assert(get(hid_map, hnew)==-1);
    assert(hid!=-1);
    put(hid_map, hnew, hid);
    put(hid_map, opposite(hnew, sm), hid);
    assert(nbe+1==sm.number_of_edges());
  }
  void after_edge_split()
  {
    assert(nbe+1==sm.number_of_edges());
  }
  void add_retriangulation_edge(halfedge_descriptor hnew, const Surface_mesh& pm)
  {
    assert(get(hid_map, hnew)==-1);
    put(hid_map, hnew, -2);
    put(hid_map, opposite(hnew, sm), -2);
    assert(&pm==&sm);
  }

  void intersection_point_detected(std::size_t /* i_id */,
                                   int /* sdim */,
                                   halfedge_descriptor h_e,
                                   halfedge_descriptor h_f,
                                   const Surface_mesh& tm_e,
                                   const Surface_mesh& tm_f,
                                   bool is_target_coplanar,
                                   bool is_source_coplanar)
  {
    assert(is_source_coplanar==false);
    assert(h_f==boost::graph_traits<Surface_mesh>::null_halfedge());
    assert(h_e!=boost::graph_traits<Surface_mesh>::null_halfedge());
    if (!is_target_coplanar)
      ++nb_new_v_on;
    else
      ++nb_input_v_on;
    assert(&tm_e==&sm);
    assert(&tm_f==&sm);
  }

  void new_vertex_added(std::size_t id, vertex_descriptor v, const Surface_mesh& pm)
  {
    assert(&pm==&sm);
    assert(get(vid_map, v)==-1);
    put(vid_map, v, nb_input_v+id-nb_input_v_on);
  }

  void before_face_copy(face_descriptor f, const Surface_mesh& src, const Surface_mesh& tgt)
  {
    assert(f==boost::graph_traits<Surface_mesh>::null_face());
    assert(&src==&sm);
    assert(&tgt==&sm);
  }

  void after_face_copy(face_descriptor fsrc, const Surface_mesh& src, face_descriptor ftgt, const Surface_mesh& tgt)
  {
    assert(fsrc==boost::graph_traits<Surface_mesh>::null_face());
    assert(ftgt!=boost::graph_traits<Surface_mesh>::null_face());
    //assert(get(fid_map, ftgt)==-1);
    put(fid_map, ftgt, -2);
    assert(&src==&sm);
    assert(&tgt==&sm);
    for (auto h : halfedges_around_face(halfedge(ftgt, sm), sm))
    {
      if (get(hid_map, h)==-1)
        put(hid_map, h, -2);
    }
  }

  void before_edge_duplicated(halfedge_descriptor h, Surface_mesh& tm)
  {
    hid = get(hid_map, h);
    assert(hid!=-1);
    assert(&tm==&sm);
  }

  void after_edge_duplicated(halfedge_descriptor h, halfedge_descriptor new_hedge, Surface_mesh& tm)
  {
    assert(hid==get(hid_map, h));
    assert(&tm==&sm);
    put(hid_map, new_hedge, hid);
    put(hid_map, opposite(new_hedge, sm), hid);
  }

  void before_vertex_copy(vertex_descriptor v, Surface_mesh& src, Surface_mesh& tgt)
  {
    assert(&src==&sm);
    assert(&tgt==&sm);
    assert(get(vid_map, v)!=-1);
  }
  void after_vertex_copy(vertex_descriptor v, Surface_mesh& src, vertex_descriptor nv, Surface_mesh& tgt)
  {
    assert(&src==&sm);
    assert(&tgt==&sm);
    assert(get(vid_map, v)!=-1);
    put(vid_map, nv, get(vid_map, v));
  }



  void check()
  {
    for (auto f :faces(sm))
    {
      if (get(fid_map, f)==-1) std::cout << sm.point(source(halfedge(f, sm), sm)) << " " << sm.point(target(halfedge(f, sm), sm)) << " " << sm.point(target(next(halfedge(f, sm), sm), sm)) << "\n";
      assert(get(fid_map, f)!=-1);
    }
    for (auto h :halfedges(sm))
    {
      if (get(hid_map, h)==-1)
        std::cout << sm.point(source(h, sm)) << " " << sm.point(target(h,sm)) << "\n";
      assert(get(hid_map, h)!=-1);
    }
    std::size_t nbv_max=nb_input_v+nb_new_v_on;
    for (auto v :vertices(sm))
      assert(get(vid_map, v)!=-1 && get(vid_map, v)<(int)nbv_max);
  }

};

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

  if constexpr(std::is_same_v<Mesh, SM>){
    m_copy = m;
    Test_visitor visitor(m_copy);
    PMP::internal::clip_convex(m_copy, Pl(pl[0], pl[1], pl[2]), CGAL::parameters::visitor(visitor));

    assert(vertices(m_copy).size() == expected_nb_vertices);
    assert(edges(m_copy).size()    == expected_nb_edges);
    assert(faces(m_copy).size()    == expected_nb_faces);
  }

  m_copy = m;
  PMP::internal::clip_convex(m_copy, Pl(pl[0], pl[1], pl[2]));

  assert(vertices(m_copy).size() == expected_nb_vertices);
  assert(edges(m_copy).size()    == expected_nb_edges);
  assert(faces(m_copy).size()    == expected_nb_faces);
}

template<class Mesh, class Plane_3>
void test_clip_convex_with_triangulation_on_mesh(const Mesh &m, const Plane_3 &pl){
  using K =typename CGAL::Kernel_traits<typename Plane_3::value_type>::Kernel;
  using Pl = typename K::Plane_3;

  auto m_copy = m;
  PMP::triangulate_faces(m_copy);
  PMP::internal::clip_convex(m_copy, Pl(pl[0], pl[1], pl[2]));
  PMP::triangulate_faces(m_copy);

  std::size_t expected_nb_vertices = vertices(m_copy).size();
  std::size_t expected_nb_edges    = edges(m_copy).size();
  std::size_t expected_nb_faces    = faces(m_copy).size();

  m_copy = m;
  PMP::triangulate_faces(m_copy);
  PMP::internal::clip_convex(m_copy, Pl(pl[0], pl[1], pl[2]), CGAL::parameters::do_not_triangulate_faces(false));

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
