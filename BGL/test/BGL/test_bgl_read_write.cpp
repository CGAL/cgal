#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#if defined(CGAL_USE_OPENMESH)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#endif

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Origin.h>

#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>                                                   Kernel;

typedef CGAL::Exact_predicates_inexact_constructions_kernel                              EPICK;

typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>                     Polyhedron;

typedef CGAL::Surface_mesh<Kernel::Point_3>                                              SM;

typedef CGAL::Linear_cell_complex_traits<3, Kernel>                                      MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper<2, 3, MyTraits>::type LCC;

#if defined(CGAL_USE_OPENMESH)

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/>                                   OMesh;

#endif

const double epsilon=5e-4;
template < typename Type >
bool approx_equal_nt(const Type &t1, const Type &t2)
{
      if (t1 == t2)
              return true;
      if (CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < std::abs(t1)*epsilon)
              return true;
      std::cout << " Approximate comparison failed between : " << t1 << "  and  " << t2 << std::endl;
      return false;
}


template <typename Mesh, typename VPM1, typename VPM2, typename VIM1, typename VIM2>
bool are_equal_meshes(const Mesh& fg1, const VPM1 vpm1, const Mesh& fg2, const VPM2 vpm2,
                      const VIM1 vim1, const VIM2 vim2)
{
  typedef typename boost::property_traits<VPM1>::value_type P;

  if(num_vertices(fg1) != num_vertices(fg2) ||
     num_halfedges(fg1) != num_halfedges(fg2) ||
     num_edges(fg1) != num_edges(fg2) ||
     num_faces(fg1) != num_faces(fg2))
    return false;
  std::vector<P> fg1_points, fg2_points;
  for(auto v : vertices(fg1))
  {
    fg1_points.push_back(get(vpm1, v));
  }

  for(auto v : vertices(fg2))
    fg2_points.push_back(get(vpm2, v));

  for(std::size_t id = 0; id < fg1_points.size(); ++id)
  {
    P p1 = fg1_points[id];
    P p2 = fg2_points[id];
    if (! (approx_equal_nt(p1.x(), p2.x())
           && approx_equal_nt(p1.y(), p2.y())
           && approx_equal_nt(p1.z(), p2.z())))
      return false;
  }
  fg1_points.clear();
  fg1_points.shrink_to_fit();
  fg2_points.clear();
  fg2_points.shrink_to_fit();

  std::vector<std::size_t> fg1_faces;
  std::vector<std::size_t> fg2_faces;
  for(auto f : faces(fg1))
    for(auto v : CGAL::vertices_around_face(halfedge(f, fg1), fg1))
      fg1_faces.push_back(get(vim1, v));

  for(auto f : faces(fg2))
    for(auto v : CGAL::vertices_around_face(halfedge(f, fg2), fg2))
      fg2_faces.push_back(get(vim2, v));

  if(fg1_points != fg2_points)
    return false;

  return true;
}

template <typename Mesh>
bool are_equal_meshes(const Mesh& fg1, const Mesh& fg2)
{
  typedef typename CGAL::GetInitializedVertexIndexMap<Mesh>::const_type      VIM;
  VIM vim1 = CGAL::get_initialized_vertex_index_map(fg1);
  VIM vim2 = CGAL::get_initialized_vertex_index_map(fg2);
  return are_equal_meshes(fg1, get(CGAL::vertex_point, fg1), fg2, get(CGAL::vertex_point, fg2), vim1, vim2);
}

template<typename Mesh, typename K>
void test_bgl_OFF(const char* filename)
{
  // read with OFF
  Mesh fg;
  std::ifstream is(filename);
  bool ok = CGAL::IO::read_OFF(is, fg);
  assert(ok);
  assert(num_vertices(fg) != 0 && num_faces(fg) != 0);
  is.close();
  fg.clear();

  is.open(filename, std::ios::binary);
  ok = CGAL::IO::read_OFF(is, fg);
  assert(ok);
  assert(num_vertices(fg) != 0 && num_faces(fg) != 0);
  // write with OFF
  {
    std::ofstream os("tmp.off");
    ok = CGAL::IO::write_OFF(os, fg);
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_OFF("tmp.off", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // write with PM
  {
    ok = CGAL::IO::write_polygon_mesh("tmp.obj.off", fg);
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_polygon_mesh("tmp.obj.off", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }


  // Test [STCN]OFF
  typedef typename K::Point_2                                                                Point_2;
  typedef typename K::Vector_3                                                               Vector;
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Vector> >::type VertexNormalMap;
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Point_2> >::type VertexTextureMap;
  typedef typename boost::property_map<Mesh, CGAL::dynamic_face_property_t<CGAL::IO::Color> >::type FaceColorMap;

  // COFF
  {
    CGAL::clear(fg);
    VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg);
    FaceColorMap fcm = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg);

    ok = CGAL::IO::read_OFF("data/mesh_with_colors.off", fg, CGAL::parameters::vertex_color_map(vcm)
                                                                              .face_color_map(fcm));
    assert(ok);
    assert(num_vertices(fg) == 8 && num_faces(fg) == 4);

    for(auto v : vertices(fg))
      assert(get(vcm, v) != CGAL::IO::Color());

    for(auto f : faces(fg))
      assert(get(fcm, f) != CGAL::IO::Color());

    // write with OFF
    {
      ok = CGAL::IO::write_OFF("tmp.off", fg, CGAL::parameters::vertex_color_map(vcm)
                                                               .face_color_map(fcm));
      assert(ok);

      Mesh fg2;
      VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg2);
      FaceColorMap fcm2 = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg2);

      ok = CGAL::IO::read_polygon_mesh("tmp.off", fg2, CGAL::parameters::vertex_color_map(vcm2)
                                                                        .face_color_map(fcm2));
      assert(ok);
      assert(are_equal_meshes(fg, fg2));

      for(auto v : vertices(fg2))
        assert(get(vcm2, v) != CGAL::IO::Color());

      for(auto f : faces(fg2))
        assert(get(fcm2, f) != CGAL::IO::Color());
    }

    // write with PM
    {
      ok = CGAL::IO::write_polygon_mesh("tmp.off", fg, CGAL::parameters::vertex_color_map(vcm));
      assert(ok);

      Mesh fg2;
      VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg2);

      ok = CGAL::IO::read_polygon_mesh("tmp.off", fg2, CGAL::parameters::vertex_color_map(vcm2));
      assert(ok);
      assert(are_equal_meshes(fg, fg2));

      for(auto v : vertices(fg2))
        assert(get(vcm2, v) != CGAL::IO::Color());
    }
  }

  // NOFF
  {
    CGAL::clear(fg);
    VertexNormalMap vnm = get(CGAL::dynamic_vertex_property_t<Vector>(), fg);

    ok = CGAL::IO::read_OFF("data/mesh_with_normals.off", fg, CGAL::parameters::vertex_normal_map(vnm));
    assert(ok);

    for(auto v : vertices(fg))
      assert(get(vnm, v) != CGAL::NULL_VECTOR);

    // write with OFF
    {
      std::ofstream os("tmp.off");
      ok = CGAL::IO::write_OFF("tmp.off", fg, CGAL::parameters::vertex_normal_map(vnm));
      assert(ok);

      Mesh fg2;
      VertexNormalMap vnm2 = get(CGAL::dynamic_vertex_property_t<Vector>(), fg2);

      ok = CGAL::IO::read_polygon_mesh("tmp.off", fg2, CGAL::parameters::vertex_normal_map(vnm2));
      assert(ok);
      assert(are_equal_meshes(fg, fg2));

      for(auto v : vertices(fg2))
        assert(get(vnm2, v) != CGAL::NULL_VECTOR);
    }

    // write with PM
    {
      ok = CGAL::IO::write_polygon_mesh("tmp.off", fg, CGAL::parameters::vertex_normal_map(vnm));
      assert(ok);

      Mesh fg2;
      VertexNormalMap vnm2 = get(CGAL::dynamic_vertex_property_t<Vector>(), fg2);

      ok = CGAL::IO::read_polygon_mesh("tmp.off", fg2, CGAL::parameters::vertex_normal_map(vnm2));
      assert(ok);
      assert(are_equal_meshes(fg, fg2));

      for(auto v : vertices(fg2))
        assert(get(vnm2, v) != CGAL::NULL_VECTOR);
    }
  }

  // STCNOFF
  {
    CGAL::clear(fg);
    std::ifstream is("data/full.off", std::ios::binary);

    VertexNormalMap vnm2 = get(CGAL::dynamic_vertex_property_t<Vector>(), fg);
    VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg);
    VertexTextureMap vtm2 = get(CGAL::dynamic_vertex_property_t<Point_2>(), fg);
    FaceColorMap fcm2 = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg);

    ok = CGAL::IO::read_OFF(is, fg, CGAL::parameters::vertex_normal_map(vnm2)
                                                     .vertex_color_map(vcm2)
                                                     .vertex_texture_map(vtm2)
                                                     .face_color_map(fcm2));
    assert(ok);
    assert(num_vertices(fg) != 0 && num_faces(fg) != 0);

    for(auto v : vertices(fg))
    {
      assert(get(vnm2, v) != CGAL::NULL_VECTOR);
      assert(get(vcm2, v) != CGAL::IO::Color());
    }

    for(auto f : faces(fg))
      assert(get(fcm2, f) != CGAL::IO::Color());
    fg.clear();
    is.close();

    is.open("data/full.off");
    VertexNormalMap vnm = get(CGAL::dynamic_vertex_property_t<Vector>(), fg);
    VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg);
    VertexTextureMap vtm = get(CGAL::dynamic_vertex_property_t<Point_2>(), fg);
    FaceColorMap fcm = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg);
    ok = CGAL::IO::read_OFF(is, fg, CGAL::parameters::vertex_normal_map(vnm)
                                                     .vertex_color_map(vcm)
                                                     .vertex_texture_map(vtm)
                                                     .face_color_map(fcm));
    assert(ok);
    assert(num_vertices(fg) != 0 && num_faces(fg) != 0);

    for(auto v : vertices(fg))
    {
      assert(get(vnm, v) != CGAL::NULL_VECTOR);
      assert(get(vcm, v) != CGAL::IO::Color());
    }
    // write with OFF
    {
      std::ofstream os("tmp.off");
      ok = CGAL::IO::write_OFF("tmp.off", fg, CGAL::parameters::vertex_normal_map(vnm)
                                                               .vertex_color_map(vcm)
                                                               .vertex_texture_map(vtm)
                                                               .face_color_map(fcm));
      assert(ok);

      Mesh fg2;
      VertexNormalMap vnm2 = get(CGAL::dynamic_vertex_property_t<Vector>(), fg2);
      VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg2);
      VertexTextureMap vtm2 = get(CGAL::dynamic_vertex_property_t<Point_2>(), fg2);
      FaceColorMap fcm2 = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg2);

      ok = CGAL::IO::read_polygon_mesh("tmp.off", fg2, CGAL::parameters::vertex_normal_map(vnm2)
                                                                        .vertex_color_map(vcm2)
                                                                        .vertex_texture_map(vtm2)
                                                                        .face_color_map(fcm2));
      assert(ok);
      assert(are_equal_meshes(fg, fg2));

      for(auto v : vertices(fg2))
      {
        assert(get(vnm2, v) != CGAL::NULL_VECTOR);
        assert(get(vcm2, v) != CGAL::IO::Color());
      }

      for(auto f : faces(fg2))
        assert(get(fcm2, f) != CGAL::IO::Color());
    }

    // write with PM
    {
      ok = CGAL::IO::write_polygon_mesh("tmp.off", fg, CGAL::parameters::vertex_normal_map(vnm)
                                                                        .vertex_color_map(vcm)
                                                                        .vertex_texture_map(vtm)
                                                                        .face_color_map(fcm));
      assert(ok);

      Mesh fg2;
      VertexNormalMap vnm2 = get(CGAL::dynamic_vertex_property_t<Vector>(), fg2);
      VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg2);
      VertexTextureMap vtm2 = get(CGAL::dynamic_vertex_property_t<Point_2>(), fg2);
      FaceColorMap fcm2 = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg2);

      ok = CGAL::IO::read_polygon_mesh("tmp.off", fg2, CGAL::parameters::vertex_normal_map(vnm2)
                                                                        .vertex_color_map(vcm2)
                                                                        .vertex_texture_map(vtm2)
                                                                        .face_color_map(fcm2));
      assert(ok);
      assert(are_equal_meshes(fg, fg2));

      for(auto v : vertices(fg2))
      {
        assert(get(vnm2, v) != CGAL::NULL_VECTOR);
        assert(get(vcm2, v) != CGAL::IO::Color());
      }

      for(auto f : faces(fg2))
        assert(get(fcm2, f) != CGAL::IO::Color());
    }
  }
  //@todo test multi objects in a single file

  // test wrong inputs
  std::cerr << " ########### Error text is expected to follow." << std::endl;
  ok = CGAL::IO::read_OFF("data/mesh_that_doesnt_exist.off", fg);
  assert(!ok);
  ok = CGAL::IO::read_OFF("data/invalid_cut.off", fg); // cut in half
  assert(!ok);
  ok = CGAL::IO::read_OFF("data/invalid_nv.off", fg); // wrong number of points
  assert(!ok);
  ok = CGAL::IO::read_OFF("data/sphere.obj", fg);
  assert(!ok);
  ok = CGAL::IO::read_OFF("data/pig.stl", fg);
  assert(!ok);
  std::cerr << " ########### No more error text from here." << std::endl;
}

template<typename Mesh, typename K>
void test_bgl_OBJ(const std::string filename)
{
  std::cout << "Test OBJ; input = " << filename << " kernel = " << typeid(K).name() << std::endl;

  Mesh fg;

  std::ifstream is(filename);
  bool ok = CGAL::IO::read_OBJ(is, fg, CGAL::parameters::verbose(true));
  assert(ok);
  assert(filename != "data/sphere.obj" || (num_vertices(fg) == 272 && num_faces(fg) == 540));

  // write with OBJ
  {
    std::ofstream os("tmp.obj");
    ok = CGAL::IO::write_OBJ(os, fg);
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_OBJ("tmp.obj", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // write with PM
  {
    ok = CGAL::IO::write_polygon_mesh("tmp.obj", fg);
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_polygon_mesh("tmp.obj", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // Test NPs
  CGAL::clear(fg);
  ok = CGAL::IO::read_OBJ("data/sphere.obj", fg);
  assert(ok);
  assert(num_vertices(fg) == 272 && num_faces(fg) == 540);

  // write with OBJ
  {
    std::ofstream os("tmp.obj");
    ok = CGAL::IO::write_OBJ("tmp.obj", fg);
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_polygon_mesh("tmp.obj", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // write with PM
  {
    ok = CGAL::IO::write_polygon_mesh("tmp.obj", fg);
    assert(ok);

    Mesh fg2;

    ok = CGAL::IO::read_polygon_mesh("tmp.obj", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // test wrong inputs
  std::cerr << " ########### Error text is expected to follow." << std::endl;
  ok = CGAL::IO::read_OBJ("data/mesh_that_doesnt_exist.obj", fg);
  assert(!ok);
  ok = CGAL::IO::read_OBJ("data/invalid_cut.obj", fg); // invalid vertex ids
  assert(!ok);
  ok = CGAL::IO::read_OBJ("data/genus3.off", fg); // wrong extension
  assert(!ok);
  ok = CGAL::IO::read_OBJ("data/pig.stl", fg);
  assert(!ok);
  std::cerr << " ########### No more error text from here." << std::endl;
}

template<class Mesh>
void test_bgl_PLY(const std::string filename,
                  bool binary = false)
{
  std::cout << "Test PLY; filename = " << filename << " Mesh = " << typeid(Mesh).name() << " binary = " << binary << std::endl;

  Mesh fg;
  std::ifstream is(filename);
  bool ok = CGAL::IO::read_PLY(is, fg, CGAL::parameters::use_binary_mode(false));
  is.close();
  assert(ok);
  assert(filename != "data/colored_tetra.ply" || (num_vertices(fg) == 4 && num_faces(fg) == 4));
   if(!binary)
   {
     fg.clear();
     is.open(filename, std::ios::binary);
     bool ok = CGAL::IO::read_PLY(is, fg, CGAL::parameters::use_binary_mode(false));
     is.close();
     assert(ok);
     assert(filename != "data/colored_tetra.ply" || (num_vertices(fg) == 4 && num_faces(fg) == 4));
   }

  // write with PLY
  {
    std::ofstream os;
    if(binary)
    {
      os.open("tmp.ply", std::ios::binary);
      CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    }
    else
    {
      os.open("tmp.ply");
    }

    ok = CGAL::IO::write_PLY(os, fg);
    assert(ok);

    ok = CGAL::IO::write_PLY(os, fg, "test");
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_PLY("tmp.ply", fg2);

    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // test NPs
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<CGAL::IO::Color> >::type VertexColorMap;
  typedef typename boost::property_map<Mesh, CGAL::dynamic_face_property_t<CGAL::IO::Color> >::type FaceColorMap;

  CGAL::clear(fg);
  VertexColorMap vcm = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg);
  FaceColorMap fcm = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg);

  std::ifstream is_c("data/colored_tetra.ply"); // ASCII
  ok = CGAL::IO::read_PLY(is_c, fg, CGAL::parameters::vertex_color_map(vcm)
                                                    .face_color_map(fcm));
  assert(ok);
  assert(num_vertices(fg) == 4 && num_faces(fg) == 4);

  for(auto v : vertices(fg))
  {
    assert(get(vcm, v) != CGAL::IO::Color());
  }

  for(auto f : faces(fg))
    assert(get(fcm, f) != CGAL::IO::Color());

  // write with PLY
  {
    ok = CGAL::IO::write_PLY("tmp.ply", fg, CGAL::parameters::vertex_color_map(vcm)
                                                             .face_color_map(fcm)
                                                             .use_binary_mode(binary));
    assert(ok);

    Mesh fg2;
    VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg2);
    FaceColorMap fcm2 = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg2);

    std::ifstream is_rpm;
    if(binary)
    {
      is_rpm.open("tmp.ply", std::ios::binary);
      CGAL::IO::set_mode(is_rpm, CGAL::IO::BINARY);
    }
    else
    {
      is_rpm.open("tmp.ply");
    }
    ok = CGAL::IO::read_PLY(is_rpm, fg2, CGAL::parameters::vertex_color_map(vcm2)
                                                          .face_color_map(fcm2));
    assert(ok);
    assert(are_equal_meshes(fg, fg2));

    // @tmp
//    for(auto v : vertices(fg2))
//      assert(get(vcm2, v) != CGAL::IO::Color());

//    for(auto f : faces(fg2))
//      assert(get(fcm2, f) != CGAL::IO::Color());
  }

  // write with PM
  {
    ok = CGAL::IO::write_polygon_mesh("tmp.ply", fg, CGAL::parameters::vertex_color_map(vcm)
                                                                      .face_color_map(fcm)
                                                                      .use_binary_mode(binary));
    assert(ok);

    Mesh fg2;
    VertexColorMap vcm2 = get(CGAL::dynamic_vertex_property_t<CGAL::IO::Color>(), fg2);
    FaceColorMap fcm2 = get(CGAL::dynamic_face_property_t<CGAL::IO::Color>(), fg2);

    ok = CGAL::IO::read_polygon_mesh("tmp.ply", fg2, CGAL::parameters::vertex_color_map(vcm2)
                                                                      .face_color_map(fcm2)
                                                                      .use_binary_mode(binary));
    assert(ok);
    assert(are_equal_meshes(fg, fg2));

    // @tmp
//    for(auto v : vertices(fg2))
//      assert(get(vcm2, v) != CGAL::IO::Color());

//    for(auto f : faces(fg2))
//      assert(get(fcm2, f) != CGAL::IO::Color());
  }

  // test wrong inputs
  std::cerr << " ########### Error text is expected to follow." << std::endl;
  ok = CGAL::IO::read_PLY("data/mesh_that_doesnt_exist.ply", fg);
  assert(!ok);
  ok = CGAL::IO::read_PLY("data/invalid_cut.ply", fg); // cut in half
  assert(!ok);
  ok = CGAL::IO::read_PLY("data/invalid_nv.ply", fg); // broken formatting
  assert(!ok);
  ok = CGAL::IO::read_PLY("data/binary_cut.ply", fg); // broken binary
  assert(!ok);
  ok = CGAL::IO::read_PLY("data/cube.off", fg);
  assert(!ok);
  ok = CGAL::IO::read_PLY("data/pig.stl", fg);
  assert(!ok);
  std::cerr << " ########### No more error text from here." << std::endl;
}

template<class Mesh>
struct Custom_VPM
{
  typedef Custom_VPM<Mesh>                                      Self;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

  typedef vertex_descriptor                                     key_type;
  typedef EPICK::Point_3                                        value_type;
  typedef value_type&                                           reference;
  typedef boost::lvalue_property_map_tag                        category;

  Custom_VPM(std::map<key_type, value_type>& points) : points(points) { }

  friend void put(const Self& m, const key_type& k, const value_type& v) { m.points[k] = value_type(v.x(), v.y(), v.z()); }
  friend reference get(const Self& m, const key_type& k) { return m.points[k]; }

  std::map<key_type, value_type>& points;
};

template<class Mesh>
void test_bgl_STL(const std::string filename)
{
  Mesh fg;

  bool ok = CGAL::IO::read_STL(filename, fg);
  assert(ok);
  ok = CGAL::IO::write_STL("tmp.stl", fg);
  assert(ok);

  // write with ASCII in binary mode
  {
    ok = CGAL::IO::write_polygon_mesh("ascii.stl", fg, CGAL::parameters::use_binary_mode(false));
    assert(ok);
    std::ifstream test_ascii("ascii.stl", std::ios::binary);
    Mesh fg2;
    ok = CGAL::IO::read_STL(test_ascii, fg2, CGAL::parameters::use_binary_mode(false));
    test_ascii.close();
    assert(ok);
    assert(num_vertices(fg) == num_vertices(fg2) && num_faces(fg) == num_faces(fg2));
  }

  CGAL::clear(fg);

  std::map<typename boost::graph_traits<Mesh>::vertex_descriptor, EPICK::Point_3> cpoints;
  Custom_VPM<Mesh> cvpm(cpoints);

  std::ifstream is(filename, std::ios::binary);
  CGAL::IO::set_mode(is, CGAL::IO::BINARY);
  ok = CGAL::IO::read_STL(is, fg, CGAL::parameters::vertex_point_map(cvpm));
  assert(ok);
  assert(filename != "data/pig.stl" || (num_vertices(fg) == 8642 && num_faces(fg) == 16848));
  assert(filename != "data/pig.stl" || cpoints.size() == 8642);

  // write with STL
  {
    std::ofstream os("tmp.stl");
    ok = CGAL::IO::write_STL(os, fg, CGAL::parameters::vertex_point_map(cvpm));
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_STL("tmp.stl", fg2, CGAL::parameters::vertex_point_map(cvpm));
    assert(ok);
    assert(num_vertices(fg) == num_vertices(fg2) && num_faces(fg) == num_faces(fg2));
  }

  // write with PM
  {
    ok = CGAL::IO::write_polygon_mesh("tmp.stl", fg, CGAL::parameters::vertex_point_map(cvpm));
    assert(ok);
    cpoints.clear();
    Mesh fg2;
    ok = CGAL::IO::read_polygon_mesh("tmp.stl", fg2, CGAL::parameters::vertex_point_map(cvpm));
    assert(ok);
    assert(num_vertices(fg) == num_vertices(fg2) && num_faces(fg) == num_faces(fg2));
  }

  std::cerr << " ########### Error text is expected to follow." << std::endl;
  ok = CGAL::IO::read_STL("data/mesh_that_doesnt_exist.stl", fg);
  assert(!ok);
  ok = CGAL::IO::read_STL("data/invalid_cut.stl", fg); // cut in half
  assert(!ok);
  ok = CGAL::IO::read_STL("data/invalid_header.stl", fg); // missing solid
  assert(!ok);
  ok = CGAL::IO::read_STL("data/sphere.obj", fg);
  assert(!ok);
  ok = CGAL::IO::read_STL("data/full.off", fg);
  assert(!ok);
  std::cerr << " ########### No more error text from here." << std::endl;
}

template<class Mesh>
void test_bgl_GOCAD(const char* filename)
{
  Mesh fg;
  std::ifstream is(filename);
  bool ok = CGAL::IO::read_GOCAD(is, fg);
  assert(ok);
  assert(num_vertices(fg) != 0 && num_faces(fg) != 0);

  is.seekg(0);
  fg.clear();
  CGAL::clear(fg);
  std::pair<std::string, std::string> name_and_color;
  ok = CGAL::IO::read_GOCAD(is, name_and_color, fg);
  assert(ok);
  assert(num_vertices(fg) != 0 && num_faces(fg) != 0);

  // write with GOCAD
  {
    std::ofstream os("tmp.ts");
    bool ok = CGAL::IO::write_GOCAD(os, "tetrahedron", fg);
    assert(ok);

    Mesh fg2;
    std::pair<std::string, std::string> cnn;
    ok = CGAL::IO::read_GOCAD("tmp.ts", cnn, fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
    assert(cnn.first == "tetrahedron");
  }

  // write with PM
  {
    ok = CGAL::IO::write_polygon_mesh("tmp.ts", fg, CGAL::parameters::stream_precision(10));
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_polygon_mesh("tmp.ts", fg2);
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // test NPs
  typedef typename boost::property_map<Mesh,CGAL::vertex_point_t>::type VertexPointMap;

  VertexPointMap vpm = get(CGAL::vertex_point, fg);

  std::ostringstream out;
  ok = CGAL::IO::write_GOCAD(out, "tetrahedron", fg, CGAL::parameters::vertex_point_map(vpm));
  assert(ok);

  {
    Mesh fg2;
    VertexPointMap vpm2 = get(CGAL::vertex_point, fg2);
    std::istringstream is(out.str());
    std::pair<std::string, std::string> cnn;
    ok = CGAL::IO::read_GOCAD(is, cnn, fg2, CGAL::parameters::vertex_point_map(vpm2));
    assert(ok);
    assert(cnn.second.empty());
    assert(num_vertices(fg2) == 12491);
    assert(num_faces(fg2) == 24191);
  }

  std::cerr << " ########### Error text is expected to follow." << std::endl;
  ok = CGAL::IO::read_GOCAD("data/mesh_that_doesnt_exist.ts", fg);
  assert(!ok);
  ok = CGAL::IO::read_GOCAD("data/invalid_cut.ts", fg); // cut in half
  assert(!ok);
  ok = CGAL::IO::read_GOCAD("data/invalid_header.ts", fg); // missing header
  assert(!ok);
  ok = CGAL::IO::read_GOCAD("data/sphere.obj", fg);
  assert(!ok);
  ok = CGAL::IO::read_GOCAD("data/full.off", fg);
  assert(!ok);
  std::cerr << " ########### No more error text from here." << std::endl;
}

#ifdef CGAL_USE_VTK

template<typename Mesh, typename K>
void test_bgl_VTP(const char* filename,
                  const bool binary = false)
{
  Mesh fg;
  bool ok = CGAL::IO::read_VTP(filename, fg);
  assert(ok);
  assert(std::string(filename) != "data/bones.vtp" ||
         (num_vertices(fg) == 2154 && num_faces(fg) == 4204));

  // write with VTP
  {
    std::ofstream os;
    if(binary){
      os.open("tmp.vtp", std::ios::binary);
      CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    }
    else
      os.open("tmp.vtp");
    ok = CGAL::IO::write_VTP(os, fg, CGAL::parameters::use_binary_mode(binary));
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_VTP("tmp.vtp", fg2, CGAL::parameters::use_binary_mode(binary));
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // write with PM
  {
    if(binary)
      ok = CGAL::IO::write_polygon_mesh("tmp.vtp", fg);
    else
      ok = CGAL::IO::write_polygon_mesh("tmp.vtp", fg, CGAL::parameters::use_binary_mode(false));
    assert(ok);

    Mesh fg2;
    ok = CGAL::IO::read_polygon_mesh("tmp.vtp", fg2, CGAL::parameters::use_binary_mode(binary));
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  typedef typename K::Point_3 Point;
  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<Point> >::type VertexPointMap;
  // Test NPs
  {
    CGAL::clear(fg);

    ok = CGAL::IO::read_VTP("data/bones.vtp", fg);
    assert(ok);
    assert(num_vertices(fg) == 2154 && num_faces(fg) == 4204);

    Mesh fg2;

    VertexPointMap vpm2 = get(CGAL::dynamic_vertex_property_t<Point>(), fg2);
    ok = CGAL::IO::read_VTP("data/bones.vtp", fg2, CGAL::parameters::vertex_point_map(vpm2));
    assert(ok);
    typedef typename CGAL::GetInitializedVertexIndexMap<Mesh>::const_type      VIM;
    VIM vim1 = CGAL::get_initialized_vertex_index_map(fg);
    VIM vim2 = CGAL::get_initialized_vertex_index_map(fg2);
    assert(are_equal_meshes(fg, get(CGAL::vertex_point, fg), fg2, vpm2, vim1, vim2));
  }

  // write with VTP
  {
    std::ofstream os;
    if(binary)
    {
      os.open("tmp.vtp", std::ios::binary);
      CGAL::IO::set_mode(os, CGAL::IO::BINARY);
    }
    else
    {
      os.open("tmp.vtp");
    }
    ok = CGAL::IO::write_VTP(os, fg, CGAL::parameters::use_binary_mode(binary));
    assert(ok);

    Mesh fg2;

    ok = CGAL::IO::read_polygon_mesh("tmp.vtp", fg2, CGAL::parameters::use_binary_mode(binary));
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // write with PM
  {
    if(binary)
      ok = CGAL::IO::write_polygon_mesh("tmp.vtp", fg);
    else
      ok = CGAL::IO::write_polygon_mesh("tmp.vtp", fg, CGAL::parameters::use_binary_mode(false));
    assert(ok);

    Mesh fg2;

    ok = CGAL::IO::read_polygon_mesh("tmp.vtp", fg2, CGAL::parameters::use_binary_mode(binary));
    assert(ok);
    assert(are_equal_meshes(fg, fg2));
  }

  // test wrong inputs
  std::cerr << " ########### Error text is expected to follow." << std::endl;
  ok = CGAL::IO::read_VTP("data/mesh_that_doesnt_exist.vtp", fg);
  assert(!ok);
  ok = CGAL::IO::read_VTP("data/invalid_cut.vtp", fg); // cut in half
  assert(!ok);
  ok = CGAL::IO::read_VTP("data/invalid_header.vtp", fg); // missing header
  assert(!ok);
  ok = CGAL::IO::read_VTP("data/wrong_nb_points.vtp", fg); // wrong number of points
  assert(!ok);
  ok = CGAL::IO::read_VTP("data/sphere.obj", fg);
  assert(!ok);
  ok = CGAL::IO::read_VTP("data/full.off", fg);
  assert(!ok);
  ok = CGAL::IO::read_VTP("corrupted_bin.vtp", fg);
  assert(!ok);
  std::cerr << " ########### No more error text from here." << std::endl;
}

#endif // CGAL_USE_VTK

int main(int argc, char** argv)
{
  // OFF

  const char* off_file = (argc > 1) ? argv[1] : "data/prim.off";

  test_bgl_OFF<Polyhedron, Kernel>(off_file);
  Polyhedron fg;
  bool ok = CGAL::IO::read_OFF("data/invalid_header.off", fg); // wrong header (NOFF but no normals)
  assert(ok);

  test_bgl_OFF<SM, Kernel>(off_file);
  SM sm;
  ok = CGAL::IO::read_OFF("data/invalid_header.off", sm); // wrong header (NOFF but no normals)
  assert(!ok);
  test_bgl_OFF<LCC, Kernel>(off_file);
  LCC lcc;
  ok = CGAL::IO::read_OFF("data/invalid_header.off", lcc); // wrong header (NOFF but no normals)
  assert(!ok);
#ifdef CGAL_USE_OPENMESH
  test_bgl_OFF<OMesh, EPICK>(off_file);
  OMesh om;
  ok = CGAL::IO::read_OFF("data/invalid_header.off", om); // wrong header (NOFF but no normals)
  assert(!ok);
#endif
  // OBJ
  const char* obj_file = (argc > 2) ? argv[2] : "data/sphere.obj";
  test_bgl_OBJ<Polyhedron, Kernel>(obj_file);
  test_bgl_OBJ<SM, Kernel>(obj_file);
  test_bgl_OBJ<LCC, Kernel>(obj_file);
#ifdef CGAL_USE_OPENMESH
  test_bgl_OBJ<OMesh, EPICK>(obj_file);
#endif

  // PLY
  const char* ply_file_ascii = (argc > 3) ? argv[3] : "data/colored_tetra.ply";
  test_bgl_PLY<Polyhedron>(ply_file_ascii, false);
  test_bgl_PLY<SM>(ply_file_ascii, false);

  const char* ply_file = (argc > 3) ? argv[3] : "data/colored_tetra.ply";
  test_bgl_PLY<Polyhedron>(ply_file, true);
  test_bgl_PLY<SM>(ply_file, true);

  // STL
  const char* stl_file = (argc > 4) ? argv[4] : "data/pig.stl";
  test_bgl_STL<Polyhedron>(stl_file);
  test_bgl_STL<SM>(stl_file);
  test_bgl_STL<LCC>(stl_file);
#ifdef CGAL_USE_OPENMESH
  test_bgl_STL<OMesh>(stl_file);
#endif

  // GOCAD
  const char* gocad_file = (argc > 5) ? argv[5] : "data/2016206_MHT_surface.ts";
  test_bgl_GOCAD<Polyhedron>(gocad_file);
  test_bgl_GOCAD<SM>(gocad_file);
 test_bgl_GOCAD<LCC>(gocad_file);
#ifdef CGAL_USE_OPENMESH
  test_bgl_GOCAD<OMesh>(gocad_file);
#endif

  // VTP
#ifdef CGAL_USE_VTK
  const char* vtp_file = (argc > 6) ? argv[6] : "data/bones.vtp";

  test_bgl_VTP<Polyhedron, Kernel>(vtp_file, false);
  test_bgl_VTP<SM, Kernel>(vtp_file, false);
  test_bgl_VTP<LCC, Kernel>(vtp_file, false);

  test_bgl_VTP<Polyhedron, Kernel>(vtp_file, true);
  test_bgl_VTP<SM, Kernel>(vtp_file, true);
  test_bgl_VTP<LCC, Kernel>(vtp_file, true);

#ifdef CGAL_USE_OPENMESH
  test_bgl_VTP<OMesh, EPICK>(vtp_file, false);
  test_bgl_VTP<OMesh, EPICK>(vtp_file, true);
#endif

#endif // CGAL_USE_VTK

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
