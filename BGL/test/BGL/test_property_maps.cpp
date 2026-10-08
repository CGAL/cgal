#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
// #include <CGAL/Seam_mesh.h>

#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>

#include <cassert>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>


using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point>;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;

// OpenMesh is optional in CGAL.
#ifdef CGAL_USE_OPENMESH
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

using OpenMeshT = OpenMesh::PolyMesh_ArrayKernelT< >;
#endif

// ---------------------------------------------------------------------------
// Vertex
// ---------------------------------------------------------------------------

template <typename Mesh>
void test_vertex_property(Mesh& mesh)
{
  using tag = CGAL::dynamic_vertex_property_t<int>;
  using vd = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::type;

  constexpr int default_value = 42;
  constexpr int assigned_value = 1234;

  pmap pm = get(tag(), mesh, default_value);
  for(vd v : vertices(mesh))
    assert(get(pm, v) == default_value);

  vd first = *vertices(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(vd v : vertices(mesh))
    if(v != first)
      assert(get(pm, v) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  int i=0;
  for(vd v: vertices(mesh))
    put(pm, v, ++i);
  i=0;
  for(vd v: vertices(mesh))
    assert(get(pm, v) == ++i);
}

template <typename Mesh>
void test_vertex_const_property(const Mesh& mesh)
{
  using tag = CGAL::dynamic_vertex_property_t<int>;
  using vd = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::const_type;

  constexpr int default_value = 42;
  constexpr int assigned_value = 1234;

  pmap pm = get(tag(), mesh, default_value);
  for(vd v : vertices(mesh))
    assert(get(pm, v) == default_value);

  vd first = *vertices(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(vd v : vertices(mesh))
    if(v != first)
      assert(get(pm, v) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  int i=0;
  for(vd v: vertices(mesh))
    put(pm, v, ++i);
  i=0;
  for(vd v: vertices(mesh))
    assert(get(pm, v) == ++i);
}




// ---------------------------------------------------------------------------
// Halfedge
// ---------------------------------------------------------------------------

template <typename Mesh>
void test_halfedge_property(Mesh& mesh)
{
  using tag = CGAL::dynamic_halfedge_property_t<std::string>;
  using hd = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::type;

  const std::string default_value = "default-halfedge";
  const std::string assigned_value = "assigned-halfedge";

  pmap pm = get(tag(), mesh, default_value);
  for(hd h : halfedges(mesh))
    assert(get(pm, h) == default_value);

  hd first = *halfedges(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(hd h : halfedges(mesh))
    if(h != first)
      assert(get(pm, h) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  int i = 0;
  for(hd h : halfedges(mesh))
    put(pm_without_default, h, std::to_string(++i));
  i = 0;
  for(hd h : halfedges(mesh))
    assert(get(pm_without_default, h) == std::to_string(++i));
}

template <typename Mesh>
void test_halfedge_const_property(const Mesh& mesh)
{
  using tag = CGAL::dynamic_halfedge_property_t<std::string>;
  using hd = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::const_type;

  const std::string default_value = "default-halfedge";
  const std::string assigned_value = "assigned-halfedge";

  pmap pm = get(tag(), mesh, default_value);
  for(hd h : halfedges(mesh))
    assert(get(pm, h) == default_value);

  hd first = *halfedges(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(hd h : halfedges(mesh))
    if(h != first)
      assert(get(pm, h) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  int i = 0;
  for(hd h : halfedges(mesh))
    put(pm_without_default, h, std::to_string(++i));
  i = 0;
  for(hd h : halfedges(mesh))
    assert(get(pm_without_default, h) == std::to_string(++i));
}


// ---------------------------------------------------------------------------
// Edge
// ---------------------------------------------------------------------------

template <typename Mesh>
void test_edge_property(Mesh& mesh)
{
  using tag = CGAL::dynamic_edge_property_t<double>;
  using ed = typename boost::graph_traits<Mesh>::edge_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::type;

  constexpr double default_value = 3.141592653589793;
  constexpr double assigned_value = 123.456;

  pmap pm = get(tag(), mesh, default_value);
  for(ed e : edges(mesh))
    assert(get(pm, e) == default_value);

  ed first = *edges(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(ed e : edges(mesh))
    if(e != first)
      assert(get(pm, e) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  int i = 0;
  for(ed e : edges(mesh))
    put(pm_without_default, e, static_cast<double>(++i));
  i = 0;
  for(ed e : edges(mesh))
    assert(get(pm_without_default, e) == static_cast<double>(++i));
}

template <typename Mesh>
void test_edge_const_property(const Mesh& mesh)
{
  using tag = CGAL::dynamic_edge_property_t<double>;
  using ed = typename boost::graph_traits<Mesh>::edge_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::const_type;

  constexpr double default_value = 3.141592653589793;
  constexpr double assigned_value = 123.456;

  pmap pm = get(tag(), mesh, default_value);
  for(ed e : edges(mesh))
    assert(get(pm, e) == default_value);

  ed first = *edges(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(ed e : edges(mesh))
    if(e != first)
      assert(get(pm, e) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  int i = 0;
  for(ed e : edges(mesh))
    put(pm_without_default, e, static_cast<double>(++i));
  i = 0;
  for(ed e : edges(mesh))
    assert(get(pm_without_default, e) == static_cast<double>(++i));
}


// ---------------------------------------------------------------------------
// Face
// ---------------------------------------------------------------------------

template <typename Mesh>
void test_face_property(Mesh& mesh)
{
  using tag = CGAL::dynamic_face_property_t<std::size_t>;
  using fd = typename boost::graph_traits<Mesh>::face_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::type;

  constexpr std::size_t default_value = 17;
  constexpr std::size_t assigned_value = 999;

  pmap pm = get(tag(), mesh, default_value);
  for(fd f : faces(mesh))
    assert(get(pm, f) == default_value);

  fd first = *faces(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(fd f : faces(mesh))
    if(f != first)
      assert(get(pm, f) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  std::size_t i = 0;
  for(fd f : faces(mesh))
    put(pm_without_default, f, ++i);
  i = 0;
  for(fd f : faces(mesh))
    assert(get(pm_without_default, f) == ++i);
}

template <typename Mesh>
void test_face_const_property(const Mesh& mesh)
{
  using tag = CGAL::dynamic_face_property_t<std::size_t>;
  using fd = typename boost::graph_traits<Mesh>::face_descriptor;
  using pmap = typename boost::property_map<Mesh, tag>::const_type;

  constexpr std::size_t default_value = 17;
  constexpr std::size_t assigned_value = 999;

  pmap pm = get(tag(), mesh, default_value);
  for(fd f : faces(mesh))
    assert(get(pm, f) == default_value);

  fd first = *faces(mesh).first;
  put(pm, first, assigned_value);
  assert(get(pm, first) == assigned_value);

  for(fd f : faces(mesh))
    if(f != first)
      assert(get(pm, f) == default_value);

  pmap pm_without_default = get(tag(), mesh);
  std::size_t i = 0;
  for(fd f : faces(mesh))
    put(pm_without_default, f, ++i);
  i = 0;
  for(fd f : faces(mesh))
    assert(get(pm_without_default, f) == ++i);
}


// ---------------------------------------------------------------------------
// Complete test for a mesh model
// ---------------------------------------------------------------------------

template <typename Mesh>
void test_mesh(const char* name)
{
  std::cout << "Testing " << name << "..." << std::endl;

  Mesh mesh;
  CGAL::make_tetrahedron(Point(0,0,0),
                         Point(0,0,1),
                         Point(0,1,0),
                         Point(1,0,0),
                         mesh);

  // assert(CGAL::num_vertices(mesh) == 4);
  // assert(CGAL::num_edges(mesh) == 6);
  // assert(CGAL::num_faces(mesh) == 4);

  test_vertex_property(mesh);
  test_halfedge_property(mesh);
  test_edge_property(mesh);
  test_face_property(mesh);

  test_vertex_const_property(mesh);
  test_halfedge_const_property(mesh);
  test_edge_const_property(mesh);
  test_face_const_property(mesh);

  std::cout << "  OK" << std::endl;
}

int main()
{
  test_mesh<Surface_mesh>("Surface_mesh");
  test_mesh<Polyhedron>("Polyhedron_3");
#ifdef CGAL_USE_OPENMESH
  test_mesh<OpenMeshT>("OpenMesh");
#endif
  return EXIT_SUCCESS;
}
