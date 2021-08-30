#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                                    K;
typedef K::Point_3                                                        Point_3;

template <typename Mesh>
void test_validity()
{
  std::cerr << "test validity" << std::endl;

  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;

  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMap;

  Mesh mesh;
  VPMap vpmap = get(CGAL::vertex_point, mesh);

  vertex_descriptor vertices[4];
  edge_descriptor edges[4];
  vertices[0] = add_vertex(mesh);
  vertices[1] = add_vertex(mesh);
  vertices[2] = add_vertex(mesh);
  vertices[3] = add_vertex(mesh);

  put(vpmap, vertices[0], Point_3(0,0,0));
  put(vpmap, vertices[1], Point_3(1,0,0));
  put(vpmap, vertices[2], Point_3(1,1,0));
  put(vpmap, vertices[3], Point_3(0,1,0));

  edges[0] = add_edge(mesh);
  edges[1] = add_edge(mesh);
  edges[2] = add_edge(mesh);
  edges[3] = add_edge(mesh);

  assert(!CGAL::is_valid_halfedge_graph(mesh));
  for(int i=0; i<4; ++i)
  {
    set_target(halfedge(edges[i], mesh), vertices[i], mesh);
    set_halfedge(vertices[i], halfedge(edges[i], mesh), mesh);
  }

  for(int i=0; i<4; ++i)
    set_target(opposite(halfedge(edges[i], mesh), mesh), vertices[(i+1)%4], mesh);
  for(int i=0; i<4; ++i)
  {
    set_next(halfedge(edges[(i+1)%4], mesh), halfedge(edges[i], mesh), mesh);
    set_next(opposite(halfedge(edges[i], mesh), mesh),
        opposite(halfedge(edges[(i+1)%4], mesh), mesh), mesh);
  }

  assert(CGAL::is_valid_halfedge_graph(mesh));
  face_descriptor faces[1];
  faces[0] = add_face(mesh);
  assert(!CGAL::is_valid_face_graph(mesh));

  for(int i=0; i<4; ++i)
  {
    set_face(opposite(halfedge(edges[i], mesh), mesh), faces[0], mesh);
  }
  set_halfedge(faces[0], opposite(halfedge(edges[0], mesh), mesh), mesh);
  assert(CGAL::is_valid_face_graph(mesh));
  assert(CGAL::is_valid_polygon_mesh(mesh));

  Mesh dummy;
  vertices[0] = add_vertex(dummy);
  vertices[1] = add_vertex(dummy);
  edges[0] = add_edge(dummy);
  set_target(halfedge(edges[0], dummy), vertices[0], dummy);
  set_halfedge(vertices[0], halfedge(edges[0], dummy), dummy);
  set_target(opposite(halfedge(edges[0], dummy), dummy), vertices[1], dummy);
  set_halfedge(vertices[1], opposite(halfedge(edges[0], dummy), dummy), dummy);
  set_next(halfedge(edges[0], dummy), opposite(halfedge(edges[0], dummy), dummy), dummy);
  set_next(opposite(halfedge(edges[0], dummy), dummy), halfedge(edges[0], dummy), dummy);
  faces[0] = add_face(dummy);
  set_halfedge(faces[0], opposite(halfedge(edges[0], dummy), dummy), dummy);
  set_face(halfedge(edges[0], dummy), faces[0], dummy);
  set_face(opposite(halfedge(edges[0], dummy), dummy), faces[0], dummy);
  assert(CGAL::is_valid_face_graph(dummy));
  assert(!CGAL::is_valid_polygon_mesh(dummy));

}

template <typename Mesh>
void test(const char *fname, bool triangle, bool quad, bool tetrahedron, bool hexahedron)
{
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

  std::cerr << "test(" << fname << ")"<< std::endl;

  Mesh m;
  std::ifstream in(fname);
  in >> m;

  halfedge_descriptor hd = *halfedges(m).first;
  assert(CGAL::is_isolated_triangle(hd, m) == triangle);
  assert(CGAL::is_isolated_quad(hd, m) == quad);
  assert(CGAL::is_tetrahedron(hd, m) == tetrahedron);
  assert(CGAL::is_hexahedron(hd, m) == hexahedron);
}

template <typename Mesh>
void test_generators()
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor      face_descriptor;

  //                                  triangle  quad    tetra   hexa
  test<Mesh>("data/triangle.off",     true,     false,  false,  false );
  test<Mesh>("data/quad.off",         false,    true,   false,  false );
  test<Mesh>("data/tetrahedron.off",  false,    false,  true,   false );
  test<Mesh>("data/cube.off",         false,    false,  false,  false );
  test<Mesh>("data/cube-quads.off",   false,    false,  false,  true );

  Point_3 a(0,0,0), b(1,0,0), c(1,1,0), d(0,1,0);
  Point_3 aa(0,0,1), bb(1,0,1), cc(1,1,1), dd(0,1,1);

  Mesh m;
  halfedge_descriptor hd;
  hd = CGAL::make_triangle(a,b,c,m);
  assert(CGAL::is_isolated_triangle(hd,m));
  assert(CGAL::is_valid_polygon_mesh(m));

  m.clear();
  hd = CGAL::make_quad(a,b,c,d,m);
  assert(CGAL::is_isolated_quad(hd,m));
  assert(CGAL::is_valid_polygon_mesh(m));
  assert(CGAL::is_quad_mesh(m));

  m.clear();
  hd = CGAL::make_tetrahedron(a,b,c,d,m);
  assert(CGAL::is_tetrahedron(hd,m));
  assert(CGAL::is_triangle_mesh(m));
  assert(CGAL::is_valid_polygon_mesh(m));
  m.clear();
  hd = CGAL::make_hexahedron(a,b,c,d,dd,aa,bb,cc,m);
  assert(CGAL::is_hexahedron(hd,m));
  assert(CGAL::is_quad_mesh(m));
  assert(CGAL::is_valid_polygon_mesh(m));

  m.clear();
  CGAL::make_icosahedron<Mesh, Point_3>(m);
  assert(num_faces(m) == 20);
  assert(CGAL::is_triangle_mesh(m));
  assert(CGAL::is_valid_polygon_mesh(m));

  m.clear();
  hd = CGAL::make_pyramid<Mesh, Point_3>(3, m);
  assert(num_faces(m) == 6);
  assert(CGAL::is_triangle_mesh(m));
  assert(CGAL::is_valid_polygon_mesh(m));

  m.clear();
  hd = CGAL::make_regular_prism<Mesh, Point_3>(4, m);
  assert(num_faces(m) == 16);
  assert(CGAL::is_triangle_mesh(m));
  assert(CGAL::is_valid_polygon_mesh(m));

  m.clear();
  CGAL::make_grid(3,3,m);
  assert(num_faces(m) == 9);
  assert(CGAL::is_quad_mesh(m));
  assert(CGAL::is_valid_polygon_mesh(m));

  // -----------------------------------------------------------------------------------------------

  std::cerr << "test random element generators" << std::endl;

  vertex_descriptor v;
  halfedge_descriptor h;
  edge_descriptor e;
  face_descriptor f;

  CGAL::Random rnd;

  // ---------------------------------------------------------------------------
  v = CGAL::internal::random_vertex_in_mesh(m, rnd);
  assert(v != boost::graph_traits<Mesh>::null_vertex());

  h = CGAL::internal::random_halfedge_in_mesh(m, rnd);
  assert(h != boost::graph_traits<Mesh>::null_halfedge());

  e = CGAL::internal::random_edge_in_mesh(m, rnd);

  f = CGAL::internal::random_face_in_mesh(m, rnd);
  assert(f != boost::graph_traits<Mesh>::null_face());

  // ---------------------------------------------------------------------------
  h = CGAL::internal::random_halfedge_in_face(f, m, rnd);
  assert(h != boost::graph_traits<Mesh>::null_halfedge());
  assert(face(h, m) == f);

  v = CGAL::internal::random_vertex_in_face(f, m, rnd);
  assert(v != boost::graph_traits<Mesh>::null_vertex());

  // could use vertices_around_face, but it's the point is not to
  bool has_vertex = false;
  halfedge_descriptor done = h;
  do
  {
    if(target(h, m) == v)
    {
      has_vertex = true;
      break;
    }

    h = next(h, m);
  }
  while(h != done);
  assert(has_vertex);
}

int main()
{
  typedef CGAL::Surface_mesh<Point_3> Mesh;

  test_validity<Mesh>();
  test_generators<Mesh>();

  std::cerr << "done" << std::endl;
  return EXIT_SUCCESS;
}
