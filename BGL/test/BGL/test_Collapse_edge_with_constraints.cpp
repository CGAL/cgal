#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/generators.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef Mesh::Property_map<Mesh::Edge_index, bool> ECM;


bool test_one_side(Mesh::Halfedge_index h, Mesh m)
{
  Mesh::Vertex_index vkept=target(h, m);
  CGAL::Euler::collapse_edge(edge(h, m), m, m.property_map<Mesh::Edge_index, bool>("ecm").value());
  return (!m.is_removed(vkept) && CGAL::is_valid_polygon_mesh(m));
}

bool test(Mesh::Halfedge_index h, Mesh& m)
{
  return test_one_side(h, m) && test_one_side(opposite(h, m), m);
}

int main()
{
  // ---------------------------------------------- //
  // two faces incident to the edge to be collapsed //
  // ---------------------------------------------- //
  Point_3 p1(0,0,0), p2(1,0,0), p3(0,1,0);
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(next(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
// duplicate block + add one border (1)
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(next(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  // duplicate block + add one border (2)
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(next(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
// duplicate block + add one border (1)
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c,m), m), true);
  put(ecm, edge(prev(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index h2c = CGAL::Euler::add_center_vertex(h, m);
  Mesh::Halfedge_index hb1 = opposite(prev(h2c, m), m);
  assert(is_border(hb1, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2 = opposite(next(opposite(h2c,m), m), m);
  assert(is_border(hb2, m));
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c,m), m), true);
  put(ecm, edge(next(opposite(h2c,m),m), m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  // ---------------------------------------------- //
  // one face incident to the edge to be collapsed  //
  // ---------------------------------------------- //
  // center triangle with 2 border edges (1)
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb1=opposite(next(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb1=opposite(next(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c, m),m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb1=opposite(next(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c, m),m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  // center triangle with 2 border edges (2)
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb2=opposite(prev(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb2=opposite(prev(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c, m),m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb2=opposite(prev(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c, m),m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  // center triangle with 1 border edges
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb1=opposite(next(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2=opposite(prev(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  /* ECM ecm =  */m.add_property_map<Mesh::Edge_index, bool>("ecm", false);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb1=opposite(next(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2=opposite(prev(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(next(h2c, m),m), true);
  bool res = test(h2c, m);
  assert(res);
  }
  {
  Mesh m;
  Mesh::Halfedge_index h2c = CGAL::make_triangle(p1,p2,p3,m);
  Mesh::Halfedge_index hb1=opposite(next(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb1, m), hb1, m);
  Mesh::Halfedge_index hb2=opposite(prev(h2c,m),m);
  CGAL::Euler::add_vertex_and_face_to_border(prev(hb2, m), hb2, m);
  ECM ecm = m.add_property_map<Mesh::Edge_index, bool>("ecm", false).first;
  put(ecm, edge(prev(h2c, m),m), true);
  bool res = test(h2c, m);
  assert(res);
  }
}
