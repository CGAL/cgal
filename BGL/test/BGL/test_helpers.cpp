#include <iostream>
#include <fstream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Euler_operations.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;


template <typename Mesh>
void
test(const char *fname, bool triangle, bool quad, bool tetrahedron, bool hexahedron)
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

int main()
{
  typedef CGAL::Surface_mesh<Point_3> Mesh;
  //                                  triangle  quad    tetra   hexa
  test<Mesh>("data/triangle.off",     true,     false,  false,  false );
  test<Mesh>("data/quad.off",         false,    true,   false,  false );
  test<Mesh>("data/tetrahedron.off",  false,    false,  true,   false );
  test<Mesh>("data/cube.off",         false,    false,  false,  false );
  test<Mesh>("data/cube-quads.off",   false,    false,  false,  true );

  typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

  Point_3 a(0,0,0), b(1,0,0), c(1,1,0), d(0,1,0);
  Point_3 aa(0,0,1), bb(1,0,1), cc(1,1,1), dd(0,1,1);
  Mesh m;
  halfedge_descriptor hd;
  hd = CGAL::make_triangle(a,b,c,m);
  assert(CGAL::is_isolated_triangle(hd,m));
  assert(CGAL::is_valid(m));
  m.clear();
  hd = CGAL::make_quad(a,b,c,d,m);
  assert(CGAL::is_isolated_quad(hd,m));
  assert(CGAL::is_valid(m));
  assert(CGAL::is_quad_mesh(m));
  m.clear();
  hd = CGAL::make_tetrahedron(a,b,c,d,m);
  assert(CGAL::is_tetrahedron(hd,m));
  assert(CGAL::is_triangle_mesh(m));
  assert(CGAL::is_valid(m));
  m.clear();
  hd = CGAL::make_hexahedron(a,b,c,d,aa,bb,cc,dd,m);
  assert(CGAL::is_hexahedron(hd,m));
  assert(CGAL::is_quad_mesh(m));
  assert(CGAL::is_valid(m));

  std::cerr << "done" << std::endl; 
  return 0;
}
