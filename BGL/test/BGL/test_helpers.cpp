#include <iostream>
#include <fstream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/helpers.h>

typedef CGAL::Simple_cartesian<double> K;

typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

void
test(char *fname, bool triangle, bool quad, bool tetrahedron, bool hexahedron)
{
  std::cerr << "test(" << fname << ")"<< std::endl;
  Mesh m;
  std::ifstream in(fname);
  in >> m;
  halfedge_descriptor hd = *halfedges(m).first;
  assert(CGAL::is_triangle(hd, m) == triangle);
  assert(CGAL::is_quad(hd, m) == quad);
  assert(CGAL::is_tetrahedron(hd, m) == tetrahedron);
  assert(CGAL::is_hexahedron(hd, m) == hexahedron);
  }

int main()
{
  //                            triangle  quad    tetra   hexa
  test("data/triangle.off",     true,     false,  false,  false );
  test("data/quad.off",         false,    true,   false,  false );
  test("data/tetrahedron.off",  false,    false,  true,   false );
  test("data/cube.off",         false,    false,  false,  false );
  test("data/cube-quads.off",   false,    false,  false,  true );

  std::cerr << "done" << std::endl; 
  return 0;
}
