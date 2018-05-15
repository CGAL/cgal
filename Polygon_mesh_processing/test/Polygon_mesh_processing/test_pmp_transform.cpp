#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <iostream>
#include <fstream>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Aff_transformation_3.h>



namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;



int main()
{
  Mesh m;
  CGAL::make_tetrahedron(Point(0,0,0), Point(0,0,1), Point(0,1,2), Point(1,0,3), m);
  CGAL::Aff_transformation_3<K> trans(0,0,0,1,0,1,0,0,0,0,1,1);
  PMP::transform(trans, m, params::all_default());
  bool ok = true;
  for(Mesh::size_type i = 0; i<m.vertices().size(); ++i)
  {
    double z = m.point(Mesh::Vertex_index(i)).z();
    ok &= z == i+1;
  }
  CGAL_assertion(ok);
  return 0;
}
