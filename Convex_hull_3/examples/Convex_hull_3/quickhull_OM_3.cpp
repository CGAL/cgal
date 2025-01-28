#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/io.h>

#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3
                              Point_3;
typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> Mesh;

int main(int argc, char* argv[])
{
  std::ifstream in((argc>1) ? argv[1] : CGAL::data_file_path("points_3/cube.xyz"));
  std::vector<Point_3> points;
  Point_3 p, n;
  while(in >> p >> n)
    points.push_back(p);

  // define polyhedron to hold convex hull
  Mesh poly;

  // compute convex hull of non-collinear points
  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  Mesh sm;
  CGAL::convex_hull_3(points.begin(), points.end(), sm);

  CGAL::IO::write_polygon_mesh((argc>2)?argv[2]:"out.off", sm);

  return 0;
}
