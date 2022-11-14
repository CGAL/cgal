#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3_to_face_graph.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_3                                              Point_3;
typedef CGAL::Delaunay_triangulation_3<K>                       Delaunay;
typedef Delaunay::Vertex_handle                                 Vertex_handle;

typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> Surface_mesh;

int main()
{
  CGAL::Random_points_in_sphere_3<Point_3> gen(100.0);
  std::list<Point_3> points;

  // generate 250 points randomly in a sphere of radius 100.0
  // and insert them into the triangulation
  std::copy_n(gen, 250, std::back_inserter(points));
  Delaunay T;
  T.insert(points.begin(), points.end());

  std::list<Vertex_handle>  vertices;
  T.incident_vertices(T.infinite_vertex(), std::back_inserter(vertices));
  std::cout << "This convex hull of the 250 points has "
            << vertices.size() << " points on it." << std::endl;

  // remove 25 of the input points
  std::list<Vertex_handle>::iterator v_set_it = vertices.begin();
  for (int i = 0; i < 25; i++)
  {
     T.remove(*v_set_it);
     v_set_it++;
  }

  //copy the convex hull of points into a polyhedron and use it
  //to get the number of points on the convex hull
  Surface_mesh chull;
  CGAL::convex_hull_3_to_face_graph(T, chull);

  std::cout << "After removal of 25 points, there are "
            << num_vertices(chull) << " points on the convex hull." << std::endl;

  return 0;
}
