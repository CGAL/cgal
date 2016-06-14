#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <boost/foreach.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

template <typename K>
int
test_triangulate_faces()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/elephant.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }
  
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh, 
                                                   CGAL::Polygon_mesh_processing::parameters::all_default());
  CGAL::Polygon_mesh_processing::triangulate_faces(faces(mesh), mesh);
  CGAL::Polygon_mesh_processing::triangulate_faces(faces(mesh), mesh, 
                                                   CGAL::Polygon_mesh_processing::parameters::all_default());


  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
    if (next(next(halfedge(fit, mesh), mesh), mesh)
        !=   prev(halfedge(fit, mesh), mesh)) {
      CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh);
      CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh,
                                                      CGAL::Polygon_mesh_processing::parameters::all_default());
    }
  return 0;
}

int main()
{
  int r = test_triangulate_faces<Epic>();
  assert(r==0);
  r = test_triangulate_faces<Epec>();
  assert(r==0);
  return 0;
}
