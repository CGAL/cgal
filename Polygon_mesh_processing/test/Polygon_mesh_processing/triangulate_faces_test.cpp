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
bool
test_triangulate_faces()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/cube_quad.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }
  
  bool success =
  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
  assert(CGAL::is_triangle_mesh(mesh));

  return success;
}

template <typename K>
bool
test_triangulate_face_range()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/cube_quad.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  bool success =
    CGAL::Polygon_mesh_processing::triangulate_faces(faces(mesh), mesh);
  assert(CGAL::is_triangle_mesh(mesh));

  return success;
}

template <typename K>
bool
test_triangulate_face()
{
  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input("data/cube_quad.off");

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  unsigned int nb = 0;
  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
  {
    if (nb > 4)
      break;
    else if (next(next(halfedge(fit, mesh), mesh), mesh)
             !=   prev(halfedge(fit, mesh), mesh))
    {
      if(CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh))
        ++nb;
      else assert(false);
    }
  }
  return true;
}

int main()
{
  assert(test_triangulate_faces<Epic>());
  assert(test_triangulate_face_range<Epic>());
  assert(test_triangulate_face<Epic>());

  assert(test_triangulate_faces<Epec>());
  assert(test_triangulate_face_range<Epec>());
  assert(test_triangulate_face<Epec>());

  return EXIT_SUCCESS;
}
