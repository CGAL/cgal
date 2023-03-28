#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#include <boost/graph/filtered_graph.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

template <typename K>
bool
test_triangulate_faces()
{
  std::cout << "\n--- test_triangulate_faces(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/cube_quad.off"));

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
  assert(CGAL::is_triangle_mesh(mesh));

  return success;
}

template <typename K>
bool
test_triangulate_faces_with_named_parameters()
{
  std::cout << "\n--- test_triangulate_faces_with_named_parameters(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Epic::Point_3>      Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/cube_quad.off"));

  typedef Surface_mesh::Property_map<Surface_mesh::Vertex_index, Point>  VertexPointMap;
  VertexPointMap custom_vpm = mesh.add_property_map<Surface_mesh::Vertex_index, Point>("exact_vpm", Point()).first;

  typedef typename CGAL::property_map_selector<Surface_mesh, boost::vertex_point_t>::const_type CVPM;
  CVPM cvpm = CGAL::get_const_property_map(CGAL::vertex_point, mesh);

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  CGAL::Cartesian_converter<Epic, K> to_custom;

  for(boost::graph_traits<Surface_mesh>::vertex_descriptor vd : vertices(mesh))
  {
    put(custom_vpm, vd, to_custom(get(cvpm, vd)));
  }

  bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh,
                                                                  CGAL::parameters::vertex_point_map(custom_vpm)
                                                                                   .geom_traits(K()));
  assert(CGAL::is_triangle_mesh(mesh));

  return success;
}

template <typename K>
bool
test_triangulate_face_range(const std::string& filename)
{
  std::cout << "\n--- test_triangulate_face_range(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(filename);

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  bool success = CGAL::Polygon_mesh_processing::triangulate_faces(faces(mesh), mesh);

  for(auto f : faces(mesh))
  {
    if(!is_triangle(halfedge(f, mesh), mesh))
    {
      std::cout << "non triangular face:" << std::endl;
      for(auto h : halfedges_around_face(halfedge(f, mesh), mesh))
        std::cout << "  " << mesh.point(target(h, mesh)) << std::endl;
      assert(false);
    }
  }

  assert(CGAL::is_triangle_mesh(mesh));

  // For compilation
  CGAL::Polygon_mesh_processing::triangulate_faces(faces(mesh), mesh,
                                                   CGAL::parameters::geom_traits(K()));

  return success;
}

template <typename K>
bool
test_triangulate_face()
{
  std::cout << "\n--- test_triangulate_face(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/cube_quad.off"));

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  for(typename boost::graph_traits<Surface_mesh>::face_descriptor fit : faces(mesh))
  {
    if (next(next(halfedge(fit, mesh), mesh), mesh) !=   prev(halfedge(fit, mesh), mesh))
    {
      if(!CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh))
        assert(false);
    }
  }

  assert(CGAL::is_triangle_mesh(mesh));

  return true;
}

template <typename K>
bool
test_triangulate_triangle_face()
{
  std::cout << "\n--- test_triangulate_triangle_face(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/reference_tetrahedron.off"));

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  for(typename boost::graph_traits<Surface_mesh>::face_descriptor fit : faces(mesh))
  {
    if(!CGAL::Polygon_mesh_processing::triangulate_face(fit, mesh, CGAL::parameters::geom_traits(K())))
      assert(false);
  }

  assert(CGAL::is_triangle_mesh(mesh));

  return true;
}

// todo: add this in Dual.h
template <class SurfaceMesh, class Point, class Primal_map>
struct Dual_vpm
{
  typedef typename boost::graph_traits<SurfaceMesh>::face_descriptor key_type;
  typedef Point reference;
  typedef Point value_type;
  typedef boost::readable_property_map_tag category;

  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;

  Dual_vpm(const SurfaceMesh& primal,
           const Primal_map& primal_map)
    : primal_(primal)
    , primal_map_(primal_map) {}

  const Primal_map& primal_map() const
  {
    return primal_map_;
  }

  const SurfaceMesh& primal() const
  {
    return primal_;
  }

  friend
  value_type get(Dual_vpm& map, key_type f)
  {
    std::vector<Point> face_points;

    for(vertex_descriptor v :
                  CGAL::vertices_around_face(halfedge(f, map.primal()), map.primal()))
    {
      face_points.push_back( get(map.primal_map(), v) );
    }

    // temp extra copy
    Point centroid = CGAL::centroid(face_points.begin(), face_points.end(),
                                    CGAL::Dimension_tag<0>());

    return centroid;
  }

  const SurfaceMesh& primal_;
  Primal_map primal_map_;
};

template <typename K>
bool
test_dual_with_various_faces()
{
  std::cout << "\n--- test_dual_with_various_faces(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/elephant.off"));

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  typedef typename boost::property_map<Surface_mesh, boost::vertex_point_t>::type Pmap;
  Pmap vpmap = get_property_map(boost::vertex_point, mesh);

  CGAL::Dual<Surface_mesh> dual(mesh);
  // copy dual to a sm
  Surface_mesh sm_dual;
  CGAL::copy_face_graph(dual, sm_dual,
                        CGAL::parameters::vertex_point_map(
                          Dual_vpm<Surface_mesh, Point, Pmap>(mesh, vpmap)));

  for(typename boost::graph_traits<Surface_mesh>::face_descriptor fit : faces(sm_dual))
  {
    if(!CGAL::Polygon_mesh_processing::triangulate_face(fit, sm_dual,
                                                        CGAL::parameters::use_2d_constrained_delaunay_triangulation(true)))
      assert(false);
  }

  assert(CGAL::is_triangle_mesh(sm_dual));

  return true;
}

template <typename K>
bool
test_triangulate_soup()
{
  std::cout << "\n--- test_triangulate_soup(" << typeid(K).name() << ") ---" << std::endl;

  typedef typename K::Point_3                    Point;
  typedef CGAL::Surface_mesh<Point>          Surface_mesh;

  Surface_mesh mesh;
  std::ifstream input(CGAL::data_file_path("meshes/elephant.off"));

  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return false;
  }

  typedef typename boost::property_map<Surface_mesh, boost::vertex_point_t>::type Pmap;
  Pmap vpmap = get_property_map(boost::vertex_point, mesh);

  CGAL::Dual<Surface_mesh> dual(mesh);
  // copy dual to a sm
  Surface_mesh sm_dual;
  CGAL::copy_face_graph(dual, sm_dual,
                        CGAL::parameters::vertex_point_map(
                          Dual_vpm<Surface_mesh, Point, Pmap>(mesh, vpmap)));

  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > polygons;
  CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(sm_dual, points, polygons);

  bool success = CGAL::Polygon_mesh_processing::triangulate_polygons(points, polygons,
                                                                     CGAL::parameters::geom_traits(K())
                                                                                      .use_2d_constrained_delaunay_triangulation(false));
  for(std::size_t i = 0; i < polygons.size(); ++i)
  {
    assert(polygons[i].size() == 3);
  }

  // For compilation
  success = CGAL::Polygon_mesh_processing::triangulate_polygons(points, polygons);

  return success;
}

int main(int argc, char** argv)
{
  if(argc > 1)
  {
    assert(test_triangulate_face_range<Epic>(argv[1]));
  }

  assert(test_triangulate_faces<Epic>());
  assert(test_triangulate_faces_with_named_parameters<Epic>());
  assert(test_triangulate_face_range<Epic>(CGAL::data_file_path("meshes/cube_quad.off")));
  assert(test_triangulate_face<Epic>());
  assert(test_triangulate_triangle_face<Epic>());
  assert(test_dual_with_various_faces<Epic>());
  assert(test_triangulate_soup<Epic>());

  assert(test_triangulate_faces<Epec>());
  assert(test_triangulate_faces_with_named_parameters<Epec>());
  assert(test_triangulate_face_range<Epec>(CGAL::data_file_path("meshes/cube_quad.off")));
  assert(test_triangulate_face<Epec>());
  assert(test_triangulate_triangle_face<Epec>());
  assert(test_dual_with_various_faces<Epec>());
  assert(test_triangulate_soup<Epec>());

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
