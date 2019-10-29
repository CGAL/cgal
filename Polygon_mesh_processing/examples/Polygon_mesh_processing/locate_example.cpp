#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::FT                                                           FT;
typedef K::Point_2                                                      Point_2;
typedef K::Ray_2                                                        Ray_2;
typedef K::Point_3                                                      Point_3;
typedef K::Ray_3                                                        Ray_3;

typedef CGAL::Surface_mesh<Point_3>                                     Mesh;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor;
typedef typename boost::graph_traits<Mesh>::face_descriptor             face_descriptor;

namespace CP = CGAL::parameters;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
typedef PMP::Face_location<Mesh, FT>                                    Face_location;

int main(int /*argc*/, char** /*argv*/)
{
  // Generate a simple 3D triangle mesh that with vertices on the plane xOy
  Mesh tm;
  CGAL::make_grid(10, 10, tm);
  PMP::triangulate_faces(tm);

  // Basic usage
  Face_location random_location = PMP::random_location_on_mesh<FT>(tm);
  const face_descriptor f = random_location.first;
  const Barycentric_coordinates& coordinates = random_location.second;

  std::cout << "Random location on the mesh: face " << f
            << " and with coordinates [" << coordinates[0] << "; "
                                         << coordinates[1] << "; "
                                         << coordinates[2] << "]\n";
  std::cout << "It corresponds to point (" << PMP::construct_point(random_location, tm) << ")\n\n";

  // Locate a known 3D point in the mesh
  const Point_3 query(1.2, 7.4, 0);
  Face_location query_location = PMP::locate(query, tm);

  std::cout << "Point (" << query << ") is located in face " << query_location.first
            << " with barycentric coordinates [" << query_location.second[0] << "; "
                                                 << query_location.second[1] << "; "
                                                 << query_location.second[2] << "]\n\n";

  // Locate a 3D point in the mesh as the intersection of the mesh and a 3D ray.
  // The AABB tree can be cached in case many queries are performed (otherwise, it is rebuilt
  // on each call, which is expensive).
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh>                AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>               AABB_face_graph_traits;

  CGAL::AABB_tree<AABB_face_graph_traits> tree;
  PMP::build_AABB_tree(tm, tree);

  const Ray_3 ray_3(Point_3(4.2, 6.8, 2.4), Point_3(7.2, 2.3, -5.8));
  Face_location ray_location = PMP::locate_with_AABB_tree(ray_3, tree, tm);

  std::cout << "Intersection of the 3D ray and the mesh is in face " << ray_location.first
            << " with barycentric coordinates [" << ray_location.second[0] << " "
                                                 << ray_location.second[1] << " "
                                                 << ray_location.second[2] << "]\n";
  std::cout << "It corresponds to point (" << PMP::construct_point(ray_location, tm) << ")\n";
  std::cout << "Is it on the face's border? " << (PMP::is_on_face_border(ray_location, tm) ? "Yes" : "No") << "\n\n";

  // -----------------------------------------------------------------------------------------------
  // Now, we artifically project the mesh to the natural 2D dimensional plane, with a little translation
  // via a custom vertex point property map

  typedef CGAL::dynamic_vertex_property_t<Point_2>                      Point_2_property;
  typedef typename boost::property_map<Mesh, Point_2_property>::type    Projection_pmap;
  Projection_pmap projection_pmap = get(Point_2_property(), tm);

  for(vertex_descriptor v : vertices(tm))
  {
    const Point_3& p = tm.point(v);
    put(projection_pmap, v, Point_2(p.x() + 1, p.y())); // simply ignoring the z==0 coordinate and translating along Ox
  }

  // Locate the same 3D point but in a 2D context
  const Point_2 query_2(query.x() + 1, query.y());
  Face_location query_location_2 = PMP::locate(query_2, tm, CP::vertex_point_map(projection_pmap));

  std::cout << "Point (" << query_2 << ") is located in face " << query_location_2.first
            << " with barycentric coordinates [" << query_location_2.second[0] << "; "
                                                 << query_location_2.second[1] << "; "
                                                 << query_location_2.second[2] << "]\n\n";

  // Shoot a 2D ray and locate the intersection with the mesh in 2D
  const Ray_2 ray_2(Point_2(-10, -10), Point_2(10, 10));
  Face_location ray_location_2 = PMP::locate(ray_2, tm, CP::vertex_point_map(projection_pmap)); // This rebuilds an AABB tree on each call

  std::cout << "Intersection of the 2D ray and the mesh is in face " << ray_location_2.first
            << " with barycentric coordinates [" << ray_location_2.second[0] << "; "
                                                 << ray_location_2.second[1] << "; "
                                                 << ray_location_2.second[2] << "]\n";
  std::cout << "It corresponds to point (" << PMP::construct_point(ray_location_2, tm, CP::vertex_point_map(projection_pmap)) << ")\n";

  if(PMP::is_on_mesh_border(ray_location_2, tm))
    std::cout << "It is on the border of the mesh!\n" << std::endl;

  return EXIT_SUCCESS;
}
