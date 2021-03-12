#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>
#include <iostream>

namespace CGAL{
namespace Polygon_mesh_processing{
namespace internal {
template<typename TriangleMesh, typename FaceRange, typename FacePairRange>
void diff(const TriangleMesh& m1, const TriangleMesh& m2, bool compute_common = false)
{
  std::map<Point_3, std::size_t> point_id_map;
  std::vector<std::size_t> m1_vertex_id(num_vertices(m1), -1);
  std::vector<std::size_t> m2_vertex_id(num_vertices(m2), -1);

  //iterate both meshes to set ids to all points, and set vertex/point_id maps.
  std::size_t id =0;
  for(auto v : vertices(m1))
  {
    Point_3 p = m1->point(v);
    auto res = point_id_map.insert(std::make_pair(p, id));
    if(res.second)
      id++;
    m1_vertex_id[(std::size_t)v]=res.first->second;
  }
  for(auto v : vertices(m2))
  {
    Point_3 p = m2->point(v);
    auto res = point_id_map.insert(std::make_pair(p, id));
    if(res.second)
      id++;
    m2_vertex_id[(std::size_t)v]=res.first->second;
  }

  //fill a set with the "faces point-ids" of m1 and then iterate faces of m2 to compare.
  std::set<std::vector<std::size_t> > m1_faces;
  for(auto f : m1->faces())
  {
    std::vector<std::size_t> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m1), m1))
    {
      ids.push_back(m1_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    m1_faces.insert(ids);
  }

  std::vector<face_index> common_faces;
  id = 0;

  for(auto f : faces(m2))
  {
    std::vector<std::size_t> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m2), m2))
    {
      ids.push_back(m2_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    if(!((m1_faces.find(ids) != m1_faces.end()) ^ compute_common))
    {
      common_faces.push_back(f);
    }
  }
}

}
template<typename TriangleMesh, typename FaceRange, typename FacePairRange>
void compare_meshes(const TriangleMesh& tm1, const TriangleMesh& tm2,
                    FacePairRange& common, FaceRange tm1_only, FaceRange tm2_only)
{

  SMesh* m1_over_m2 = diff(m1, m2);
  SMesh* m2_over_m1 = diff(m2, m1);
  SMesh* common = diff(m2, m1, true);

  Scene_surface_mesh_item* m1_over_m2_item = new Scene_surface_mesh_item(m1_over_m2);
  m1_over_m2_item->setColor(QColor(Qt::blue));
  m1_over_m2_item->setName(QString("%2 - %1").arg(m1_item->name()).arg(m2_item->name()));
  CGAL::Three::Three::scene()->addItem(m1_over_m2_item);
  Scene_surface_mesh_item* m2_over_m1_item = new Scene_surface_mesh_item(m2_over_m1);
  m2_over_m1_item->setColor(QColor(Qt::red));
  m2_over_m1_item->setName(QString("%1 - %2").arg(m1_item->name()).arg(m2_item->name()));
  CGAL::Three::Three::scene()->addItem(m2_over_m1_item);
  Scene_surface_mesh_item* common_item = new Scene_surface_mesh_item(common);
  common_item->setColor(QColor(Qt::green));
  CGAL::Three::Three::scene()->addItem(common_item);
  common_item->setName(QString("%1 && %2").arg(m1_item->name()).arg(m2_item->name()));
}

}
}
typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef K::Point_3                                                Point;

typedef CGAL::Surface_mesh<Point>                                 Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor        face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/eight.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight-bis.off";

  Surface_mesh mesh1, mesh2;
  if(!PMP::read_polygon_mesh(filename1, mesh1))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  if(!PMP::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  return 0;
}
