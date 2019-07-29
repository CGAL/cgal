#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>

#include "Scene_surface_mesh_item.h"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Three.h>
#include "Messages_interface.h"

using namespace CGAL::Three;
class Polyhedron_demo_diff_between_meshes_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface, Messages_interface*);

  bool applicable(QAction*) const {
    if(scene->selectionIndices().size() != 2)
      return false;
    for(int id : scene->selectionIndices())
    {
      if(!qobject_cast<Scene_surface_mesh_item*>(scene->item(id)))
      {
        return false;
      }
    }
    return true;
  }

  QList<QAction*> actions() const;

public Q_SLOTS:
  void diff();

private:
  CGAL::Three::Scene_interface* scene;
  QAction* actionDiff;
  SMesh* diff(SMesh* m1, SMesh* m2, bool compute_common);

}; // end Polyhedron_demo_diff_between_meshes_plugin

void Polyhedron_demo_diff_between_meshes_plugin::init(QMainWindow* mainWindow,
                                                          CGAL::Three::Scene_interface* scene_interface,
                                                          Messages_interface*)
{
  scene = scene_interface;
  actionDiff = new QAction(tr("&Differences between Meshes"), mainWindow);
  actionDiff->setObjectName("actionDiff");
  connect(actionDiff, SIGNAL(triggered()),
          this, SLOT(diff()));
}

QList<QAction*> Polyhedron_demo_diff_between_meshes_plugin::actions() const {
  return QList<QAction*>() << actionDiff;
}

SMesh* Polyhedron_demo_diff_between_meshes_plugin::diff(SMesh* m1, SMesh* m2, bool compute_common = false)
{
  std::map<Point_3, std::size_t> point_id_map;
  std::vector<std::size_t> m1_vertex_id(num_vertices(*m1), -1);
  std::vector<std::size_t> m2_vertex_id(num_vertices(*m2), -1);

  //iterate both meshes to set ids to all points, and set vertex/point_id maps.
  std::size_t id =0;
  for(auto v : m1->vertices())
  {
    Point_3 p = m1->point(v);
    auto res = point_id_map.insert(std::make_pair(p, id));
    if(res.second)
      id++;
    m1_vertex_id[(std::size_t)v]=res.first->second;
  }
  for(auto v : m2->vertices())
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
    for(auto v : CGAL::vertices_around_face(halfedge(f, *m1), *m1))
    {
      ids.push_back(m1_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    m1_faces.insert(ids);
  }

  std::vector<SMesh::Face_index> common_faces;
  std::vector<Point_3> common_points;
  std::map<SMesh::Vertex_index, std::size_t> id_map;
  id = 0;

  for(auto f : m2->faces())
  {
    std::vector<std::size_t> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, *m2), *m2))
    {
      ids.push_back(m2_vertex_id[(std::size_t)v]);
    }
    std::sort(ids.begin(), ids.end());
    if(!((m1_faces.find(ids) != m1_faces.end()) ^ compute_common))
    {
      common_faces.push_back(f);
      for(auto v : CGAL::vertices_around_face(halfedge(f, *m2), *m2))
      {
        auto res = id_map.insert(std::make_pair(v,id));
        if(res.second)
        {
          common_points.push_back(m2->point(v));
          id++;
        }
      }
    }
  }

  //iterate m1_faces and fill a polygon vector using the id_map previously filled.
  std::vector<std::vector<std::size_t> > polygons(common_faces.size());
  id = 0;
  for(auto f : common_faces)
  {
    for(auto v : vertices_around_face(halfedge(f, *m2),*m2))
    {
      polygons[id].push_back(id_map[v]);
    }
    ++id;
  }

  SMesh* common = new SMesh();
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<SMesh>(
    common_points, polygons, *common);
  return common;

}


void Polyhedron_demo_diff_between_meshes_plugin::diff()
{

  QCursor c(Qt::WaitCursor);
  CGAL::Three::Three::CursorScopeGuard guard(c);

  //Get the two meshes. No need to check their existance, applicable()
  //is not permissive enough to let it crash.
  Scene_surface_mesh_item* m1_item = qobject_cast<Scene_surface_mesh_item*>(
        scene->item(scene->selectionIndices().front())),
      *m2_item = qobject_cast<Scene_surface_mesh_item*>(
        scene->item(scene->selectionIndices().back()));

  SMesh* m1=m1_item->face_graph(),
      *m2=m2_item->face_graph();

  SMesh* m1_over_m2 = diff(m1, m2);
  SMesh* m2_over_m1 = diff(m2, m1);
  SMesh* common = diff(m2, m1, true);

  Scene_surface_mesh_item* m1_over_m2_item = new Scene_surface_mesh_item(m1_over_m2);
  m1_over_m2_item->setColor(QColor(Qt::blue));
  m1_over_m2_item->setName(QString("%1 - %2").arg(m1_item->name()).arg(m2_item->name()));
  CGAL::Three::Three::scene()->addItem(m1_over_m2_item);
  Scene_surface_mesh_item* m2_over_m1_item = new Scene_surface_mesh_item(m2_over_m1);
  m2_over_m1_item->setColor(QColor(Qt::red));
  m2_over_m1_item->setName(QString("%2 - %1").arg(m1_item->name()).arg(m2_item->name()));
  CGAL::Three::Three::scene()->addItem(m2_over_m1_item);
  Scene_surface_mesh_item* common_item = new Scene_surface_mesh_item(common);
  common_item->setColor(QColor(Qt::green));
  CGAL::Three::Three::scene()->addItem(common_item);
  common_item->setName(QString("%1 && %2").arg(m1_item->name()).arg(m2_item->name()));

  Scene_group_item* group = new Scene_group_item();
  group->setName("Diff result");
  CGAL::Three::Three::scene()->addItem(group);
  CGAL::Three::Three::scene()->changeGroup(m1_over_m2_item, group);
  CGAL::Three::Three::scene()->changeGroup(m2_over_m1_item, group);
  CGAL::Three::Three::scene()->changeGroup(common_item, group);

  m1_item->setVisible(false);
  m2_item->setVisible(false);

}


#include "Diff_between_meshes_plugin.moc"
