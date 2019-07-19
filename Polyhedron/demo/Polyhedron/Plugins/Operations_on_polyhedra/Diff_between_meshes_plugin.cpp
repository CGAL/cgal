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
  SMesh* diff(SMesh* m1, SMesh* m2);
  SMesh* common(SMesh* m1, SMesh* m2);

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
SMesh* Polyhedron_demo_diff_between_meshes_plugin::diff(SMesh* m1, SMesh* m2)
{
  //Collect points of both meshes in separate sets to easily detect if a point is
  // contained by a mesh or not.
std::set<Point_3> m2_verts;
  for(auto v : m2->vertices())
  {
    m2_verts.insert(m2->point(v));
  }

  //Get Vertices that are in m1 but not in m2
  std::vector<Point_3> m1_points;

  std::map<SMesh::Vertex_index, std::size_t> id_map;
  std::size_t id = 0;
  std::vector<SMesh::Face_index> m1_faces;

  // parse faces of m1. Select each face that has at least one point that is not
  // in m2. Fill a points vector with them.
  for(auto f : m1->faces())
  {
    bool take = false;
    for(auto v : CGAL::vertices_around_face(halfedge(f, *m1), *m1))
    {
      if(m2_verts.find(m1->point(v)) == m2_verts.end())
      {
        take = true;
        break;
      }
    }
    if(take)
    {
      m1_faces.push_back(f);
      for(auto v : CGAL::vertices_around_face(halfedge(f, *m1), *m1))
      {
        if(id_map.find(v) == id_map.end())
        {
          m1_points.push_back(m1->point(v));
          id_map[v] = id++;
        }
      }
    }
  }

  //iterate m1_faces and fill a polygon vector using the id_map previously filled.
  std::vector<std::vector<std::size_t> > polygons(m1_faces.size());
  id = 0;
  for(auto f : m1_faces)
  {
    for(auto v : vertices_around_face(halfedge(f, *m1),*m1))
    {
      polygons[id].push_back(id_map[v]);
    }
    ++id;
  }

  SMesh* m1_over_m2 = new SMesh();
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<SMesh>(
    m1_points, polygons, *m1_over_m2);
  return m1_over_m2;
}

SMesh* Polyhedron_demo_diff_between_meshes_plugin::common(SMesh* m1, SMesh* m2)
{
  std::set<Point_3> m1_verts;
  for(auto v : m1->vertices())
  {
    m1_verts.insert(m1->point(v));
  }
  std::set<Point_3> m2_verts;
  for(auto v : m2->vertices())
  {
    m2_verts.insert(m2->point(v));
  }

  std::vector<Point_3> m1_points;

  std::map<SMesh::Vertex_index, std::size_t> id_map;
  std::size_t id = 0;
  std::vector<SMesh::Face_index> m1_faces;
  //take all faces of m1 whose points are all in m2
  for(auto f : m1->faces())
  {
    for(auto v : CGAL::vertices_around_face(halfedge(f, *m1), *m1))
    {
      if(m2_verts.find(m1->point(v)) == m2_verts.end())
      {
        m1_faces.push_back(f);
        break;
      }
    }
  }
  //take all faces of m2 whose points are all in m1
  for(auto f : m2->faces())
  {
    for(auto v : CGAL::vertices_around_face(halfedge(f, *m2), *m2))
    {
    }
  }


  std::vector<std::vector<std::size_t> > polygons(m1_faces.size());
  id = 0;
  for(auto f : m1_faces)
  {
    for(auto v : vertices_around_face(halfedge(f, *m1),*m1))
    {
      polygons[id].push_back(id_map[v]);
    }
    ++id;
  }

  SMesh* m1_over_m2 = new SMesh();
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh<SMesh>(
    m1_points, polygons, *m1_over_m2);
  return m1_over_m2;
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
  SMesh* common = common(m2, m1);

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
