#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>

#include "Scene_surface_mesh_item.h"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Three.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
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


void Polyhedron_demo_diff_between_meshes_plugin::diff()
{

  typedef CGAL::Face_filtered_graph<SMesh> Filtered_graph;

  QCursor c(Qt::WaitCursor);
  CGAL::Three::Three::CursorScopeGuard guard(c);

  //Get the two meshes. No need to check their existence, applicable()
  //is not permissive enough to let it crash.
  Scene_surface_mesh_item* m1_item = qobject_cast<Scene_surface_mesh_item*>(
        scene->item(scene->selectionIndices().front())),
      *m2_item = qobject_cast<Scene_surface_mesh_item*>(
        scene->item(scene->selectionIndices().back()));

  SMesh m1=*m1_item->face_graph(),
      m2=*m2_item->face_graph();
  std::vector<face_descriptor> m1_only, m2_only;
  std::vector<std::pair<face_descriptor, face_descriptor> > common;
  CGAL::Polygon_mesh_processing::match_faces(
        m1,
        m2,
        std::back_inserter(common),
        std::back_inserter(m1_only),
        std::back_inserter(m2_only));

  Filtered_graph filter1(m1, m1_only);
  SMesh mesh1_only, mesh2_only, common_mesh;
  CGAL::copy_face_graph(filter1, mesh1_only);
  Scene_surface_mesh_item* mesh1_only_item = nullptr;
  if(mesh1_only.faces().size() > 0)
  {
    mesh1_only_item = new Scene_surface_mesh_item(mesh1_only);
    mesh1_only_item->setColor(QColor(Qt::blue));
    mesh1_only_item->setName(QString("%1_only").arg(m1_item->name()));
    CGAL::Three::Three::scene()->addItem(mesh1_only_item);
  }

  Filtered_graph filter2(m2, m2_only);
  CGAL::copy_face_graph(filter2, mesh2_only);
  Scene_surface_mesh_item* mesh2_only_item = nullptr;
  if(mesh2_only.faces().size() > 0)
  {
    mesh2_only_item = new Scene_surface_mesh_item(mesh2_only);
    mesh2_only_item->setColor(QColor(Qt::red));
    mesh2_only_item->setName(QString("%1_only").arg(m2_item->name()));
    CGAL::Three::Three::scene()->addItem(mesh2_only_item);
  }
  m1_only.clear();
  m1_only.reserve(common.size());
  for(const auto& f_pair : common)
  {
    m1_only.push_back(f_pair.first);
  }
  Filtered_graph filter_common(m1, m1_only);
  CGAL::copy_face_graph(filter_common, common_mesh);
  Scene_surface_mesh_item* common_item = nullptr;
  if(common_mesh.faces().size() > 0)
  {
    common_item = new Scene_surface_mesh_item(common_mesh);
    common_item->setColor(QColor(Qt::green));
    CGAL::Three::Three::scene()->addItem(common_item);
    common_item->setName(QString("%1 && %2").arg(m1_item->name()).arg(m2_item->name()));
  }
  Scene_group_item* group = new Scene_group_item();
  group->setName("Diff result");
  CGAL::Three::Three::scene()->addItem(group);
  if(mesh1_only_item)
    CGAL::Three::Three::scene()->changeGroup(mesh1_only_item, group);
  if(mesh2_only_item)
    CGAL::Three::Three::scene()->changeGroup(mesh2_only_item, group);
  if(common_item)
    CGAL::Three::Three::scene()->changeGroup(common_item, group);

  m1_item->setVisible(false);
  m2_item->setVisible(false);

}


#include "Diff_between_meshes_plugin.moc"
