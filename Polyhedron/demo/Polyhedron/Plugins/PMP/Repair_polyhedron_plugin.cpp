#include <QtCore/qglobal.h>

#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Scene_interface.h>
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Messages_interface.h"
#include <CGAL/gl.h>

#include <QAction>
#include <QMainWindow>
#include <QObject>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

using namespace CGAL::Three;
class Polyhedron_demo_repair_polyhedron_plugin :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // To silent a warning -Woverloaded-virtual
  // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
  using Polyhedron_demo_plugin_helper::init;

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface* m)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->messages = m;

    actionRemoveIsolatedVertices = new QAction(tr("Remove isolated vertices"), mw);
    if (actionRemoveIsolatedVertices){
      connect(actionRemoveIsolatedVertices, SIGNAL(triggered()),
              this, SLOT(on_actionRemoveIsolatedVertices_triggered()));
    }

    actionRemoveDegenerateFaces = new QAction(tr("Remove degenerate faces"), mw);
    if (actionRemoveDegenerateFaces){
      connect(actionRemoveDegenerateFaces, SIGNAL(triggered()),
        this, SLOT(on_actionRemoveDegenerateFaces_triggered()));
    }
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionRemoveIsolatedVertices
      << actionRemoveDegenerateFaces;
  }

  bool applicable(QAction*) const
  {
    int item_id = scene->mainSelectionIndex();
    return qobject_cast<Scene_polyhedron_item*>(
      scene->item(item_id));
  }

public Q_SLOTS:
  void on_actionRemoveIsolatedVertices_triggered();
  void on_actionRemoveDegenerateFaces_triggered();

private:
  QAction* actionRemoveIsolatedVertices;
  QAction* actionRemoveDegenerateFaces;

  Messages_interface* messages;
}; // end Polyhedron_demo_repair_polyhedron_plugin


void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveIsolatedVertices_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* poly_item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if (poly_item)
  {
    std::size_t nbv =
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(
        *poly_item->polyhedron());
    messages->information(tr(" %1 isolated vertices have been removed.")
      .arg(nbv));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}

void Polyhedron_demo_repair_polyhedron_plugin::on_actionRemoveDegenerateFaces_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* poly_item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if (poly_item)
  {
    std::size_t nbv =
      CGAL::Polygon_mesh_processing::remove_degenerate_faces(
      *poly_item->polyhedron());
    messages->information(tr(" %1 degenerate faces have been removed.")
      .arg(nbv));
    poly_item->invalidateOpenGLBuffers();
    Q_EMIT poly_item->itemChanged();
  }
}


#include "Repair_polyhedron_plugin.moc"
