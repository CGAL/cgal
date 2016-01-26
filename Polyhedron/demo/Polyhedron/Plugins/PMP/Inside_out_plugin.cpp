#include <QApplication>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
using namespace CGAL::Three;
class Polyhedron_demo_inside_out_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionInsideOut";
  }

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface)
  {
      mw = mainWindow;
      scene = scene_interface;
      actions_map["actionInsideOut"] = getActionFromMainWindow(mw, "actionInsideOut");
      actions_map["actionInsideOut"]->setProperty("subMenuName", "Polygon Mesh Processing");
      autoConnectActions();

  }
  bool applicable(QAction*) const { 
    const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
    return qobject_cast<Scene_polyhedron_item*>(scene->item(index)) 
      || qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
  }

public Q_SLOTS:
  void on_actionInsideOut_triggered();

}; // end Polyhedron_demo_inside_out_plugin

void Polyhedron_demo_inside_out_plugin::on_actionInsideOut_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_polygon_soup_item* soup_item = 
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(poly_item || soup_item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(poly_item) {
      Polyhedron* pMesh = poly_item->polyhedron();
      if(!pMesh) return;
  
      // inside out
      CGAL::Polygon_mesh_processing::reverse_face_orientations(*pMesh);
      poly_item->invalidateOpenGLBuffers();
    }
    else {
      soup_item->inside_out();
      soup_item->invalidateOpenGLBuffers();
    }

    // update scene
    scene->itemChanged(index);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Inside_out_plugin.moc"
