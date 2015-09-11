#include <QApplication>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include <CGAL/Polygon_mesh_processing/orientation.h>

class Polyhedron_demo_inside_out_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionInsideOut";
  }

  bool applicable(QAction*) const { 
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    return qobject_cast<Scene_polyhedron_item*>(scene->item(index)) 
      || qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
  }

public Q_SLOTS:
  void on_actionInsideOut_triggered();

}; // end Polyhedron_demo_inside_out_plugin

void Polyhedron_demo_inside_out_plugin::on_actionInsideOut_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
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
      poly_item->invalidate_buffers();
    }
    else {
      soup_item->inside_out();
      soup_item->invalidate_buffers();
    }

    // update scene
    scene->itemChanged(index);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Polyhedron_demo_inside_out_plugin.moc"
