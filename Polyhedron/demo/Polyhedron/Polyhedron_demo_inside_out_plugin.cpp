#include <QApplication>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

class Polyhedron_demo_inside_out_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionInsideOut";
  }

  bool applicable() const { 
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    return qobject_cast<Scene_polyhedron_item*>(scene->item(index)) 
      || qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
  }

public slots:
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
      pMesh->inside_out();
    }
    else {
      soup_item->inside_out();
    }

    // update scene
    scene->itemChanged(index);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_inside_out_plugin, Polyhedron_demo_inside_out_plugin)

#include "Polyhedron_demo_inside_out_plugin.moc"
