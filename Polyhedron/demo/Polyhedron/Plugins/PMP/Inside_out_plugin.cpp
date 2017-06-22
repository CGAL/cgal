#include <QApplication>
#include <QAction>
#include <QStringList>
#include <QMainWindow>
#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
using namespace CGAL::Three;
class Polyhedron_demo_inside_out_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
      scene = scene_interface;
      QAction* actionInsideOut = new QAction(tr("Inside Out"), mw);

      actionInsideOut->setProperty("subMenuName", "Polygon Mesh Processing");
      connect(actionInsideOut, SIGNAL(triggered()), this, SLOT(on_actionInsideOut_triggered()));
      _actions << actionInsideOut;

  }
  bool applicable(QAction*) const { 
    const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
    return qobject_cast<Scene_polyhedron_item*>(scene->item(index)) 
      || qobject_cast<Scene_polygon_soup_item*>(scene->item(index))
      || qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  }

public Q_SLOTS:
  void on_actionInsideOut_triggered();

private:
  QList<QAction*> _actions;
  Scene_interface *scene;
}; // end Polyhedron_demo_inside_out_plugin

void Polyhedron_demo_inside_out_plugin::on_actionInsideOut_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_polygon_soup_item* soup_item = 
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item = 
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(poly_item || soup_item || sm_item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(poly_item) {
      Polyhedron* pMesh = poly_item->polyhedron();
      if(pMesh){
        CGAL::Polygon_mesh_processing::reverse_face_orientations(*pMesh);
        poly_item->invalidateOpenGLBuffers();
      }
    }
    else if(sm_item) {
      SMesh* pMesh = sm_item->polyhedron();
      if(pMesh){
        CGAL::Polygon_mesh_processing::reverse_face_orientations(*pMesh);
        sm_item->invalidateOpenGLBuffers();
      }
    }else{
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
