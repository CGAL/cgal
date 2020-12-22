#include <QApplication>
#include <QAction>
#include <QStringList>
#include <QMainWindow>
#include <QInputDialog>
#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
using namespace CGAL::Three;
class Polyhedron_demo_inside_out_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "inside_out_plugin.json")

public:

  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
      scene = scene_interface;
      this->mw = mw;
      actionInsideOut = new QAction(tr("Inside Out"), mw);

      actionInsideOut->setProperty("subMenuName", "Polygon Mesh Processing");
      connect(actionInsideOut, SIGNAL(triggered()), this, SLOT(on_actionInsideOut_triggered()));
      _actions << actionInsideOut;

      actionOrientCC = new QAction(tr("Orient Connected Components"), mw);

      actionOrientCC->setProperty("subMenuName", "Polygon Mesh Processing");
      connect(actionOrientCC, SIGNAL(triggered()), this, SLOT(on_actionOrientCC_triggered()));
      _actions << actionOrientCC;


  }
  bool applicable(QAction* action) const {
    const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
    if(action == actionInsideOut)
      return qobject_cast<Scene_polygon_soup_item*>(scene->item(index))
          || qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
    else if(action == actionOrientCC)
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
    return false;
  }

public Q_SLOTS:
  void on_actionInsideOut_triggered();
  void on_actionOrientCC_triggered();

private:
  QAction* actionInsideOut;
  QAction* actionOrientCC;
  QList<QAction*> _actions;
  Scene_interface *scene;
  QMainWindow* mw;
}; // end Polyhedron_demo_inside_out_plugin

void Polyhedron_demo_inside_out_plugin::on_actionInsideOut_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polygon_soup_item* soup_item =
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(soup_item || sm_item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(sm_item) {
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

void Polyhedron_demo_inside_out_plugin::on_actionOrientCC_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(sm_item)
  {
    QStringList items;
    items << tr("Outward") << tr("Inward");

    bool ok;
    QString item = QInputDialog::getItem(mw, tr("QInputDialog::getItem()"),
                                         tr("The connected components should be oriented:"), items, 0, false, &ok);
    if (!ok )
      return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    if(sm_item) {
      SMesh* pMesh = sm_item->polyhedron();
      if(pMesh){
        if(is_closed(*pMesh))
          CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(*pMesh);
        else
          CGAL::Polygon_mesh_processing::orient(*pMesh);
        sm_item->invalidateOpenGLBuffers();
      }
    }

    // update scene
    scene->itemChanged(index);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Inside_out_plugin.moc"
