#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>
using namespace CGAL::Three;
class Surface_mesh_to_polyhedron_plugin :
    public QObject,
    public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface* )
  {
    scene = scene_interface;
    this->mw = mw;
    QAction *actionToPoly = new QAction("Convert to Polyhedron", mw);
    QAction *actionToSM = new QAction("Convert to Surface_mesh", mw);
    actionToPoly->setProperty("subMenuName",
                              "BGL");
    actionToSM->setProperty("subMenuName",
                            "BGL");

    connect(actionToPoly, SIGNAL(triggered()), this, SLOT(on_actionToPoly_triggered()));
    connect(actionToSM, SIGNAL(triggered()), this, SLOT(on_actionToSM_triggered()));
    _actions<< actionToPoly
            << actionToSM;

  }

  bool applicable(QAction* action) const {
    if(action == _actions.first())
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    else
      return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return _actions;
  }

private:
  Scene_interface* scene;
  QList<QAction*> _actions;
  QMainWindow* mw;

public Q_SLOTS:
  void on_actionToPoly_triggered()
  {
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!sm_item)
      return;
    Polyhedron* polyhedron = new Polyhedron();
    CGAL::copy_face_graph(*sm_item->polyhedron(), *polyhedron);
    Scene_polyhedron_item* poly_item = new Scene_polyhedron_item(polyhedron);
    poly_item->setColor(sm_item->color());
    poly_item->setName(sm_item->name());
    sm_item->setVisible(false);
    scene->addItem(poly_item);
  }

  void on_actionToSM_triggered()
  {
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
    if(!poly_item)
      return;
    SMesh* sm = new SMesh();
    CGAL::copy_face_graph(*poly_item->polyhedron(), *sm);
    Scene_surface_mesh_item* sm_item = new Scene_surface_mesh_item(sm);
    sm_item->setColor(poly_item->color());
    sm_item->setName(poly_item->name());
    poly_item->setVisible(false);
    scene->addItem(sm_item);
  }
}; // end class Surface_mesh_to_polyhedron_plugin

#include "Surface_mesh_to_polyhedron_plugin.moc"
