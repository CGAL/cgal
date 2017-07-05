#include <QTime>
#include <QApplication>
#include <QMainWindow>
#include <QAction>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/subdivision_method_3.h>
using namespace CGAL::Three;
class Polyhedron_demo_subdivision_methods_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  // used by Polyhedron_demo_plugin_helper
  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface* m)
  {
      messages = m;
      scene = scene_interface;
      QAction *actionLoop = new QAction("Loop", mw);
      QAction *actionCatmullClark = new QAction("Catmull Clark", mw);
      QAction *actionSqrt3 = new QAction("Sqrt3", mw);
      QAction *actionDooSabin = new QAction("Doo Sabin", mw);
      actionLoop->setObjectName("actionLoop");
      actionCatmullClark->setObjectName("actionCatmullClark");
      actionSqrt3->setObjectName("actionSqrt3");
      actionDooSabin->setObjectName("actionDooSabin");
      _actions << actionLoop
               << actionCatmullClark
               << actionSqrt3
               << actionDooSabin;
      Q_FOREACH(QAction* action, _actions)
        action->setProperty("subMenuName", "3D Surface Subdivision Methods");
      autoConnectActions();

  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }
public Q_SLOTS:
  void on_actionLoop_triggered();
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();
  void on_actionDooSabin_triggered();
private :
  Messages_interface* messages;
  QList<QAction*> _actions;
}; // end Polyhedron_demo_subdivision_methods_plugin

void Polyhedron_demo_subdivision_methods_plugin::on_actionLoop_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  QTime time;
  time.start();
  messages->information("Loop subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Loop_subdivision(*poly, 1);
  messages->information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void Polyhedron_demo_subdivision_methods_plugin::on_actionCatmullClark_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  if(!poly) return;
  QTime time;
  time.start();
  messages->information("Catmull-Clark subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::CatmullClark_subdivision(*poly, 1);
  messages->information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void Polyhedron_demo_subdivision_methods_plugin::on_actionSqrt3_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  if(!poly) return;
  QTime time;
  time.start();
  messages->information("Sqrt3 subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Sqrt3_subdivision(*poly, 1);
  messages->information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

void Polyhedron_demo_subdivision_methods_plugin::on_actionDooSabin_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  if(!poly) return;
  QTime time;
  time.start();
  messages->information("Doo Sabin subdivision...");
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::DooSabin_subdivision(*poly, 1);
  messages->information(QString("ok (%1 ms)").arg(time.elapsed()));
  QApplication::restoreOverrideCursor();
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}

#include "Subdivision_methods_plugin.moc"
