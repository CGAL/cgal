#include <QTime>
#include <QApplication>

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Subdivision_method_3.h>

class Polyhedron_demo_subdivision_methods_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionLoop"
                         << "actionCatmullClark"
                         << "actionSqrt3";
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }
public Q_SLOTS:
  void on_actionLoop_triggered();
  void on_actionCatmullClark_triggered();
  void on_actionSqrt3_triggered();
}; // end Polyhedron_demo_subdivision_methods_plugin

void Polyhedron_demo_subdivision_methods_plugin::on_actionLoop_triggered()
{
  Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  QTime time;
  time.start();
  std::cout << "Loop subdivision...";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Loop_subdivision(*poly, 1);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
  item->invalidate_buffers();
  scene->itemChanged(item);
}

void Polyhedron_demo_subdivision_methods_plugin::on_actionCatmullClark_triggered()
{
  Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  if(!poly) return;
  QTime time;
  time.start();
  std::cout << "Catmull-Clark subdivision...";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::CatmullClark_subdivision(*poly, 1);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
  item->invalidate_buffers();
  scene->itemChanged(item);
}

void Polyhedron_demo_subdivision_methods_plugin::on_actionSqrt3_triggered()
{
  Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;
  Polyhedron* poly = item->polyhedron();
  if(!poly) return;
  QTime time;
  time.start();
  std::cout << "Sqrt3 subdivision...";
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Subdivision_method_3::Sqrt3_subdivision(*poly, 1);
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  QApplication::restoreOverrideCursor();
  item->invalidate_buffers();
  scene->itemChanged(item);
}

#include "Polyhedron_demo_subdivision_methods_plugin.moc"
