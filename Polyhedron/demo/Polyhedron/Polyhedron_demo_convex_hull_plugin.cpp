#include <QTime>
#include <QApplication>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/convex_hull_3.h>

class Polyhedron_demo_convex_hull_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionConvexHull";
  }

public slots:
  void on_actionConvexHull_triggered();

}; // end Polyhedron_demo_convex_hull_plugin


void Polyhedron_demo_convex_hull_plugin::on_actionConvexHull_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    QTime time;
    time.start();
    std::cout << "Convex hull...";

    // add convex hull as new polyhedron
    Polyhedron *pConvex_hull = new Polyhedron;
    CGAL::convex_hull_3(pMesh->points_begin(),pMesh->points_end(),*pConvex_hull);
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pConvex_hull);
    new_item->setName(tr("%1 (convex hull)").arg(item->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(item->renderingMode());
    scene->addItem(new_item);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_convex_hull_plugin, Polyhedron_demo_convex_hull_plugin)

#include "Polyhedron_demo_convex_hull_plugin.moc"
