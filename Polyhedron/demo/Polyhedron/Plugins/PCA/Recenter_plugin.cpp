#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/centroid.h>


#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

using namespace CGAL::Three;
class Polyhedron_demo_recenter_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionRecenter;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionRecenter = new QAction(tr("Recenter"), mainWindow);
    actionRecenter->setObjectName("actionRecenter");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRecenter;
  }
  
  //! Applicable if the currently selected item is a
  //! points_with_normal_item (can be extended to other items).
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionRecenter_triggered();
}; // end Polyhedron_demo_recenter_plugin

void Polyhedron_demo_recenter_plugin::on_actionRecenter_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    Kernel::Point_3 centroid = CGAL::centroid (points->points().begin(),
                                               points->points().end());
    Kernel::Vector_3 recenter (centroid, CGAL::ORIGIN);
    std::cerr << "Recenter by " << recenter << std::endl;
    for (Point_set::iterator it = points->begin(); it != points->end(); ++ it)
      points->point(*it) = points->point(*it) + recenter;

    QApplication::restoreOverrideCursor();
    points->invalidate_bounds();
    item->invalidateOpenGLBuffers();
    scene->itemChanged(index);
  }
}


#include "Recenter_plugin.moc"
