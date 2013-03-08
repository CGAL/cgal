#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Polyhedron_demo_point_set_simplification_plugin.h"

class Polyhedron_demo_point_set_simplification_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionSimplify;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionSimplify = new QAction(tr("Point set simplification selection"), mainWindow);
    actionSimplify->setObjectName("actionSimplify");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable() const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSimplify;
  }

public slots:
  void on_actionSimplify_triggered();

}; // end Polyhedron_demo_point_set_simplification_plugin

class Point_set_demo_point_set_simplification_dialog : public QDialog, private Ui::PointSetSimplificationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_point_set_simplification_dialog(QWidget * /*parent*/ = 0)
    {
      setupUi(this);
    }

    QString simplificationMethod() const { return m_simplificationMethod->currentText(); }
    double randomSimplificationPercentage() const { return m_randomSimplificationPercentage->value(); }
    double gridCellSize() const { return m_gridCellSize->value(); }
};

void Polyhedron_demo_point_set_simplification_plugin::on_actionSimplify_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    Point_set_demo_point_set_simplification_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    CGAL::Timer task_timer; task_timer.start();

    // First point to delete
    Point_set::iterator first_point_to_remove = points->end();

    if (dialog.simplificationMethod() == "Random")
    {
      std::cerr << "Random point cloud simplification (" << dialog.randomSimplificationPercentage() <<"%)...\n";

      // Computes points to remove by random simplification
      first_point_to_remove =
        CGAL::random_simplify_point_set(points->begin(), points->end(),
                                        dialog.randomSimplificationPercentage());
    }
    else if (dialog.simplificationMethod() == "Grid Clustering")
    {
      std::cerr << "Point cloud simplification by clustering (cell size = " << dialog.gridCellSize() <<" * average spacing)...\n";

      // Computes average spacing
      double average_spacing = CGAL::compute_average_spacing(
                                      points->begin(), points->end(),
                                      6 /* knn = 1 ring */);

      // Computes points to remove by Grid Clustering
      first_point_to_remove =
        CGAL::grid_simplify_point_set(points->begin(), points->end(),
                                      dialog.gridCellSize()*average_spacing);
    }

    std::size_t nb_points_to_remove = std::distance(first_point_to_remove, points->end());
    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Simplification: " << nb_points_to_remove << " point(s) are selected for removal ("
                                    << task_timer.time() << " seconds, "
                                    << (memory>>20) << " Mb allocated)"
                                    << std::endl;

    // Selects points to delete
    points->select(points->begin(), points->end(), false);
    points->select(first_point_to_remove, points->end(), true);

    // Updates scene
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();

    // Warns user
    if (nb_points_to_remove > 0)
    {
      QMessageBox::information(NULL,
                               tr("Points selected for removal"),
                               tr("%1 point(s) are selected for removal.\nYou may delete or reset the selection using the item context menu.")
                               .arg(nb_points_to_remove));
    }
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_point_set_simplification_plugin, Polyhedron_demo_point_set_simplification_plugin)

#include "Polyhedron_demo_point_set_simplification_plugin.moc"
