#include "config.h"
#include "Point_set_scene_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

#include "ui_Point_set_demo_point_set_simplification_plugin.h"

class Point_set_demo_point_set_simplification_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);
  QAction* actionSimplify;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionSimplify = this->getActionFromMainWindow(mw, "actionSimplify");
    if(actionSimplify) {
      connect(actionSimplify, SIGNAL(triggered()),
              this, SLOT(on_actionSimplify_triggered()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSimplify;
  }

public slots:
  void on_actionSimplify_triggered();

}; // end Point_set_demo_point_set_simplification_plugin

class Point_set_demo_point_set_simplification_dialog : public QDialog, private Ui::PointSetSimplificationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_point_set_simplification_dialog(QWidget *parent = 0)
    {
      setupUi(this);
    }

    QString simplificationMethod() const { return m_simplificationMethod->currentText(); }
    float randomSimplificationPercentage() const { return m_randomSimplificationPercentage->value(); }
    float gridCellSize() const { return m_gridCellSize->value(); }
};

void Point_set_demo_point_set_simplification_plugin::on_actionSimplify_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

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
      
    // First point to delete
    Point_set::iterator first_point_to_remove = points->end();

    if (dialog.simplificationMethod() == "Random")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Random point cloud simplification (" << dialog.randomSimplificationPercentage() <<"%)...\n";

      // Computes points to remove by random simplification
      first_point_to_remove =
        CGAL::random_simplify_point_set(points->begin(), points->end(),
                                        dialog.randomSimplificationPercentage());

      long memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Simplification: " << task_timer.time() << " seconds, "
                                      << (memory>>20) << " Mb allocated"
                                      << std::endl;
    }
    else if (dialog.simplificationMethod() == "Grid Clustering")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Point cloud simplification by clustering (cell size = " << dialog.gridCellSize() <<" * point set radius)...\n";

      // Gets point set's radius
      Sphere bounding_sphere = points->bounding_sphere();
      FT radius = sqrt(bounding_sphere.squared_radius());

      // Computes points to remove by Grid Clustering
      first_point_to_remove =
        CGAL::grid_simplify_point_set(points->begin(), points->end(),
                                      dialog.gridCellSize()*radius);

      long memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Simplification: " << task_timer.time() << " seconds, "
                                      << (memory>>20) << " Mb allocated"
                                      << std::endl;
    }

    // Selects points to delete
    points->select(points->begin(), points->end(), false);
    points->select(first_point_to_remove, points->end(), true);

    // Warns user
    if (first_point_to_remove != points->end())
    {
      int nb_selected_points = std::distance(first_point_to_remove, points->end());
      QMessageBox::information(NULL,
                               tr("Points selected from removal"),
                               tr("%1 point(s) are selected for removal.\nYou may remove them with the \"Delete selection\" menu item.")
                               .arg(nb_selected_points));
    }

    // Updates scene
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Point_set_demo_point_set_simplification_plugin, Point_set_demo_point_set_simplification_plugin);

#include "Point_set_demo_point_set_simplification_plugin.moc"
