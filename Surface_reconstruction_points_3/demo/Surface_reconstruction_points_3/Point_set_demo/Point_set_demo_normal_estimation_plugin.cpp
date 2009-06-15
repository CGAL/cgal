#include "config.h"
#include "Point_set_scene_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include <CGAL/Timer.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

#include "ui_Point_set_demo_normal_estimation_plugin.h"

class Point_set_demo_normal_estimation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);
  QAction* actionNormalEstimation;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionNormalEstimation = this->getActionFromMainWindow(mw, "actionNormalEstimation");
    if(actionNormalEstimation) {
      connect(actionNormalEstimation, SIGNAL(triggered()),
              this, SLOT(on_actionNormalEstimation_triggered()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionNormalEstimation;
  }

public slots:
  void on_actionNormalEstimation_triggered();

}; // end Point_set_demo_smoothing_plugin

class Point_set_demo_normal_estimation_dialog : public QDialog, private Ui::NormalEstimationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_normal_estimation_dialog(QWidget *parent = 0)
    {
      setupUi(this);
    }

    QString directionMethod() const { return m_inputDirection->currentText(); }
    int directionNbNeighbors() const { return m_inputNbNeighborsDirection->value(); }

    QString orientationMethod() const { return m_inputOrientation->currentText(); }
    int orientationNbNeighbors() const { return m_inputNbNeighborsOrientation->value(); }
};

void Point_set_demo_normal_estimation_plugin::on_actionNormalEstimation_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

  if(item)
  {
    // get point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    Point_set_demo_normal_estimation_dialog dialog;
    if(!dialog.exec())
      return;

    // normal estimation
    if (dialog.directionMethod() == "plane")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Estimates normals directions by PCA (k=" << dialog.directionNbNeighbors() <<")...\n";

      // Estimates normals direction.
      // Note: pca_estimate_normals() requires an iterator over points
      // + property maps to access each point's position and normal.
      // The position property map can be omitted here as we use iterators over Point_3 elements.
      CGAL::pca_estimate_normals(points->begin(), points->end(),
                                CGAL::make_normal_vector_property_map(points->begin()),
                                dialog.directionNbNeighbors());

      // Mark all normals as unoriented
      points->unoriented_points_begin() = points->begin();

      long memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "done: " << task_timer.time() << " seconds, "
                            << (memory>>20) << " Mb allocated"
                            << std::endl;
    }
    else if (dialog.directionMethod() == "quadric")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Estimates normal directions by Jet Fitting (k=" << dialog.directionNbNeighbors() <<")...\n";

      // Estimates normals direction.
      // Note: jet_estimate_normals() requires an iterator over points
      // + property maps to access each point's position and normal.
      // The position property map can be omitted here as we use iterators over Point_3 elements.
      CGAL::jet_estimate_normals(points->begin(), points->end(),
                                CGAL::make_normal_vector_property_map(points->begin()),
                                dialog.directionNbNeighbors());

      // Mark all normals as unoriented
      points->unoriented_points_begin() = points->begin();

      long memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "done: " << task_timer.time() << " seconds, "
                            << (memory>>20) << " Mb allocated"
                            << std::endl;
    }

    // normal orientation
    if (dialog.orientationMethod() == "MST")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Orient normals with a Minimum Spanning Tree (k=" << dialog.orientationNbNeighbors() << ")...\n";

      // Orient normals.
      // Note: mst_orient_normals() requires an iterator over points
      // + property maps to access each point's position and normal.
      // The position property map can be omitted here as we use iterators over Point_3 elements.
      points->unoriented_points_begin() =
        CGAL::mst_orient_normals(points->begin(), points->end(),
                                CGAL::make_normal_vector_property_map(points->begin()),
                                dialog.orientationNbNeighbors());

      // Delete points with an unoriented normal (required by APSS and Poisson)
      points->erase(points->unoriented_points_begin(), points->end());
      points->unoriented_points_begin() = points->end();

      // After erase(), use Scott Meyer's "swap trick" to trim excess capacity
      Point_set(*points).swap(*points);

      long memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "done: " << task_timer.time() << " seconds, "
                            << (memory>>20) << " Mb allocated"
                            << std::endl;
    }

    points->invalidate_bounds();
    item->changed();

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Point_set_demo_normal_estimation_plugin, Point_set_demo_normal_estimation_plugin);

#include "Point_set_demo_normal_estimation_plugin.moc"
