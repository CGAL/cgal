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

#include <CGAL/remove_outliers.h>

#include "ui_Point_set_demo_outlier_removal_plugin.h"

class Point_set_demo_cleaning_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);
  QAction* actionOutlierRemoval;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionOutlierRemoval = this->getActionFromMainWindow(mw, "actionOutlierRemoval");
    if(actionOutlierRemoval) {
      connect(actionOutlierRemoval, SIGNAL(triggered()),
              this, SLOT(on_actionOutlierRemoval_triggered()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionOutlierRemoval;
  }

public slots:
  void on_actionOutlierRemoval_triggered();

}; // end Point_set_demo_cleaning_plugin

class Point_set_demo_outlier_removal_dialog : public QDialog, private Ui::OutlierRemovalDialog
{
  Q_OBJECT
  public:
    Point_set_demo_outlier_removal_dialog(QWidget *parent = 0)
    {
      setupUi(this);
    }

    double percentage() const { return m_inputPercentage->value(); }
    int nbNeighbors() const { return m_inputNbNeighbors->value(); }
};

void Point_set_demo_cleaning_plugin::on_actionOutlierRemoval_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

  if(item)
  {
    Point_set* points = item->point_set();

    if(!points) return;

    Point_set_demo_outlier_removal_dialog dialog;
    if(!dialog.exec())
      return;

    bool areOriented = points->unoriented_points_begin() == points->end();
    const double removed_percentage = dialog.percentage(); // percentage of points to remove
    const int nb_neighbors = dialog.nbNeighbors(); // considers 7 nearest neighbor points
    points->erase(CGAL::remove_outliers(points->begin(), points->end(),
                                        nb_neighbors,
                                        removed_percentage),
                  points->end());

    // Scott Meyer's "swap trick" to trim excess capacity
    Point_set(points).swap(points);

    points->invalidate_bounds();
    if (areOriented)
      points->unoriented_points_begin() = points->end();
    else
      points->unoriented_points_begin() = points->begin();
    item->changed();

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Point_set_demo_cleaning_plugin, Point_set_demo_cleaning_plugin);

#include "Point_set_demo_cleaning_plugin.moc"
