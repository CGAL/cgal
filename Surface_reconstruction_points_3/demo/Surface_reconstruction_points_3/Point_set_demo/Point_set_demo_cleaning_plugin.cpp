#include "config.h"
#include "Point_set_scene_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/remove_outliers.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

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
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    Point_set_demo_outlier_removal_dialog dialog;
    if(!dialog.exec())
      return;
    const double removed_percentage = dialog.percentage(); // percentage of points to remove
    const int nb_neighbors = dialog.nbNeighbors(); 
      
    // First point to delete
    Point_set::iterator first_point_to_remove = points->end();

    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Remove outliers (" << removed_percentage <<"%)...\n";

    // Computes outliers
    first_point_to_remove =
      CGAL::remove_outliers(points->begin(), points->end(),
                            nb_neighbors,
                            removed_percentage);

    long memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Remove outliers: " << task_timer.time() << " seconds, "
                                     << (memory>>20) << " Mb allocated"
                                     << std::endl;
                          
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

Q_EXPORT_PLUGIN2(Point_set_demo_cleaning_plugin, Point_set_demo_cleaning_plugin);

#include "Point_set_demo_cleaning_plugin.moc"
