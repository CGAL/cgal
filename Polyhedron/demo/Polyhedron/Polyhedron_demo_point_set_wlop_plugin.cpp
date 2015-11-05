#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Polyhedron_demo_point_set_wlop_plugin.h"

class Polyhedron_demo_point_set_wlop_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  
  QAction* actionSimplifyAndRegularize;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionSimplifyAndRegularize = new QAction(tr("Point set WLOP simplification and regularization selection"), mainWindow);
    actionSimplifyAndRegularize->setObjectName("actionSimplifyAndRegularize");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSimplifyAndRegularize;
  }

public Q_SLOTS:
  void on_actionSimplifyAndRegularize_triggered();

}; // end Polyhedron_demo_point_set_wlop_plugin

class Point_set_demo_point_set_wlop_simplification_and_regularization_dialog : public QDialog, private Ui::WLOPRegularizationAndSimplificationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_point_set_wlop_simplification_and_regularization_dialog(QWidget * /*parent*/ = 0)
    {
      setupUi(this);
    }

    double retainedPercentage() const { return m_retainedPercentage->value(); }
    double neighborhoodRadius() const { return m_neighborhoodRadius->value(); }
};

void Polyhedron_demo_point_set_wlop_plugin::on_actionSimplifyAndRegularize_triggered()
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
    Point_set_demo_point_set_wlop_simplification_and_regularization_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    CGAL::Timer task_timer; task_timer.start();

    std::cerr << "Point cloud simplification and regularization by WLOP ("
	      << dialog.retainedPercentage () << "% retained points, neighborhood radius = "
	      << dialog.neighborhoodRadius() <<" * average spacing)...\n";

    // Computes average spacing
    double average_spacing = CGAL::compute_average_spacing(points->begin(), points->end(),
							   6 /* knn = 1 ring */);

    Scene_points_with_normal_item* new_item
      = new Scene_points_with_normal_item();
    new_item->setName (tr("%1 (WLOP processed)").arg(item->name()));
    new_item->setVisible(true);
    scene->addItem(new_item);
      
    CGAL::wlop_simplify_and_regularize_point_set<CGAL::Parallel_tag> // parallel version
      (points->begin(),
       points->end(),
       std::back_inserter(*(new_item->point_set ())),
       dialog.retainedPercentage (),
       dialog.neighborhoodRadius()*average_spacing
       );

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Simplification and regularization: "
	      << new_item->point_set ()->size () << " point(s) created ("
	      << task_timer.time() << " seconds, "
	      << (memory>>20) << " Mb allocated)"
	      << std::endl;

    QApplication::restoreOverrideCursor();
  }
}

#include "Polyhedron_demo_point_set_wlop_plugin.moc"
