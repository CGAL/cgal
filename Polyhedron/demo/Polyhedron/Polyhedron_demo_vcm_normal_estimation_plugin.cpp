#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

#include "ui_Polyhedron_demo_vcm_normal_estimation_plugin.h"

class Polyhedron_demo_vcm_normal_estimation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionVCMNormalEstimation;
  QAction* actionNormalInversion;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionVCMNormalEstimation = new QAction(tr("Normal estimation of point set using VCM"), mainWindow);
    actionVCMNormalEstimation->setObjectName("actionVCMNormalEstimation");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionVCMNormalEstimation;
  }

  bool applicable() const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public slots:
  void on_actionVCMNormalEstimation_triggered();

}; // end Polyhedron_demo_vcm_normal_estimation_plugin

class Point_set_demo_normal_estimation_dialog : public QDialog, private Ui::VCMNormalEstimationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_normal_estimation_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
    }

    float offsetRadius() const { return m_inputOffsetRadius->value(); }
    float convolveRadius() const { return m_inputConvolveRadius->value(); }
};

void Polyhedron_demo_vcm_normal_estimation_plugin::on_actionVCMNormalEstimation_triggered()
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
    Point_set_demo_normal_estimation_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    //***************************************
    // VCM normal estimation
    //***************************************

    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Estimates Normals Direction using VCM (R="
        << dialog.offsetRadius() << " and r=" << dialog.convolveRadius() << ")...\n";

    // Estimates normals direction.
    CGAL::vcm_estimate_normals(points->begin(), points->end(),
                               CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                               dialog.offsetRadius(), dialog.convolveRadius());

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
        << (memory>>20) << " Mb allocated"
        << std::endl;

    // Updates scene
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_vcm_normal_estimation_plugin, Polyhedron_demo_vcm_normal_estimation_plugin)

#include "Polyhedron_demo_vcm_normal_estimation_plugin.moc"
