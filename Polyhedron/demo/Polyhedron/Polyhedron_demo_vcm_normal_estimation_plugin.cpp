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

#if BOOST_VERSION == 105700
#if !defined(BOOST_NO_CXX11_RVALUE_REFERENCES) && !defined(BOOST_NO_CXX11_DEFAULTED_FUNCTIONS)
#  define CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN 1
#endif
#endif

class Polyhedron_demo_vcm_normal_estimation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  QAction* actionVCMNormalEstimation;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionVCMNormalEstimation = new QAction(tr("VCM normal estimation"), mainWindow);
    actionVCMNormalEstimation->setObjectName("actionVCMNormalEstimation");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionVCMNormalEstimation;
  }

  bool applicable(QAction*) const {
#if CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN
    return false;
#else
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
#endif
  }

public Q_SLOTS:
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
    unsigned int convolveNeighbors() const { return m_inputConvolveNeighbors->value(); }
    bool convolveUsingRadius() const { return m_inputUseConvRad->isChecked(); }
};

void Polyhedron_demo_vcm_normal_estimation_plugin::on_actionVCMNormalEstimation_triggered()
{
#if !CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN
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

    // First point to delete
    Point_set::iterator first_unoriented_point = points->end();

    //***************************************
    // VCM normal estimation
    //***************************************

    CGAL::Timer task_timer; task_timer.start();

    // Estimates normals direction.
    if (dialog.convolveUsingRadius())
    {
      std::cerr << "Estimates Normals Direction using VCM (R="
          << dialog.offsetRadius() << " and r=" << dialog.convolveRadius() << ")...\n";

      CGAL::vcm_estimate_normals(points->begin(), points->end(),
                                 CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                                 dialog.offsetRadius(), dialog.convolveRadius());
    }
    else
    {
      std::cerr << "Estimates Normals Direction using VCM (R="
        << dialog.offsetRadius() << " and k=" << dialog.convolveNeighbors() << ")...\n";

      CGAL::vcm_estimate_normals(points->begin(), points->end(),
                                 CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                                 dialog.offsetRadius(), dialog.convolveNeighbors());
    }
      // Mark all normals as unoriented
      first_unoriented_point = points->begin();

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    task_timer.stop();
    std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
        << (memory>>20) << " Mb allocated"
        << std::endl;

    //***************************************
    // normal orientation
    //***************************************

    unsigned int neighbors = 18;
    task_timer.reset();
    task_timer.start();
    std::cerr << "Orient normals with a Minimum Spanning Tree (k=" << neighbors << ")...\n";

    // Tries to orient normals
    first_unoriented_point = CGAL::mst_orient_normals(points->begin(), points->end(),
                                                      CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                                                      neighbors);

    std::size_t nb_unoriented_normals = std::distance(first_unoriented_point, points->end());
    memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Orient normals: " << nb_unoriented_normals << " point(s) with an unoriented normal are selected ("
                                    << task_timer.time() << " seconds, "
                                    << (memory>>20) << " Mb allocated)"
                                    << std::endl;

    // Updates scene
    item->invalidate_buffers();
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();
  }
#endif // !CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN
}

#include "Polyhedron_demo_vcm_normal_estimation_plugin.moc"
