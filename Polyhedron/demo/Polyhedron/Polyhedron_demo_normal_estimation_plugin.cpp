#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
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

#include "ui_Polyhedron_demo_normal_estimation_plugin.h"

#if BOOST_VERSION == 105700
#if !defined(BOOST_NO_CXX11_RVALUE_REFERENCES) && !defined(BOOST_NO_CXX11_DEFAULTED_FUNCTIONS)
#  define CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN 1
#endif
#endif

class Polyhedron_demo_normal_estimation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionNormalEstimation;
  QAction* actionNormalInversion;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionNormalEstimation = new QAction(tr("Normal estimation of point set"), mainWindow);
    actionNormalEstimation->setObjectName("actionNormalEstimation");

    actionNormalInversion = new QAction(tr("Inverse normal orientation"), mainWindow);
    actionNormalInversion->setObjectName("actionNormalInversion");
    
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionNormalEstimation << actionNormalInversion;
  }

  bool applicable(QAction* action) const {
    Scene_points_with_normal_item* item = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));

    if (action==actionNormalEstimation)
#if CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN
    return false;
#else
    return item;
#endif
    else
      return item && item->has_normals();
  }

public Q_SLOTS:
  void on_actionNormalEstimation_triggered();
  void on_actionNormalInversion_triggered();

}; // end PS_demo_smoothing_plugin

class Point_set_demo_normal_estimation_dialog : public QDialog, private Ui::NormalEstimationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_normal_estimation_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
    }

    QString directionMethod() const { return m_inputDirection->currentText(); }
    int directionNbNeighbors() const { return m_inputNbNeighborsDirection->value(); }

    QString orientationMethod() const { return m_inputOrientation->currentText(); }
    int orientationNbNeighbors() const { return m_inputNbNeighborsOrientation->value(); }
};


void Polyhedron_demo_normal_estimation_plugin::on_actionNormalInversion_triggered()
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
  
    for(Point_set::iterator it = points->begin(); it != points->end(); ++it){
      it->normal() = -1 * it->normal();
    }
    item->invalidate_buffers();
    scene->itemChanged(item);
  }
}

void Polyhedron_demo_normal_estimation_plugin::on_actionNormalEstimation_triggered()
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
    // normal estimation
    //***************************************

    if (dialog.directionMethod() == "plane")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Estimates normal direction by PCA (k=" << dialog.directionNbNeighbors() <<")...\n";

      // Estimates normals direction.
      CGAL::pca_estimate_normals(points->begin(), points->end(),
                                CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                                dialog.directionNbNeighbors());

      // Mark all normals as unoriented
      first_unoriented_point = points->begin();

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
                                                  << (memory>>20) << " Mb allocated"
                                                  << std::endl;
    }
    else if (dialog.directionMethod() == "quadric")
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Estimates normal direction by Jet Fitting (k=" << dialog.directionNbNeighbors() <<")...\n";

      // Estimates normals direction.
      CGAL::jet_estimate_normals(points->begin(), points->end(),
                                CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                                dialog.directionNbNeighbors());

      // Mark all normals as unoriented
      first_unoriented_point = points->begin();

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
                                                  << (memory>>20) << " Mb allocated"
                                                  << std::endl;
    }

    //***************************************
    // normal orientation
    //***************************************

    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Orient normals with a Minimum Spanning Tree (k=" << dialog.orientationNbNeighbors() << ")...\n";

    // Tries to orient normals
    first_unoriented_point =
      CGAL::mst_orient_normals(points->begin(), points->end(),
                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
                              dialog.orientationNbNeighbors());

    //indicates that the point set has normals
    if (first_unoriented_point!=points->begin()){
      item->set_has_normals(true);
      item->setRenderingMode(PointsPlusNormals);
    }

    std::size_t nb_unoriented_normals = std::distance(first_unoriented_point, points->end());
    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Orient normals: " << nb_unoriented_normals << " point(s) with an unoriented normal are selected ("
                                    << task_timer.time() << " seconds, "
                                    << (memory>>20) << " Mb allocated)"
                                    << std::endl;

    // Selects points with an unoriented normal
    points->select(points->begin(), points->end(), false);
    points->select(first_unoriented_point, points->end(), true);

    // Updates scene
    item->invalidate_buffers();
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();

    // Warns user
    if (nb_unoriented_normals > 0)
    {
      QMessageBox::information(NULL,
                               tr("Points with an unoriented normal"),
                               tr("%1 point(s) with an unoriented normal are selected.\nPlease orient them or remove them before running Poisson reconstruction.")
                               .arg(nb_unoriented_normals));
    }
  }
#endif // !CGAL_DISABLE_NORMAL_ESTIMATION_PLUGIN
}

#include "Polyhedron_demo_normal_estimation_plugin.moc"
