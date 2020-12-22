#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include <sstream>

#include "run_with_qprogressdialog.h"

#include "ui_Point_set_bilateral_smoothing_plugin.h"

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

struct Bilateral_smoothing_functor
  : public Functor_with_signal_callback
{
  Point_set* points;
  unsigned int neighborhood_size;
  unsigned int sharpness_angle;
  boost::shared_ptr<double> result;

  Bilateral_smoothing_functor  (Point_set* points,
                                unsigned int neighborhood_size,
                                unsigned int sharpness_angle)
    : points (points), neighborhood_size (neighborhood_size)
    , sharpness_angle (sharpness_angle), result (new double(0)) { }

  void operator()()
  {
    *result = CGAL::bilateral_smooth_point_set<Concurrency_tag>
      (points->all_or_selection_if_not_empty(),
       neighborhood_size,
       points->parameters().
       sharpness_angle(sharpness_angle).
       callback (*(this->callback())));
  }
};

using namespace CGAL::Three;
class Polyhedron_demo_point_set_bilateral_smoothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionBilateralSmoothing;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionBilateralSmoothing = new QAction(tr("Bilateral Smoothing"), mainWindow);
    actionBilateralSmoothing->setProperty("subMenuName","Point Set Processing");
    actionBilateralSmoothing->setObjectName("actionBilateralSmoothing");
    autoConnectActions();
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionBilateralSmoothing;
  }

public Q_SLOTS:
  void on_actionBilateralSmoothing_triggered();

}; // end Polyhedron_demo_point_set_bilateral_smoothing_plugin

class Point_set_demo_point_set_bilateral_smoothing_dialog : public QDialog, private Ui::PointSetBilateralSmoothingDialog
{
  Q_OBJECT
  public:
    Point_set_demo_point_set_bilateral_smoothing_dialog(QWidget * /*parent*/ = 0)
    {
      setupUi(this);
    }

    unsigned int iterations() const { return m_iterations->value(); }
    unsigned int neighborhood_size () const { return m_neighborhoodSize->value(); }
    unsigned int sharpness_angle () const { return m_sharpnessAngle->value(); }
};

void Polyhedron_demo_point_set_bilateral_smoothing_plugin::on_actionBilateralSmoothing_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    if(!item->has_normals())
    {
      QMessageBox::warning(mw,
                           tr("Cannot smooth"),
                           tr("%1 does not have normals.")
                           .arg(item->name()));
      return;
    }
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    Point_set_demo_point_set_bilateral_smoothing_dialog dialog;
    if(!dialog.exec())
      return;

    std::cerr << "Bilateral smoothing using "
              << dialog.iterations () << " iteration(s), neighborhood size of "
              << dialog.neighborhood_size () << " and sharpness angle of "
              << dialog.sharpness_angle () << "... ";
    QApplication::setOverrideCursor(Qt::BusyCursor);

    CGAL::Timer task_timer; task_timer.start();

    for (unsigned int i = 0; i < dialog.iterations (); ++i)
      {
        std::ostringstream oss;
        oss << "Bilateral smoothing (iteration " << i+1 << "/" << dialog.iterations() << ")";

        Bilateral_smoothing_functor functor (points,
                                             dialog.neighborhood_size (),
                                             dialog.sharpness_angle ());
        run_with_qprogressdialog (functor, oss.str().c_str(), mw);
        double error = *functor.result;

        if (std::isnan(error)) // NaN return means algorithm was interrupted
          break;
      }

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << task_timer.time() << " seconds, "
              << (memory>>20) << " Mb allocated)"
              << std::endl;

    // Updates scene
    item->invalidateOpenGLBuffers();
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();

  }
}

#include "Point_set_bilateral_smoothing_plugin.moc"
