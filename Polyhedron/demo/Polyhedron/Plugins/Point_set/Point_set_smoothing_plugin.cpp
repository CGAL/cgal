#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include <CGAL/jet_smooth_point_set.h>

#include "run_with_qprogressdialog.h"

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

struct Jet_smoothing_functor
  : public Functor_with_signal_callback
{
  Point_set* points;
  const int nb_neighbors;

  Jet_smoothing_functor (Point_set* points, const int nb_neighbors)
    : points (points), nb_neighbors (nb_neighbors) { }

  void operator()()
  {
    CGAL::jet_smooth_point_set<Concurrency_tag>(points->all_or_selection_if_not_empty(),
                                                nb_neighbors,
                                                points->parameters().
                                                callback (*(this->callback())));
  }
};

using namespace CGAL::Three;
class Polyhedron_demo_point_set_smoothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionJetSmoothing;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionJetSmoothing = new QAction(tr("Jet Smoothing"), mainWindow);
    actionJetSmoothing->setProperty("subMenuName","Point Set Processing");
    actionJetSmoothing->setObjectName("actionJetSmoothing");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionJetSmoothing;
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionJetSmoothing_triggered();

}; // end Polyhedron_demo_point_set_smoothing_plugin

void Polyhedron_demo_point_set_smoothing_plugin::on_actionJetSmoothing_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    Point_set* points = item->point_set();
    if(!points) return;

    // Gets options
    bool ok;

    const unsigned int nb_neighbors =
      QInputDialog::getInt((QWidget*)mw,
                           tr("Jet Smoothing"), // dialog title
                           tr("Number of neighbors:"), // field label
                           24, // default value = fast
                           6, // min
                           1000, // max
                           1, // step
                           &ok);
    if(!ok) return;

    QApplication::setOverrideCursor(Qt::BusyCursor);
    QApplication::processEvents();

    Jet_smoothing_functor functor (points, nb_neighbors);
    run_with_qprogressdialog (functor, "Smoothing point set...", mw);

    points->invalidate_bounds();

    // update scene
    item->invalidateOpenGLBuffers();
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();
  }
}

#include "Point_set_smoothing_plugin.moc"
