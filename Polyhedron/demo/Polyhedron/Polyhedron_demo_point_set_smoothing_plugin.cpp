#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include <CGAL/jet_smooth_point_set.h>

class Polyhedron_demo_point_set_smoothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionJetSmoothing;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionJetSmoothing = new QAction(tr("Point set jet smoothing"), mainWindow);
    actionJetSmoothing->setObjectName("actionJetSmoothing");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionJetSmoothing;
  }

  bool applicable() const { 
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public slots:
  void on_actionJetSmoothing_triggered();

}; // end Polyhedron_demo_point_set_smoothing_plugin

void Polyhedron_demo_point_set_smoothing_plugin::on_actionJetSmoothing_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    Point_set* points = item->point_set();
    if(!points) return;

    // Gets options
    bool ok;
    const unsigned int nb_neighbors =
      QInputDialog::getInteger((QWidget*)mw,
                               tr("Jet Smoothing"), // dialog title
                               tr("Number of neighbors:"), // field label
                               24, // default value = fast
                               6, // min
                               1000, // max
                               1, // step
                               &ok);
    if(!ok) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    CGAL::jet_smooth_point_set(points->begin(), points->end(), nb_neighbors);

    points->invalidate_bounds();

    // update scene
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_point_set_smoothing_plugin, Polyhedron_demo_point_set_smoothing_plugin)

#include "Polyhedron_demo_point_set_smoothing_plugin.moc"
