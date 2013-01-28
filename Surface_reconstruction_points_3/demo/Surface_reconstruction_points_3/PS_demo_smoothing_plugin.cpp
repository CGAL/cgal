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

#include <CGAL/jet_smooth_point_set.h>

class PS_demo_smoothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionJetSmoothing;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionJetSmoothing = this->getActionFromMainWindow(mw, "actionJetSmoothing");
    if(actionJetSmoothing) {
      connect(actionJetSmoothing, SIGNAL(triggered()),
              this, SLOT(on_actionJetSmoothing_triggered()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionJetSmoothing;
  }

public slots:
  void on_actionJetSmoothing_triggered();

}; // end PS_demo_smoothing_plugin

void PS_demo_smoothing_plugin::on_actionJetSmoothing_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

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

Q_EXPORT_PLUGIN2(PS_demo_smoothing_plugin, PS_demo_smoothing_plugin)

#include "PS_demo_smoothing_plugin.moc"
