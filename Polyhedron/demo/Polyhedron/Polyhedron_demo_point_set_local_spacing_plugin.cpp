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

class Polyhedron_demo_local_spacing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionRadiusFromDensity = new QAction(tr("Point set local spacing"), mainWindow);
    actionRadiusFromDensity->setObjectName("actionRadiusFromDensity");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRadiusFromDensity;
  }

  //! Applicable if the currently selected item is a
  //! points_with_normal_item.
  bool applicable() const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }
public slots:
  void on_actionRadiusFromDensity_triggered();

private:
  QAction* actionRadiusFromDensity;
}; // end Polyhedron_demo_local_spacing_plugin

void Polyhedron_demo_local_spacing_plugin::on_actionRadiusFromDensity_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Check there is a point set
    if(item->point_set() == NULL)
        return;

    // Gets options
    bool ok;
    const int k =
      QInputDialog::getInteger((QWidget*)mw,
                              tr("Local spacing"), // dialog title
                              tr("Number of neighbors:"), // field label
                              6, // default value = small
                              1, // min
                              1000, // max
                              1, // step
                              &ok);
    if(!ok) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    item->computes_local_spacing(k);

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_local_spacing_plugin, Polyhedron_demo_local_spacing_plugin)

#include "Polyhedron_demo_point_set_local_spacing_plugin.moc"
