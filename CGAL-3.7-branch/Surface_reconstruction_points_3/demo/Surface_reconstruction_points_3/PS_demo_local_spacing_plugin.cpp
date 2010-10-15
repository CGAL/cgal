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

class PS_demo_local_spacing_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionRadiusFromDensity = this->getActionFromMainWindow(mw, "actionRadiusFromDensity");
    if(actionRadiusFromDensity) {
      connect(actionRadiusFromDensity, SIGNAL(triggered()),
              this, SLOT(on_actionRadiusFromDensity_triggered()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRadiusFromDensity;
  }

public slots:
  void on_actionRadiusFromDensity_triggered();

private:
  QAction* actionRadiusFromDensity;
}; // end PS_demo_local_spacing_plugin

void PS_demo_local_spacing_plugin::on_actionRadiusFromDensity_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

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

Q_EXPORT_PLUGIN2(PS_demo_local_spacing_plugin, PS_demo_local_spacing_plugin)

#include "PS_demo_local_spacing_plugin.moc"
