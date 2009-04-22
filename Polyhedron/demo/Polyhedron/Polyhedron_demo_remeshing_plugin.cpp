#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include "Scene_polyhedron_item.h"
#include <QInputDialog>

// declare the CGAL function
Polyhedron* cgal_code_remesh(const Polyhedron*,
                             const double angle,
                             const double sizing,
                             const double approx);

class Polyhedron_demo_remeshing_plugin : 
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);
public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionRemeshing = this->getActionFromMainWindow(mw, "actionRemeshing");
    if(actionRemeshing) {
      connect(actionRemeshing, SIGNAL(triggered()),
              this, SLOT(remesh()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRemeshing;
  }
public slots:
  void remesh();

private:
  QAction* actionRemeshing;
}; // end class Polyhedron_demo_remeshing_plugin

void Polyhedron_demo_remeshing_plugin::remesh()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    // TODO: get parameters using ONE dialog box
    // sizing and approximation parameters should be expressed as ratio of 
    // scene bbox diagonal.

    double diag = scene->len_diagonal();

    bool ok;
    const double angle = 
      QInputDialog::getDouble(mw, tr("Min triangle angle"),
                              tr("Angle:"),
                              25., // default value
                              1., // min
                              30., // max
                              2, // decimals
                              &ok);
    if(!ok) return;

    const double sizing = 
      QInputDialog::getDouble(mw, "Sizing",
      "Size:",
      diag * 0.05, // default value
      diag * 10e-6, // min
      diag, // max
      4, // decimals
      &ok);
    if(!ok) return;

    const double approx = 
      QInputDialog::getDouble(mw, "Approximation error",
      "Error:",
      diag * 0.005, // default value
      diag * 10e-7, // min
      diag, // max
      6, // decimals
      &ok);
    if(!ok) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    Polyhedron* pRemesh = cgal_code_remesh(pMesh, angle, sizing, approx);

    if(pRemesh)
    {
      Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pRemesh);
      new_item->setName(tr("%1 remeshed (%2 %3 %4)")
                         .arg(item->name())
                         .arg(angle)
                         .arg(sizing)
                         .arg(approx));
      new_item->setColor(Qt::magenta);
      new_item->setRenderingMode(item->renderingMode());
      item->setVisible(false);
      scene->itemChanged(index);
      scene->addItem(new_item);
    }

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_remeshing_plugin, Polyhedron_demo_remeshing_plugin);

#include "Polyhedron_demo_remeshing_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
