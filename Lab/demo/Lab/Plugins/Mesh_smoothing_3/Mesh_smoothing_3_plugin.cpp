
#include <QtCore/qglobal.h>

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>

#include "Scene_c3t3_item.h"
#include "C3t3_type.h"

#include <unordered_map>
#include <memory>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

using namespace CGAL::Three;
class CGAL_Lab_mesh_smoothing_3_plugin :
  public QObject,
  public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "mesh_smoothing_3_plugin.json")

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionMeshSmoothing_ = new QAction("Mesh smoothing", mw);
    if (actionMeshSmoothing_) {
      connect(actionMeshSmoothing_, SIGNAL(triggered()),
        this, SLOT(mesh_smoothing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMeshSmoothing_;
  }

  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_c3t3_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void mesh_smoothing()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_c3t3_item* c3t3_item =
      qobject_cast<Scene_c3t3_item*>(scene->item(index));
    const auto& c3t3 = c3t3_item->c3t3();

    if (c3t3_item)
    {
      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QElapsedTimer time;
      time.start();

      //todo : add smoothing code here

      std::cout << "Smoothing done (" << time.elapsed() << " ms)" << std::endl;

      c3t3_item->invalidateOpenGLBuffers();
      c3t3_item->c3t3_changed();
      this->scene->itemChanged(index);

      // default cursor
      QApplication::restoreOverrideCursor();
    }
    else
    {
      std::cout << "Can't smooth that type of thing" << std::endl;
    }
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;


private:
  QAction* actionMeshSmoothing_;

}; // end CGAL_Lab_mesh_smoothing_3_plugin

#include "Mesh_smoothing_3_plugin.moc"
