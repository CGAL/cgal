#include <QtCore/qglobal.h>
#include "Messages_interface.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_io_plugin_interface.h"

#include "Camera_positions_list.h"
#include "Viewer_interface.h"

#include <QMainWindow>

class Polyhedron_demo_camera_positions_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface Polyhedron_demo_io_plugin_interface)
public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const;

  QString name() const { return "camera_positions_plugin"; }
  QString nameFilters() const { return "Camera positions (*.camera.txt)"; }
  bool canLoad() const { return true; }
  Scene_item* load(QFileInfo fileinfo) { cpl->load(fileinfo.filePath()); return 0; }

  bool canSave(const Scene_item*) { return false; }
  bool save(const Scene_item*, QFileInfo ) {return false; }
  bool applicable() const {return false;}
private:
  Camera_positions_list* cpl;
};

void Polyhedron_demo_camera_positions_plugin::init(QMainWindow* mainWindow, Scene_interface*)
{
  cpl = new Camera_positions_list(mainWindow);
  cpl->setViewer(mainWindow->findChild<Viewer_interface*>("viewer"));
  mainWindow->addDockWidget(Qt::LeftDockWidgetArea, cpl);
}

QList<QAction*> 
Polyhedron_demo_camera_positions_plugin::actions() const
{
  return QList<QAction*>();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_camera_positions_plugin, Polyhedron_demo_camera_positions_plugin)

#include "Polyhedron_demo_camera_positions_plugin.moc"
