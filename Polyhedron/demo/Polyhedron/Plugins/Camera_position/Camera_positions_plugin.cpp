#include <QtCore/qglobal.h>
#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include "Camera_positions_list.h"
#include <CGAL/Three/Viewer_interface.h>

#include <QMainWindow>
using namespace CGAL::Three;
class Polyhedron_demo_camera_positions_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* );
  QList<QAction*> actions() const;

  QString name() const { return "camera_positions_plugin"; }
  QString nameFilters() const { return "Camera positions (*.camera.txt)"; }
  bool canLoad() const { return true; }
  Scene_item* load(QFileInfo fileinfo) { cpl->load(fileinfo.filePath()); return 0; }

  bool canSave(const Scene_item*) { return false; }
  bool save(const Scene_item*, QFileInfo ) {return false; }
  bool applicable(QAction*) const {return false;}
private:
  Camera_positions_list* cpl;
};

void Polyhedron_demo_camera_positions_plugin::init(QMainWindow* mainWindow, Scene_interface*, Messages_interface *)
{

  cpl = new Camera_positions_list(mainWindow);
  Viewer_interface* viewer = mainWindow->findChild<Viewer_interface*>("viewer");
  cpl->setViewer(viewer);
  mainWindow->addDockWidget(Qt::LeftDockWidgetArea, cpl);
}

QList<QAction*> 
Polyhedron_demo_camera_positions_plugin::actions() const
{
  return QList<QAction*>();
}

#include "Camera_positions_plugin.moc"
