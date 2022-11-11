#include <QtCore/qglobal.h>
#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include "Camera_positions_list.h"

#include <CGAL/Three/Three.h>
#include <CGAL/Three/Viewer_interface.h>

#include <QMainWindow>
using namespace CGAL::Three;
class Polyhedron_demo_camera_positions_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90")

public:
  void init() override;

  QString name() const override { return "camera_positions_plugin"; }
  QString nameFilters() const override { return "Camera positions (*.camera.txt)"; }
  bool canLoad(QFileInfo) const override { return true; }
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override
  {
    Q_UNUSED(add_to_scene);
    ok = true;
    cpl->load(fileinfo.filePath());
    return QList<Scene_item*>();
  }

  bool canSave(const Scene_item*) override { return true; }
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& ) override
  {
    return cpl->save(fileinfo.filePath());
  }
private:
  Camera_positions_list* cpl;
};

void Polyhedron_demo_camera_positions_plugin::init()
{

  cpl = new Camera_positions_list(CGAL::Three::Three::mainWindow());
  CGAL::Three::Three::mainWindow()->addDockWidget(Qt::LeftDockWidgetArea, cpl);
}

#include "Camera_positions_plugin.moc"
