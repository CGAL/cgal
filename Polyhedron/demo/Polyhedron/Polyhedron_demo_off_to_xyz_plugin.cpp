#include "Scene_points_with_normal_item.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

class Polyhedron_demo_off_to_xyz_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QString name() const { return "off_to_xyz_plugin"; }
  QString nameFilters() const { return "OFF files as Point set (*.off)"; }
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

bool Polyhedron_demo_off_to_xyz_plugin::canLoad() const {
  return true;
}

Scene_item* 
Polyhedron_demo_off_to_xyz_plugin::load(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8());

  if(!in)
    std::cerr << "Error!\n";
  Scene_points_with_normal_item* item = new Scene_points_with_normal_item();

  if(!item->read_off_point_set(in))
  {
    delete item;
    return 0;
  }

  item->setName(fileinfo.baseName());
  return item;
}

bool Polyhedron_demo_off_to_xyz_plugin::canSave(const Scene_item*)
{
  return false;
}

bool Polyhedron_demo_off_to_xyz_plugin::save(const Scene_item*, QFileInfo)
{
  return false;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Polyhedron_demo_off_to_xyz_plugin, Polyhedron_demo_off_to_xyz_plugin)
#include "Polyhedron_demo_off_to_xyz_plugin.moc"
