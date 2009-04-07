#include "Point_set_scene_item.h"
#include "Kernel_type.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

class Point_set_demo_xyz_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface);

public:
  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QStringList Point_set_demo_xyz_plugin::nameFilters() const {
  return QStringList() << "XYZ files (*.xyz)"
                       << "Point Sets with Normal (*.pwn)";
};

bool Point_set_demo_xyz_plugin::canLoad() const {
  return true;
}


Scene_item* 
Point_set_demo_xyz_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }

  // Try to read .xyz in a point set
  Point_set_scene_item* point_set_item = new Point_set_scene_item;
  point_set_item->setName(fileinfo.baseName());
  if(!point_set_item->read_xyz_point_set(in)) {
    delete point_set_item;
    return NULL;
  }
  return point_set_item;
}

bool Point_set_demo_xyz_plugin::canSave(const Scene_item* item)
{
  // This plugin supports point sets
  return qobject_cast<const Point_set_scene_item*>(item);
}

bool Point_set_demo_xyz_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  // This plugin supports point sets
  const Point_set_scene_item* point_set_item = 
    qobject_cast<const Point_set_scene_item*>(item);
  if(!point_set_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());

  return point_set_item->write_xyz_point_set(out);
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Point_set_demo_xyz_plugin, Point_set_demo_xyz_plugin);
#include "Point_set_demo_xyz_plugin.moc"
