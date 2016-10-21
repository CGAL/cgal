#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include "read_las_point_set.h"

#include <fstream>

class Polyhedron_demo_las_to_xyz_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  QString name() const { return "las_to_xyz_plugin"; }
  QString nameFilters() const { return "LAS files as Point set (*.las);;Compressed LAZ files as Point set (*.laz)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
};

bool Polyhedron_demo_las_to_xyz_plugin::canLoad() const {
  return true;
}

CGAL::Three::Scene_item*
Polyhedron_demo_las_to_xyz_plugin::load(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  Scene_points_with_normal_item* item;
  item = new Scene_points_with_normal_item();
  if(!read_las_point_set(in, *(item->point_set())) || item->isEmpty())
    {
      delete item;
      return 0;
    }
  item->point_set()->check_colors();
  std::cerr << item->point_set()->info();
  
  item->setName(fileinfo.completeBaseName());
  return item;
}

bool Polyhedron_demo_las_to_xyz_plugin::canSave(const CGAL::Three::Scene_item*)
{
  // No save function yet
  return false;
}

bool Polyhedron_demo_las_to_xyz_plugin::save(const CGAL::Three::Scene_item*, QFileInfo)
{
  return false;
}


#include "LAS_to_xyz_io_plugin.moc"
