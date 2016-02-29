#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include <CGAL/IO/Polyhedron_iostream.h>

#include <fstream>
using namespace CGAL::Three;
class Polyhedron_demo_off_to_xyz_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  QString name() const { return "off_to_xyz_plugin"; }
  QString nameFilters() const { return "OFF files as Point set (*.off)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
};

bool Polyhedron_demo_off_to_xyz_plugin::canLoad() const {
  return true;
}

CGAL::Three::Scene_item*
Polyhedron_demo_off_to_xyz_plugin::load(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8());

  if(!in)
    std::cerr << "Error!\n";

  Scene_points_with_normal_item* item;

  Polyhedron p;
  in >> p;
  if (in && !p.empty())
    item = new Scene_points_with_normal_item(p);
  else{
    in.close();
    in.open(fileinfo.filePath().toUtf8());
    item = new Scene_points_with_normal_item();
    if(!item->read_off_point_set(in))
    {
      delete item;
      return 0;
    }
  }

  item->setName(fileinfo.completeBaseName());
  return item;
}

bool Polyhedron_demo_off_to_xyz_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports point sets
  return qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool Polyhedron_demo_off_to_xyz_plugin::save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "off" && extension != "OFF")
    return false;

  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  // Save point set as .xyz
  std::ofstream out(fileinfo.filePath().toUtf8().data());
  out.precision (std::numeric_limits<double>::digits10 + 2);
  return point_set_item->write_off_point_set(out);
}


#include "OFF_to_xyz_io_plugin.moc"
