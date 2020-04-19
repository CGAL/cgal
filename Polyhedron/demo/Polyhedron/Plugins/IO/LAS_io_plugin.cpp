#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Scene_item.h>
#include <fstream>
using namespace CGAL::Three;
class Polyhedron_demo_las_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "las_io_plugin.json")

public:
  QString name() const { return "las_plugin"; }
  QString nameFilters() const { return "LAS files (*.las);;Compressed LAS files (*.laz)"; }
  bool canLoad(QFileInfo fileinfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& );
};

bool Polyhedron_demo_las_plugin::canLoad(QFileInfo ) const {
  return true;
}

QList<Scene_item*> Polyhedron_demo_las_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  Scene_points_with_normal_item* item;
  item = new Scene_points_with_normal_item();
  if(!item->read_las_point_set(in))
    {
      delete item;
      ok = false;
      return QList<Scene_item*>();
    }

  item->setName(fileinfo.completeBaseName());
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>()<<item;
}

bool Polyhedron_demo_las_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool Polyhedron_demo_las_plugin::save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "las" && extension != "LAS")
    return false;

  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8().data());
  items.pop_front();
  return point_set_item->write_las_point_set(out);
}


#include "LAS_io_plugin.moc"
