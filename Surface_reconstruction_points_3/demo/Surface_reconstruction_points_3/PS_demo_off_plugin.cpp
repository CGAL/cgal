#include "Point_set_scene_item.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

class PS_demo_off_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QStringList PS_demo_off_plugin::nameFilters() const {
  return QStringList() << "OFF files (*.off)";
}

bool PS_demo_off_plugin::canLoad() const {
  return true;
}


Scene_item*
PS_demo_off_plugin::load(QFileInfo fileinfo) {

  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "off" && extension != "OFF")
    return false;

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8().data());
  if(!in) {
    std::cerr << "Error! Cannot open file " << fileinfo.filePath().toStdString() << std::endl;
    return NULL;
  }

  // Try to read .off in a polyhedron
  Scene_polyhedron_item* item = new Scene_polyhedron_item();
  item->setName(fileinfo.completeBaseName());
  if(!item->load(in))
  {
    delete item;

    // Try to read .off in a point set
    Point_set_scene_item* point_set_item = new Point_set_scene_item;
    point_set_item->setName(fileinfo.completeBaseName());
    in.close();
    std::ifstream in2(fileinfo.filePath().toUtf8().data());
    if(!point_set_item->read_off_point_set(in2)) {
      delete point_set_item;
      return 0;
    }
    return point_set_item;
  }

  return item;
}

bool PS_demo_off_plugin::canSave(const Scene_item* item)
{
  // This plugin supports polyhedrons and point sets
  return qobject_cast<const Scene_polyhedron_item*>(item) ||
         qobject_cast<const Point_set_scene_item*>(item);
}

bool PS_demo_off_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "off" && extension != "OFF")
    return false;

  // This plugin supports polyhedrons and point sets
  const Scene_polyhedron_item* poly_item =
    qobject_cast<const Scene_polyhedron_item*>(item);
  const Point_set_scene_item* point_set_item =
    qobject_cast<const Point_set_scene_item*>(item);
  if(!poly_item && !point_set_item)
    return false;

  // Save polyhedron/point set as .off
  std::ofstream out(fileinfo.filePath().toUtf8().data());
  return (poly_item && poly_item->save(out)) ||
         (point_set_item && point_set_item->write_off_point_set(out));
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(PS_demo_off_plugin, PS_demo_off_plugin)
#include "PS_demo_off_plugin.moc"
