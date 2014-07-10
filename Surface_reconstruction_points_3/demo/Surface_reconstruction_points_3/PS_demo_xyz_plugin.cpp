#include "Point_set_scene_item.h"
#include "Kernel_type.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

class PS_demo_xyz_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

  #if QT_VERSION >= 0x050000
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")//New for Qt5 version !
  #endif

public:
  QStringList nameFilters() const;
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QStringList PS_demo_xyz_plugin::nameFilters() const {
  return QStringList() << "XYZ files (*.xyz)"
                       << "Point Sets with Normal (*.pwn)";
}

bool PS_demo_xyz_plugin::canLoad() const {
  return true;
}


Scene_item*
PS_demo_xyz_plugin::load(QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "xyz" && extension != "XYZ" &&
      extension != "pwn" && extension != "PWN")
    return NULL;

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8().data());
  if(!in) {
    std::cerr << "Error! Cannot open file " << fileinfo.filePath().toStdString() << std::endl;
    return NULL;
  }

  // Read .xyz in a point set
  Point_set_scene_item* point_set_item = new Point_set_scene_item;
  point_set_item->setName(fileinfo.completeBaseName());
  if(!point_set_item->read_xyz_point_set(in)) {
    delete point_set_item;
    return NULL;
  }
  return point_set_item;
}

bool PS_demo_xyz_plugin::canSave(const Scene_item* item)
{
  // This plugin supports point sets
  return qobject_cast<const Point_set_scene_item*>(item);
}

bool PS_demo_xyz_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "xyz" && extension != "XYZ" &&
      extension != "pwn" && extension != "PWN")
    return false;

  // This plugin supports point sets
  const Point_set_scene_item* point_set_item =
    qobject_cast<const Point_set_scene_item*>(item);
  if(!point_set_item)
    return false;

  // Save point set as .xyz
  std::ofstream out(fileinfo.filePath().toUtf8().data());
  return point_set_item->write_xyz_point_set(out);
}

#if QT_VERSION < 0x050000
#include <QtPlugin>
Q_EXPORT_PLUGIN2(PS_demo_xyz_plugin, PS_demo_xyz_plugin)
#endif

#include "PS_demo_xyz_plugin.moc"
