#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <QInputDialog>
#include <fstream>

class Polyhedron_demo_ply_to_xyz_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  QString name() const { return "ply_to_xyz_plugin"; }
  QString nameFilters() const { return "PLY files as Point set (*.ply)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
};

bool Polyhedron_demo_ply_to_xyz_plugin::canLoad() const {
  return true;
}

CGAL::Three::Scene_item*
Polyhedron_demo_ply_to_xyz_plugin::load(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  Scene_points_with_normal_item* item;
  item = new Scene_points_with_normal_item();
  if(!item->read_ply_point_set(in))
    {
      delete item;
      return 0;
    }

  item->setName(fileinfo.completeBaseName());
  return item;
}

bool Polyhedron_demo_ply_to_xyz_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports point sets
  return qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool Polyhedron_demo_ply_to_xyz_plugin::save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "ply" && extension != "PLY")
    return false;

  QStringList list;
  list << tr("Binary");
  list << tr("Ascii");
  bool ok = false;
  QString choice
    = QInputDialog::getItem(NULL, tr("Save PLY file"), tr("Format"), list, 0, false, &ok);

  if (!ok)
    return false;
  
  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  // Save point set as .xyz
  std::ofstream out(fileinfo.filePath().toUtf8().data());
  out.precision (std::numeric_limits<double>::digits10 + 2);
  return point_set_item->write_ply_point_set(out, (choice == tr("Binary")));
}


#include "PLY_to_xyz_io_plugin.moc"
