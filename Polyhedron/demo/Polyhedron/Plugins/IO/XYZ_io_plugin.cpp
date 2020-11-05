#include "Scene_points_with_normal_item.h"
#include "Kernel_type.h"

#include <QMainWindow>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <fstream>
#include <QMessageBox>
#include <QMenu>
using namespace CGAL::Three;

class Polyhedron_demo_xyz_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "xyz_io_plugin.json")

public:

  QString name() const { return "xyz_plugin"; }
  QString nameFilters() const { return "XYZ as Point Set (*.xyz);;Point Set with Normal (*.pwn)"; }
  bool canLoad(QFileInfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&);
};

bool Polyhedron_demo_xyz_plugin::canLoad(QFileInfo) const {
  return true;
}


QList<Scene_item*>
Polyhedron_demo_xyz_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene)
{
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8().data());
  if(!in) {
    std::cerr << "Error! Cannot open file " << fileinfo.filePath().toStdString() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }


  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_points_with_normal_item* item =
        new Scene_points_with_normal_item();
    item->setName(fileinfo.completeBaseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }
  // Read .xyz in a point set
  Scene_points_with_normal_item* point_set_item = new Scene_points_with_normal_item;
  point_set_item->setName(fileinfo.completeBaseName());
  if(!point_set_item->read_xyz_point_set(in)) {
    delete point_set_item;
    ok = false;
    return QList<Scene_item*>();
  }
  if(point_set_item->has_normals())
    point_set_item->setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());

  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(point_set_item);
  return QList<Scene_item*>()<<point_set_item;
}

bool Polyhedron_demo_xyz_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports point sets
  return qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool Polyhedron_demo_xyz_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "xyz" && extension != "XYZ" &&
      extension != "pwn" && extension != "PWN")
    return false;

  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  // Save point set as .xyz
  std::ofstream out(fileinfo.filePath().toUtf8().data());
  out.precision (std::numeric_limits<double>::digits10 + 2);
  bool res = point_set_item->write_xyz_point_set(out);
  if(res)
    items.pop_front();
  return res;
}

#include "XYZ_io_plugin.moc"
