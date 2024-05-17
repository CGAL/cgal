#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Point_set_3/IO.h>
#include <fstream>
using namespace CGAL::Three;
class CGAL_Lab_las_plugin :
  public QObject,
  public CGAL::Three::CGAL_Lab_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90" FILE "las_io_plugin.json")

public:
  QString name() const { return "las_plugin"; }
  QString nameFilters() const { return "LAS files (*.las);;Compressed LAS files (*.laz)"; }
  bool canLoad(QFileInfo fileinfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& );
};

bool CGAL_Lab_las_plugin::canLoad(QFileInfo ) const {
  return true;
}

QList<Scene_item*> CGAL_Lab_las_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  Scene_points_with_normal_item* item;
  item = new Scene_points_with_normal_item();
  Q_ASSERT(item->point_set() != nullptr);

  item->point_set()->clear();

  ok = CGAL::IO::read_LAS (in, *(item->point_set())) && !item->isEmpty();
  if(ok)
  {
    std::cerr << item->point_set()->info();

    if (!item->point_set()->has_normal_map())
    {
      item->setRenderingMode(Points);
    }
    else{
      item->setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
    }
    if (item->point_set()->check_colors())
      std::cerr << "-> Point set has colors" << std::endl;

    item->invalidateOpenGLBuffers();
  }
  if(!ok)
  {
    delete item;
    return QList<Scene_item*>();
  }

  item->setName(fileinfo.completeBaseName());
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>()<<item;
}

bool CGAL_Lab_las_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool CGAL_Lab_las_plugin::save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "las" && extension != "LAS")
    return false;

  // This plugin supports point sets
  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8().data());
  items.pop_front();
  Q_ASSERT(point_set_item->point_set() != nullptr);

  point_set_item->point_set()->reset_indices();

  return CGAL::IO::write_LAS(out, *(point_set_item->point_set()));
}


#include "LAS_io_plugin.moc"
