#include "Scene_nef_polyhedron_item.h"
#include "Nef_type.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <fstream>
using namespace CGAL::Three;
class Polyhedron_demo_off_to_nef_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0" FILE "off_to_nef_io_plugin.json")

public:
  QString name() const { return "off_to_nef_plugin"; }
  QString nameFilters() const { return "OFF files, into nef (*.off)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
};

bool Polyhedron_demo_off_to_nef_plugin::canLoad() const {
  return true;
}

CGAL::Three::Scene_item*
Polyhedron_demo_off_to_nef_plugin::load(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8());

  if(!in)
    std::cerr << "Error!\n";
  Scene_nef_polyhedron_item* item = new Scene_nef_polyhedron_item();
  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    item->setName(fileinfo.completeBaseName());
    return item;
  }
  if(!item->load_from_off(in))
  {
    delete item;
    return 0;
  }

  item->setName(fileinfo.baseName());
  return item;
}

bool Polyhedron_demo_off_to_nef_plugin::canSave(const CGAL::Three::Scene_item*)
{
  return false;
}

bool Polyhedron_demo_off_to_nef_plugin::save(const CGAL::Three::Scene_item*, QFileInfo)
{
  return false;
}

#include "OFF_to_nef_io_plugin.moc"
