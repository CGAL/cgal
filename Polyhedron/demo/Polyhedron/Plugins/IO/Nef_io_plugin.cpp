#include "Scene_nef_polyhedron_item.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <fstream>
#include <limits>

using namespace CGAL::Three;
class Polyhedron_demo_io_nef_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "nef_io_plugin.json")

public:
  QString nameFilters() const override;
  QString name() const override { return "io_nef_plugin"; }
  bool canLoad(QFileInfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override;

  bool canSave(const CGAL::Three::Scene_item*) override;
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items) override;
  bool isDefaultLoader(const Scene_item* item) const override{
    if(qobject_cast<const Scene_nef_polyhedron_item*>(item))
      return true;
    return false;
  }
};

QString Polyhedron_demo_io_nef_plugin::nameFilters() const {
  return "nef files (*.nef3)";
}

bool Polyhedron_demo_io_nef_plugin::canLoad(QFileInfo) const {
  return true;
}


QList<Scene_item*> Polyhedron_demo_io_nef_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {
  //do not try file with extension different from nef3
  if (fileinfo.suffix() != "nef3")
  {
    ok = false;
    return QList<Scene_item*>();
  }

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }

  // Try to read .nef3 in a polyhedron
  Scene_nef_polyhedron_item* item = new Scene_nef_polyhedron_item();
  item->setName(fileinfo.baseName());
  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }
  if(!item->load(in))
  {
    delete item;
    ok = false;
    return QList<Scene_item*>();
  }

  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>()<<item;
}

bool Polyhedron_demo_io_nef_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports polyhedrons and polygon soups
  return qobject_cast<const Scene_nef_polyhedron_item*>(item);
}

bool Polyhedron_demo_io_nef_plugin::save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // This plugin supports polyhedrons and polygon soups
  const Scene_nef_polyhedron_item* nef_item =
    qobject_cast<const Scene_nef_polyhedron_item*>(item);

  if(!nef_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());
  out.precision (std::numeric_limits<double>::digits10 + 2);
  items.pop_front();
  return (nef_item && nef_item->save(out));
}

#include "Nef_io_plugin.moc"
