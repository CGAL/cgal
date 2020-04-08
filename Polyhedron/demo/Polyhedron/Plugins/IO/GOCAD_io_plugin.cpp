#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Kernel_type.h"
#include "Scene.h"
#include "SMesh_type.h"
#include <CGAL/gocad_io.h>
#include <CGAL/Timer.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <fstream>

#include <QColor>
#include <QMainWindow>
using namespace CGAL::Three;

class Polyhedron_demo_gocad_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface

{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "gocad_io_plugin.json")

public:

  QString nameFilters() const;
  QString name() const { return "gocad_plugin"; }
  bool canLoad(QFileInfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo, QList<CGAL::Three::Scene_item*>& );
};

QString Polyhedron_demo_gocad_plugin::nameFilters() const {
  return "GOCAD files (*.ts)";
}

bool Polyhedron_demo_gocad_plugin::canLoad(QFileInfo) const {
  return true;
}


QList<Scene_item*>
Polyhedron_demo_gocad_plugin::load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }


  CGAL::Timer t;
  t.start();
  // Try to read GOCAD file in a surface_mesh
  Scene_surface_mesh_item* item = new Scene_surface_mesh_item(new SMesh());
  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    ok = true;
    if(add_to_scene)
       CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }
  SMesh& P = * const_cast<SMesh*>(item->polyhedron());

  std::string name, color;
  if(! read_gocad(P, in, name, color)){
    std::cerr << "Error: Invalid polyhedron" << std::endl;
    delete item;
    ok = false;
    return QList<Scene_item*>();
  }

  t.stop();
  std::cerr << "Reading took " << t.time() << " sec." << std::endl;
  if(name.size() == 0){
    item->setName(fileinfo.baseName());
  } else {
    item->setName(name.c_str());
  }
  QColor qcolor(color.c_str());
  if(qcolor.isValid())
  {
    item->setColor(qcolor);
  }
  item->invalidateOpenGLBuffers();
  ok = true;
  if(add_to_scene)
     CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>()<<item;
}

bool Polyhedron_demo_gocad_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports polyhedrons
  return  qobject_cast<const Scene_surface_mesh_item*>(item);
}

bool Polyhedron_demo_gocad_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // This plugin supports polyhedrons
  const Scene_surface_mesh_item* sm_item =
      qobject_cast<const Scene_surface_mesh_item*>(item);

  if(!sm_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());
  out.precision (std::numeric_limits<double>::digits10 + 2);
  SMesh* poly = const_cast<SMesh*>(sm_item->polyhedron());
  write_gocad(*poly, out, qPrintable(fileinfo.baseName()));
  items.pop_front();
  return true;

}

#include "GOCAD_io_plugin.moc"
