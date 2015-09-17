#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include <CGAL/gocad_io.h>
#include <CGAL/Timer.h>
#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

#include <QColor>

 

class Polyhedron_demo_gocad_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  QString nameFilters() const;
  QString name() const { return "gocad_plugin"; }
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QString Polyhedron_demo_gocad_plugin::nameFilters() const {
  return "GOCAD files (*.ts *.xyz)";
}

bool Polyhedron_demo_gocad_plugin::canLoad() const {
  return true;
}


Scene_item* 
Polyhedron_demo_gocad_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
    
  // Try to read GOCAD file in a polyhedron

  CGAL::Timer t;
  t.start();
  Scene_polyhedron_item* item = new Scene_polyhedron_item(Polyhedron());
  Polyhedron& P = * const_cast<Polyhedron*>(item->polyhedron());
  
  std::string name, color;
  if(! read_gocad(P, in, name, color)){
    // std::cerr << "Error: Invalid polyhedron" << std::endl;
    delete item;
    return 0;
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
    item->invalidate_buffers();
  }
  

  return item;
}

bool Polyhedron_demo_gocad_plugin::canSave(const Scene_item* item)
{
  // This plugin supports polyhedrons
  return qobject_cast<const Scene_polyhedron_item*>(item);
}

bool Polyhedron_demo_gocad_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  // This plugin supports polyhedrons
  const Scene_polyhedron_item* poly_item = 
    qobject_cast<const Scene_polyhedron_item*>(item);
 
  if(!poly_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());

  Polyhedron* poly = const_cast<Polyhedron*>(poly_item->polyhedron());

  write_gocad(*poly, out, qPrintable(fileinfo.baseName()));

  
  return true;

}

#include "Polyhedron_demo_gocad_plugin.moc"
