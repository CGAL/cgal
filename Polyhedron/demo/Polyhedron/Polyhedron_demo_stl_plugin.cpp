#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>

#include <CGAL/IO/Polyhedron_builder_from_STL.h>

#include <QColor>

class Polyhedron_demo_stl_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QString nameFilters() const;
  QString name() const { return "stl_plugin"; }
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

QString Polyhedron_demo_stl_plugin::nameFilters() const {
  return "STL files (*.stl)";
}

bool Polyhedron_demo_stl_plugin::canLoad() const {
  return true;
}


Scene_item* 
Polyhedron_demo_stl_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
    
  // Try to read STL file in a polyhedron
  Polyhedron P;
  //Scene_polyhedron_item* item = new Scene_polyhedron_item(P);
  //item->setName(fileinfo.baseName());
  
  CGAL::Polyhedron_builder_from_STL<Polyhedron::HalfedgeDS> builder(in);
  //  item->polyhedron()->

  P.delegate(builder);
 
  if(! P.is_valid()){
    std::cerr << "Error: Invalid polyhedron" << std::endl;
    return 0;
  }  

  if(P.empty()){
    return 0;
  }  

  Scene_polyhedron_item* item = new Scene_polyhedron_item(P);
  if(builder.name.size() == 0){
    item->setName(fileinfo.baseName());
  } else {
    item->setName(builder.name.c_str());
  }
  QColor color(builder.color.c_str());
  if(color.isValid()) 
  {  
    item->setColor(color);
    item->changed();
  }
  

  return item;
}

bool Polyhedron_demo_stl_plugin::canSave(const Scene_item*)
{
  return false;
}

bool Polyhedron_demo_stl_plugin::save(const Scene_item*, QFileInfo)
{
  return false;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Polyhedron_demo_stl_plugin, Polyhedron_demo_stl_plugin)
#include "Polyhedron_demo_stl_plugin.moc"
