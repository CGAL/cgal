#include "SMesh_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Scene.h"

#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include <CGAL/IO/OM.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <QColor>
#include <QString>
#include <QStringList>
#include <QMainWindow>
#include <QInputDialog>

#include <boost/property_map/function_property_map.hpp>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

using namespace CGAL::Three;
class CGAL_Lab_om_plugin :
  public QObject,
  public CGAL_Lab_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90" FILE "om_io_plugin.json")

public:
  QString nameFilters() const;
  QString name() const { return "om_plugin"; }
  bool canLoad(QFileInfo fileinfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&);
};

QString CGAL_Lab_om_plugin::nameFilters() const {
  return "om files (*.om)";
}

bool CGAL_Lab_om_plugin::canLoad(QFileInfo) const {
  return true;
}



QList<Scene_item*>
CGAL_Lab_om_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene){

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios::in | std::ios::binary);
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }
  in.close();
  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_surface_mesh_item* item = new Scene_surface_mesh_item();
    item->setName(fileinfo.completeBaseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }

  try
  {
    typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;
    std::map<vertex_descriptor,bool> sm_selected_map;
    auto sm_selected_pmap = boost::make_assoc_property_map(sm_selected_map);

    std::map<edge_descriptor,bool> sm_feature_map;
    auto sm_feature_pmap = boost::make_assoc_property_map(sm_feature_map);

    // Try building a surface_mesh
    SMesh* SM = new SMesh();
    if (CGAL::IO::read_OM((const char*)fileinfo.filePath().toUtf8(), *SM, sm_selected_pmap, sm_feature_pmap))
    {/*
      std::cout << "vertex selection values:\n";
      for(auto v : vertices(*SM)){
        std::cout << std::boolalpha << get(sm_selected_pmap, v) << std::endl;
      }

      std::cout << "edge feature values:\n";
      for(auto e : edges(*SM)){
        std::cout  << std::boolalpha << get(sm_feature_pmap, e) << std::endl;
      }
      */
    }

    if(!SM->is_valid() || SM->is_empty()){
      std::cerr << "Error: Invalid facegraph" << std::endl;
    }
    else{
      Scene_surface_mesh_item* item = new Scene_surface_mesh_item(SM);
      item->setName(fileinfo.completeBaseName());
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>()<<item;
    }
  }
  catch(...){}

  ok = false;
  return QList<Scene_item*>();
}

bool CGAL_Lab_om_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_surface_mesh_item*>(item);
}

bool CGAL_Lab_om_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
#if 0
  Scene_item* item = items.front();
  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(item);

  if(!sm_item)
    return false;

  QStringList list;
  list << tr("Binary");
  list << tr("Ascii");
  bool ok = false;
  QString choice
    = QInputDialog::getItem(nullptr, tr("Save om file"), tr("Format"), list, 0, false, &ok);

  if (!ok)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8(), std::ios::out | std::ios::binary);
  if ( choice == tr("Binary") )
    CGAL::IO::set_mode(out, CGAL::IO::BINARY);
  else
  {
    CGAL::IO::set_mode(out, CGAL::IO::ASCII);
    out.precision (std::numeric_limits<double>::digits10 + 2);
  }

  if (sm_item)
  {
    CGAL::IO::write_om(out, *sm_item->face_graph());
    items.pop_front();
    return true;
  }
  #endif
  return false;
}

#include "om_io_plugin.moc"