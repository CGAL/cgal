#include "SMesh_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Scene.h"

#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include "Scene_polyhedron_selection_item.h"

#include <CGAL/boost/graph/IO/OM.h>
#include <CGAL/boost/graph/io.h>

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
    std::map<vertex_descriptor,bool> sm_vfeature_map;
    auto sm_vfeature_pmap = boost::make_assoc_property_map(sm_vfeature_map);

    std::map<edge_descriptor,bool> sm_efeature_map;
    auto sm_efeature_pmap = boost::make_assoc_property_map(sm_efeature_map);

    // Try building a surface_mesh
    SMesh* sm = new SMesh();
    ok = CGAL::IO::read_OM((const char*)fileinfo.filePath().toUtf8(),
                           *sm,
                           CGAL::parameters::vertex_is_constrained_map(sm_vfeature_pmap)
                           .edge_is_constrained_map(sm_efeature_pmap));

    if(!ok || !sm->is_valid() || sm->is_empty()){
      std::cerr << "Error: Invalid facegraph" << std::endl;
    }
    else{
      Scene_surface_mesh_item* item = new Scene_surface_mesh_item(sm);

      Scene_polyhedron_selection_item* selection_item = new Scene_polyhedron_selection_item(item, CGAL::Three::Three::mainWindow());
      for(auto v : vertices(*sm)){
        if(get(sm_vfeature_pmap, v)){
          selection_item->selected_vertices.insert(v);
        }
      }
      for(auto e : edges(*sm)){
        if(get(sm_efeature_pmap, e)){
          selection_item->selected_edges.insert(e);
        }
      }
      item->setName(fileinfo.completeBaseName());
      ok = true;
      if(add_to_scene)
      {
        CGAL::Three::Three::scene()->addItem(item);
        if(!selection_item->isEmpty())
          CGAL::Three::Three::scene()->addItem(selection_item);
      }

      QList<Scene_item*> res;
      res << item;
      if(!selection_item->isEmpty())
        res << selection_item;
      else
        delete selection_item;
      return res;
    }
  }
  catch(...){}

  ok = false;
  return QList<Scene_item*>();
}

bool CGAL_Lab_om_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_surface_mesh_item*>(item)
      || qobject_cast<const Scene_polyhedron_selection_item*>(item);
}

bool CGAL_Lab_om_plugin::
save(QFileInfo fileinfo, QList<CGAL::Three::Scene_item*>& items)
{
  Scene_surface_mesh_item* sm_item = nullptr;
  Scene_polyhedron_selection_item* selection_item = nullptr;

  for (Scene_item* item : items)
  {
    if (sm_item == nullptr)
    {
      sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
      if (sm_item != nullptr) //surface_mesh_item found
        continue;
    }
    if (selection_item == nullptr)
    {
      selection_item = qobject_cast<Scene_polyhedron_selection_item*>(item);
    }
  }

  if (sm_item == nullptr && selection_item == nullptr)
    return false;

  if (selection_item != nullptr)
  {
    if (sm_item == nullptr)
      sm_item = selection_item->polyhedron_item();

    if (sm_item != selection_item->polyhedron_item())
    {
      std::cerr << "Warning! Selection is not associated to the surface_mesh. Ignoring selection." << std::endl;
      selection_item = nullptr;
    }
  }

  bool res = false;
  if (selection_item != nullptr)
  {
    res = CGAL::IO::write_OM((const char*)fileinfo.filePath().toUtf8()
              , *sm_item->face_graph()
              , CGAL::parameters::vertex_is_constrained_map(selection_item->constrained_vertices_pmap())
              .edge_is_constrained_map(selection_item->constrained_edges_pmap())
              .stream_precision(17));
  }
  else
  {
    res = CGAL::IO::write_OM((const char*)fileinfo.filePath().toUtf8()
              , *sm_item->face_graph()
              , CGAL::parameters::stream_precision(17));
  }

  if (res)
  {
    items.removeAll(sm_item);
    if (selection_item != nullptr)
      items.removeAll(selection_item);
  }
  return res;
}

#include "OM_io_plugin.moc"
