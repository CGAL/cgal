#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Scene.h"
#include "SMesh_type.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/IO/STL_writer.h>
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
class Polyhedron_demo_stl_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "stl_io_plugin.json")

public:
  QString nameFilters() const;
  QString name() const { return "stl_plugin"; }
  bool canLoad(QFileInfo fileinfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&);
};

QString Polyhedron_demo_stl_plugin::nameFilters() const {
  return "STL files (*.stl)";
}

bool Polyhedron_demo_stl_plugin::canLoad(QFileInfo) const {
  return true;
}



QList<Scene_item*>
Polyhedron_demo_stl_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene){

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios::in | std::ios::binary);
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }
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
  std::vector<std::array<double, 3> > points;
  std::vector<std::array<int, 3> > triangles;
  if (!CGAL::read_STL(in, points, triangles))
  {
    std::cerr << "Error: invalid STL file" << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }

  try
  {
    // Try building a surface_mesh
    SMesh* SM = new SMesh();
    if (CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles))
    {
      auto pmap = boost::make_function_property_map<std::array<double, 3>, Point_3>(
                    [](const std::array<double, 3>& a) { return Point_3(a[0], a[1], a[2]); });
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, *SM,
                                                                  CGAL::parameters::point_map(pmap),
                                                                  CGAL::parameters::all_default());
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

  Scene_polygon_soup_item* item = new Scene_polygon_soup_item();
  item->setName(fileinfo.completeBaseName());
  item->load(points, triangles);
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>()<<item;
}

bool Polyhedron_demo_stl_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_surface_mesh_item*>(item);
}

bool Polyhedron_demo_stl_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
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
    = QInputDialog::getItem(NULL, tr("Save STL file"), tr("Format"), list, 0, false, &ok);

  if (!ok)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8(), std::ios::out | std::ios::binary);
  if ( choice == tr("Binary") )
    CGAL::set_mode(out, CGAL::IO::BINARY);
  else
  {
    CGAL::set_mode(out, CGAL::IO::ASCII);
    out.precision (std::numeric_limits<double>::digits10 + 2);
  }

  if (sm_item)
  {
    CGAL::write_STL(*sm_item->face_graph(), out);
    items.pop_front();
    return true;
  }
  return false;
}

#include "STL_io_plugin.moc"
