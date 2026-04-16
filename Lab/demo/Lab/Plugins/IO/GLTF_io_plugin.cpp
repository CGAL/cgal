#include "SMesh_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Scene.h"

#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include <CGAL/IO/GLTF.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <QString>
#include <QStringList>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Three;

class CGAL_Lab_gltf_plugin :
  public QObject,
  public CGAL_Lab_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90")

public:
  QString name() const { return "gltf_plugin"; }
  QString nameFilters() const { return "GLTF files (*.gltf *.glb)"; }
  bool canLoad(QFileInfo) const { return true; }
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene = true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo, QList<CGAL::Three::Scene_item*>&);
};

QList<Scene_item*>
CGAL_Lab_gltf_plugin::load(QFileInfo fileinfo, bool& ok, bool add_to_scene)
{
  std::vector<EPICK::Point_3> points;
  std::vector<std::vector<std::size_t>> polygons;

  if(!CGAL::IO::read_GLTF(fileinfo.filePath().toUtf8().toStdString(),
                           points, polygons))
  {
    std::cerr << "Error: cannot read GLTF file "
              << fileinfo.filePath().toUtf8().toStdString() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }

  // Try to build a surface mesh from the polygon soup
  try
  {
    SMesh* sm = new SMesh();
    if(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons))
    {
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
        points, polygons, *sm);
    }

    if(sm->is_valid() && !sm->is_empty())
    {
      Scene_surface_mesh_item* item = new Scene_surface_mesh_item(sm);
      item->setName(fileinfo.completeBaseName());
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>() << item;
    }
    delete sm;
  }
  catch(...) {}

  // Fall back to polygon soup
  Scene_polygon_soup_item* item = new Scene_polygon_soup_item();
  item->setName(fileinfo.completeBaseName());
  item->load(points, polygons);
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>() << item;
}

bool CGAL_Lab_gltf_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_surface_mesh_item*>(item) ||
         qobject_cast<const Scene_polygon_soup_item*>(item);
}

bool CGAL_Lab_gltf_plugin::save(QFileInfo fileinfo,
                                 QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();

  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(item);
  const Scene_polygon_soup_item* soup_item =
    qobject_cast<const Scene_polygon_soup_item*>(item);

  if(!sm_item && !soup_item)
    return false;

  const std::string fname = fileinfo.filePath().toUtf8().toStdString();

  if(sm_item)
  {
    // Convert surface mesh to polygon soup and write
    std::vector<EPICK::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;
    const SMesh& mesh = *sm_item->face_graph();
    for(auto v : mesh.vertices())
      points.push_back(mesh.point(v));
    for(auto f : mesh.faces())
    {
      std::vector<std::size_t> poly;
      for(auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh))
        poly.push_back(static_cast<std::size_t>(v));
      polygons.push_back(poly);
    }
    if(CGAL::IO::write_GLTF(fname, points, polygons))
    {
      items.pop_front();
      return true;
    }
    return false;
  }

  if(soup_item)
  {
    if(CGAL::IO::write_GLTF(fname, soup_item->points(), soup_item->polygons()))
    {
      items.pop_front();
      return true;
    }
    return false;
  }

  return false;
}

#include "GLTF_io_plugin.moc"