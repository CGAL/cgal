#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/Three.h>
#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include <CGAL/exceptions.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <QMessageBox>
#include <QApplication>

#include <iostream>
#include <fstream>
#include <limits>
#include <memory>

using namespace CGAL::Three;

class CGAL_Lab_obj_plugin :
  public QObject,
  public CGAL_Lab_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90" FILE "obj_io_plugin.json")

public:
  bool isDefaultLoader(const Scene_item *item) const override
  {
    if(qobject_cast<const Scene_surface_mesh_item*>(item)
       || qobject_cast<const Scene_polygon_soup_item*>(item))
      return true;
    return false;
  }
  bool isDefaultLoader(const QString& name) const override
  {
    if(name == QString("obj"))
      return true;
    return false;
  }
  QString name() const override{ return "obj_plugin"; }
  QString nameFilters() const override { return "Wavefront OBJ (*.obj)"; }
  bool canLoad(QFileInfo fileinfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override;
  CGAL::Three::Scene_item* load_obj(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*) override;
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& ) override;
};

bool CGAL_Lab_obj_plugin::canLoad(QFileInfo) const {
  return true;
}

QList<Scene_item*> CGAL_Lab_obj_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {

  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_surface_mesh_item* item =
        new Scene_surface_mesh_item(SMesh());
    item->setName(fileinfo.completeBaseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }
  if(fileinfo.suffix().toLower() == "obj"){
    Scene_item* item = load_obj(fileinfo);
    if(item){
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>()<<item;
    }
    else
    {
      ok = false;
      return QList<Scene_item*>();
    }
  }
  return QList<Scene_item*>();
}

CGAL::Three::Scene_item*
CGAL_Lab_obj_plugin::load_obj(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return nullptr;
  }

  auto set_name = [&](CGAL::Three::Scene_item* item){
    item->setName(fileinfo.completeBaseName());
    return item;
  };

  auto mesh_item = std::make_unique<Scene_surface_mesh_item>();
  if(CGAL::IO::read_OBJ(in, *(mesh_item->polyhedron()), CGAL::parameters::verbose(true)))
    return set_name(mesh_item.release());

  in.clear();
  in.seekg(0, std::ios::beg);

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > polylines;
  std::vector<std::vector<std::size_t> > polygons;
  if(!CGAL::IO::internal::read_OBJ(in, points, polygons, polylines,
                                   CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(),
                                   true /*verbose*/)) {
    return nullptr;
  }

  if(!polygons.empty()) {
    auto soup_item = std::make_unique<Scene_polygon_soup_item>();
    soup_item->load(points, polygons);
    return set_name(soup_item.release());
  }

  if(!polylines.empty()) {
    std::list<std::vector<Point_3> > item_polylines;
    for(const std::vector<std::size_t>& pl : polylines) {
      std::vector<Point_3> item_pl;
      item_pl.reserve(pl.size());
      for(const std::size_t& pi : pl) {
        item_pl.push_back(points[pi]);
      }
      item_polylines.push_back(std::move(item_pl));
    }
    auto polyline_item = std::make_unique<Scene_polylines_item>();
    polyline_item->polylines = std::move(item_polylines);
    return set_name(polyline_item.release());
  }

  return nullptr;
}

bool CGAL_Lab_obj_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_surface_mesh_item*>(item) ||
    qobject_cast<const Scene_polygon_soup_item*>(item);
}

bool
CGAL_Lab_obj_plugin::save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(item);
  const Scene_polygon_soup_item* soup_item =
    qobject_cast<const Scene_polygon_soup_item*>(item);

  if(!sm_item && !soup_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());
  out.precision (std::numeric_limits<double>::digits10 + 2);

  if(fileinfo.suffix().toLower() == "obj"){
    bool res = false;
    if(sm_item)
    {
      std::optional<SMesh::template Property_map<SMesh::Vertex_index, EPICK::Vector_3>> vnormals
        = sm_item->face_graph()->template property_map<SMesh::Vertex_index, EPICK::Vector_3>("v:normal");
      if(vnormals.has_value())
        res = CGAL::IO::write_OBJ(out, *(sm_item->face_graph()), CGAL::parameters::vertex_normal_map(vnormals.value()));
      else
        res = CGAL::IO::write_OBJ(out, *(sm_item->face_graph()));
    }
    else if(soup_item)
    {
      res = CGAL::IO::write_OBJ(out, soup_item->points(), soup_item->polygons());
    }

    if(res)
    {
      items.pop_front();
      return true;
    }
    return false;
  }
  return false;
}

#include "OBJ_io_plugin.moc"
