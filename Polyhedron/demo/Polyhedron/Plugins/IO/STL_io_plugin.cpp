#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Kernel_type.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <fstream>

#include <CGAL/IO/Polyhedron_builder_from_STL.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>


#include <QColor>
#include <QMainWindow>
using namespace CGAL::Three;
class Polyhedron_demo_stl_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*) {
    //get the references
    this->scene = scene_interface;
    this->mw = mainWindow;
  }
  QList<QAction*> actions() const {
    return QList<QAction*>();
  }
  bool applicable(QAction*) const { return false;}
  QString nameFilters() const;
  QString name() const { return "stl_plugin"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
};

QString Polyhedron_demo_stl_plugin::nameFilters() const {
  return "STL files (*.stl)";
}

bool Polyhedron_demo_stl_plugin::canLoad() const {
  return true;
}


CGAL::Three::Scene_item*
Polyhedron_demo_stl_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios::in | std::ios::binary);
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }

  std::vector<CGAL::cpp11::array<double, 3> > points;
  std::vector<CGAL::cpp11::array<int, 3> > triangles;
  if (!CGAL::read_STL(in, points, triangles))
  {
    std::cerr << "Error: invalid STL file" << std::endl;
    return NULL;
  }

  try{
    if(this->mw->property("is_polyhedron_mode").toBool())
    {
      // Try building a polyhedron
      Polyhedron P;
      if (CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles))
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, P);

      if(! P.is_valid() || P.empty()){
        std::cerr << "Error: Invalid facegraph" << std::endl;
      }
      else{
        Scene_polyhedron_item* item = new Scene_polyhedron_item(P);
        item->setName(fileinfo.completeBaseName());
        return item;
      }
    }
    else
    {
      // Try building a surface_mesh
      SMesh* SM = new SMesh();
      if (CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles))
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, *SM);
      if(!SM->is_valid() || SM->is_empty()){
        std::cerr << "Error: Invalid facegraph" << std::endl;
      }
      else{
        Scene_surface_mesh_item* item = new Scene_surface_mesh_item(SM);
        item->setName(fileinfo.completeBaseName());
        return item;
      }
    }
  }
  catch(...){}

  Scene_polygon_soup_item* item = new Scene_polygon_soup_item();
  item->setName(fileinfo.completeBaseName());
  item->load(points, triangles);
  return item;
}

bool Polyhedron_demo_stl_plugin::canSave(const CGAL::Three::Scene_item*)
{
  return false;
}

bool Polyhedron_demo_stl_plugin::save(const CGAL::Three::Scene_item*, QFileInfo)
{
  return false;
}

#include "STL_io_plugin.moc"
