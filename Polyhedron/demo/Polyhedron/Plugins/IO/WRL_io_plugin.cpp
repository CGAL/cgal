#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Kernel_type.h"
#include "Scene.h"
#include "SMesh_type.h"
#include <CGAL/IO/WRL_reader.h>
#include <CGAL/Timer.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <fstream>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <QColor>
#include <QMainWindow>
using namespace CGAL::Three;

class WRL_io_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface,
  public Polyhedron_demo_plugin_helper

{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "gocad_io_plugin.json")
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0" )

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
  QString name() const { return "wrl_plugin"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
};

QString WRL_io_plugin::nameFilters() const {
  return "VRML files (*.wrl)";
}

bool WRL_io_plugin::canLoad() const {
  return true;
}


CGAL::Three::Scene_item*
WRL_io_plugin::load(QFileInfo fileinfo) {
  
  typedef CGAL::cpp11::array<double, 3>       Point;
  typedef std::vector<std::size_t>            Face;
  
  
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
  
  // Try to read WRL file in a surface_mesh
  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_surface_mesh_item* item = new Scene_surface_mesh_item();
    item->setName(fileinfo.baseName());
    return item;
  }
  
  std::string name, color;
  std::vector<Point> points;
  std::vector<Face> facets;
  CGAL::IO::internal::read_WRL_and_merge_shapes(in, points, facets);
  if(points.empty())
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_surface_mesh_item* item = new Scene_surface_mesh_item();
    item->setName(fileinfo.baseName());
    return item;
  }
  Scene_item* item = nullptr;
  if(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(facets))
  {
    SMesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points,facets,
                                                                mesh);
    item = new Scene_surface_mesh_item(mesh);
  }
  else
  {
    item = new Scene_polygon_soup_item();
    std::vector<Point_3> cgal_points;
    BOOST_FOREACH(const Point &p, points)
    {
      cgal_points.push_back(Point_3(p[0], p[1],p[2]));
    }
    qobject_cast<Scene_polygon_soup_item*>(item)->load(cgal_points, facets);
  }
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
  return item;
}

bool WRL_io_plugin::canSave(const CGAL::Three::Scene_item* )
{
  return false;
}

bool WRL_io_plugin::save(const CGAL::Three::Scene_item* , QFileInfo )
{
  return false;
}

#include "WRL_io_plugin.moc"
