#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <QMainWindow>
#include <QObject>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/IO/read_surf_trianglemesh.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/array.h>

#include "Color_map.h"
#include <fstream>
#include <boost/container/flat_set.hpp>

using namespace CGAL::Three;
class Surf_io_plugin:
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
  QString name() const { return "surf_io_plugin"; }
  QString nameFilters() const { return "Amira files (*.surf)"; }
  bool canLoad() const{ return true; }
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*) { return false; }
  bool save(const CGAL::Three::Scene_item*, QFileInfo) { return false; }
};


CGAL::Three::Scene_item* Surf_io_plugin::load(QFileInfo fileinfo)
{
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }

  std::vector<Polyhedron> patches;
  std::vector<MaterialData> material_data;
  CGAL::Bbox_3 grid_box;
  CGAL::cpp11::array<unsigned int, 3> grid_size = {{1, 1, 1}};
  boost::container::flat_set<Polyhedron::Point_3> duplicated_points;
  read_surf(in, patches, material_data, grid_box, grid_size
          , std::inserter(duplicated_points, duplicated_points.end()));

  for(std::size_t i=0; i<material_data.size(); ++i)
  {
   std::cout<<"The patch #"<<i<<":\n  -inner region : material's id = "<<material_data[i].innerRegion.first<<" material's name = "
           <<material_data[i].innerRegion.second<<"\n  -outer region: material's id = "<<material_data[i].outerRegion.first<<" material's name = "
             <<material_data[i].outerRegion.second<<std::endl;
  }
  if (!duplicated_points.empty())
    std::cout << duplicated_points.size() << " points have been duplicated." << std::endl;
  
  std::vector<QColor> colors_;
  compute_color_map(QColor(100, 100, 255), static_cast<unsigned>(patches.size()),
                    std::back_inserter(colors_));
  Scene_group_item* group = new Scene_group_item(fileinfo.completeBaseName());
  for(std::size_t i=0; i<patches.size(); ++i)
  {
    Scene_polyhedron_item *patch = new Scene_polyhedron_item(patches[i]);
    patch->setName(QString("Patch #%1").arg(i));
    patch->setColor(colors_[i]);
    scene->addItem(patch);
    scene->changeGroup(patch, group);
  }
  return group;
}

#include "Surf_io_plugin.moc"
