#include "Scene_surface_mesh_item.h"

#include <QMainWindow>
#include <QObject>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include <CGAL/Three/Three.h>

#include <CGAL/IO/read_surf_trianglemesh.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/array.h>

#include "Color_map.h"
#include <fstream>
#include <boost/container/flat_set.hpp>

using namespace CGAL::Three;
class Surf_io_plugin:
    public QObject,
    public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "surf_io_plugin.json")

public:

  QString name() const { return "surf_io_plugin"; }
  QString nameFilters() const { return "Amira files (*.surf);;Amira binary files (*.surf.am)"; }
  bool canLoad(QFileInfo) const{ return true; }
  template<class FaceGraphItem>
  CGAL::Three::Scene_item* actual_load(QFileInfo fileinfo);
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*) { return false; }
  bool save(QFileInfo ,QList<CGAL::Three::Scene_item*>& ) { return false; }
};


QList<Scene_item*>
Surf_io_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene)
{
  Scene_item* item =
      actual_load<Scene_surface_mesh_item>(fileinfo);
  if(item)
  {
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
template< class FaceGraphItem>
CGAL::Three::Scene_item* Surf_io_plugin::actual_load(QFileInfo fileinfo)
{
  typedef typename FaceGraphItem::Face_graph FaceGraph;
  typedef typename boost::property_traits<
      typename boost::property_map<FaceGraph, boost::vertex_point_t>::type
      >::value_type Point_3;
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }

  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_group_item* item =
        new Scene_group_item(fileinfo.completeBaseName());
    return item;
  }
  std::vector<FaceGraph> patches;
  std::vector<MaterialData> material_data;
  CGAL::Bbox_3 grid_box;
  std::array<unsigned int, 3> grid_size = {{1, 1, 1}};
  boost::container::flat_set<Point_3> duplicated_points;
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
    FaceGraphItem *patch = new FaceGraphItem(patches[i]);
    patch->setName(QString("Patch #%1").arg(i));
    patch->setColor(colors_[i]);
    CGAL::Three::Three::scene()->addItem(patch);
    CGAL::Three::Three::scene()->changeGroup(patch, group);
  }
  return group;
}

#include "Surf_io_plugin.moc"
