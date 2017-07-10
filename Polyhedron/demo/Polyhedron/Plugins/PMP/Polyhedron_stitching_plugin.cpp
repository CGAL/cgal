#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>

template <typename G>
struct Is_border {
  const G& g;
  Is_border(const G& g)
    : g(g)
  {}

 template <typename Descriptor>
  bool operator()(const Descriptor& d) const {
   return is_border(d,g);
  }

  bool operator()(typename boost::graph_traits<G>::vertex_descriptor d) const {
    return is_border(d,g) != boost::none;
  }

};


using namespace CGAL::Three;
class Polyhedron_demo_polyhedron_stitching_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionDetectBorders;
  QAction* actionStitchBorders;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionDetectBorders << actionStitchBorders; }
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* /* m */)
  {
    scene = scene_interface;
    actionDetectBorders= new QAction(tr("Detect Boundaries"), mainWindow);
    actionStitchBorders= new QAction(tr("Stitch Duplicated Boundaries"), mainWindow);
    actionDetectBorders->setObjectName("actionDetectBorders");
    actionStitchBorders->setObjectName("actionStitchBorders");
    actionStitchBorders->setProperty("subMenuName", "Polygon Mesh Processing");
    actionDetectBorders->setProperty("subMenuName", "Polygon Mesh Processing");
    autoConnectActions();
  }

  bool applicable(QAction*) const {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(index)) ||
           qobject_cast<Scene_surface_mesh_item*>(scene->item(index)) )
        return true;
    }
    return false;
  }

  template <typename Item>
  void on_actionDetectBorders_triggered(Scene_interface::Item_id index);

  template <typename Item>
  void on_actionStitchBorders_triggered(Scene_interface::Item_id index);

public Q_SLOTS:
  void on_actionDetectBorders_triggered();
  void on_actionStitchBorders_triggered();

}; // end Polyhedron_demo_polyhedron_stitching_plugin


template <typename Poly>
struct Polyline_visitor
{
  Scene_polylines_item* new_item;
  typename boost::property_map<Poly, CGAL::vertex_point_t>::const_type vpm;

  Polyline_visitor(const Poly& poly, Scene_polylines_item* new_item)
    : new_item(new_item), vpm(get(CGAL::vertex_point,poly))
  {}

  void start_new_polyline()
  {
    new_item->polylines.push_back( Scene_polylines_item::Polyline() );
  }

  void add_node(typename boost::graph_traits<Poly>::vertex_descriptor vd)
  {
    
    new_item->polylines.back().push_back(get(vpm,vd));
  }

  void end_polyline(){}
};


template <typename Item>
void Polyhedron_demo_polyhedron_stitching_plugin::on_actionDetectBorders_triggered(Scene_interface::Item_id index)
{
  typedef typename Item::Face_graph  FaceGraph;
  Item* item = qobject_cast<Item*>(scene->item(index));

  if(item)
    {
      Scene_polylines_item* new_item = new Scene_polylines_item();

      FaceGraph* pMesh = item->polyhedron();
      normalize_border(*pMesh);


      typedef boost::filtered_graph<FaceGraph,Is_border<FaceGraph>, Is_border<FaceGraph> > BorderGraph;
      
      Is_border<FaceGraph> ib(*pMesh);
      BorderGraph bg(*pMesh,ib,ib);
      Polyline_visitor<FaceGraph> polyline_visitor(*pMesh, new_item); 
      CGAL::split_graph_into_polylines( bg,
                                        polyline_visitor,
                                        CGAL::internal::IsTerminalDefault() );

      
      if (new_item->polylines.empty())
        {
          delete new_item;
        }
      else
        {
          new_item->setName(tr("Boundary of %1").arg(item->name()));
          new_item->setColor(Qt::red);
          scene->addItem(new_item);
          new_item->invalidateOpenGLBuffers();
        }
    }
}

void Polyhedron_demo_polyhedron_stitching_plugin::on_actionDetectBorders_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()){
    on_actionDetectBorders_triggered<Scene_polyhedron_item>(index);
    on_actionDetectBorders_triggered<Scene_surface_mesh_item>(index);
  }
}

template <typename Item>
void Polyhedron_demo_polyhedron_stitching_plugin::on_actionStitchBorders_triggered(Scene_interface::Item_id index)
{
  Item* item =
    qobject_cast<Item*>(scene->item(index));

  if(item){
    typename Item::Face_graph* pMesh = item->polyhedron();
    CGAL::Polygon_mesh_processing::stitch_borders(*pMesh);
    item->invalidateOpenGLBuffers();
    scene->itemChanged(item);
  }
}


void Polyhedron_demo_polyhedron_stitching_plugin::on_actionStitchBorders_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()){
    on_actionStitchBorders_triggered<Scene_polyhedron_item>(index);
    on_actionStitchBorders_triggered<Scene_surface_mesh_item>(index);
  }
}
#include "Polyhedron_stitching_plugin.moc"
