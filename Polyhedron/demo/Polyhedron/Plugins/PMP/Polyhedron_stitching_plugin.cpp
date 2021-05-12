#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>

using namespace CGAL::Three;
class Polyhedron_demo_polyhedron_stitching_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "polyhedron_stitching_plugin.json")

  QAction* actionDetectBorders;
  QAction* actionStitchBorders;
  QAction* actionStitchByCC;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionDetectBorders << actionStitchBorders << actionStitchByCC; }
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* /* m */)
  {
    scene = scene_interface;
    actionDetectBorders= new QAction(tr("Detect Boundaries"), mainWindow);
    actionDetectBorders->setObjectName("actionDetectBorders");
    actionDetectBorders->setProperty("subMenuName", "Polygon Mesh Processing");
    actionStitchBorders= new QAction(tr("Stitch Duplicated Boundaries"), mainWindow);
    actionStitchBorders->setObjectName("actionStitchBorders");
    actionStitchBorders->setProperty("subMenuName", "Polygon Mesh Processing/Repair");

    actionStitchByCC= new QAction(tr("Stitch Borders Per Connected Components"), mainWindow);
    actionStitchByCC->setObjectName("actionStitchByCC");
    actionStitchByCC->setProperty("subMenuName", "Polygon Mesh Processing/Repair");

    autoConnectActions();
  }

  bool applicable(QAction*) const {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if ( qobject_cast<Scene_surface_mesh_item*>(scene->item(index)) )
        return true;
    }
    return false;
  }

  template <typename Item>
  void on_actionDetectBorders_triggered(Scene_interface::Item_id index);

  template <typename Item>
  void on_actionStitchBorders_triggered(Scene_interface::Item_id index);

  template <typename Item>
  void on_actionStitchByCC_triggered(Scene_interface::Item_id index);

public Q_SLOTS:
  void on_actionDetectBorders_triggered();
  void on_actionStitchBorders_triggered();
  void on_actionStitchByCC_triggered();

}; // end Polyhedron_demo_polyhedron_stitching_plugin



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
      for(auto ed : edges(*pMesh))
      {
        if(pMesh->is_border(ed))
        {
          new_item->polylines.push_back(Scene_polylines_item::Polyline());
          new_item->polylines.back().push_back(pMesh->point(pMesh->source(pMesh->halfedge(ed))));
          new_item->polylines.back().push_back(pMesh->point(pMesh->target(pMesh->halfedge(ed))));
        }
      }

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
    on_actionStitchBorders_triggered<Scene_surface_mesh_item>(index);
  }
}

template <typename Item>
void Polyhedron_demo_polyhedron_stitching_plugin::on_actionStitchByCC_triggered(Scene_interface::Item_id index)
{
  Item* item =
      qobject_cast<Item*>(scene->item(index));

  if(!item)
    return;
  CGAL::Polygon_mesh_processing::stitch_borders(*item->polyhedron(),
                                                CGAL::Polygon_mesh_processing::parameters::apply_per_connected_component(true));
  item->invalidateOpenGLBuffers();
  scene->itemChanged(item);
}


void Polyhedron_demo_polyhedron_stitching_plugin::on_actionStitchByCC_triggered()
{
  Q_FOREACH(int index, scene->selectionIndices()){
    on_actionStitchByCC_triggered<Scene_surface_mesh_item>(index);
  }
}
#include "Polyhedron_stitching_plugin.moc"
