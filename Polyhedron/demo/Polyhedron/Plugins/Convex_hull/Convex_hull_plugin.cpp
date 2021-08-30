#include <QElapsedTimer>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStringList>
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "SMesh_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
using namespace CGAL::Three;
class Polyhedron_demo_convex_hull_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
    void init(QMainWindow*mw,
              Scene_interface* scene_interface,
              Messages_interface*)
    {
        scene = scene_interface;
        this->mw = mw;
        QAction *actionConvexHull = new QAction("Convex Hull", mw);
        actionConvexHull->setProperty("subMenuName","3D Convex Hulls");
        connect(actionConvexHull, SIGNAL(triggered()), this, SLOT(on_actionConvexHull_triggered()));
        _actions <<actionConvexHull;
    }

  QList<QAction*> actions()const {return _actions;}

  bool applicable(QAction*) const {
    return
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionConvexHull_triggered();
private:
  QList<QAction*> _actions;
  Scene_interface* scene;
  QMainWindow* mw;
}; // end Polyhedron_demo_convex_hull_plugin


void Polyhedron_demo_convex_hull_plugin::on_actionConvexHull_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  Scene_polylines_item* lines_item =
    qobject_cast<Scene_polylines_item*>(scene->item(index));

  Scene_polyhedron_selection_item* selection_item =
    qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if( pts_item || lines_item || selection_item || sm_item)
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    QElapsedTimer time;
    time.start();
    std::cout << "Convex hull...";

    // add convex hull as new polyhedron
    SMesh *pConvex_hull  = new SMesh;
    if(selection_item) {
      SMesh* pMesh = selection_item->polyhedron();
      boost::property_map<Face_graph,CGAL::vertex_point_t>::const_type
          vpm = get(CGAL::vertex_point, *pMesh);
      std::vector<Kernel::Point_3> all_points;
      all_points.reserve( pMesh->number_of_vertices() );

      for(const auto& v_it : selection_item->selected_vertices) {
        all_points.push_back(get(vpm,v_it));
      }

      for(const auto& e_it : selection_item->selected_edges) {
        all_points.push_back(get(vpm,target(halfedge(e_it,*pMesh),*pMesh)));
        all_points.push_back(get(vpm,target(opposite(halfedge(e_it,*pMesh),*pMesh),*pMesh)));
      }

      for(const auto& fd : selection_item->selected_facets) {
        for(fg_halfedge_descriptor he : halfedges_around_face(halfedge(fd,*pMesh),*pMesh)){
          all_points.push_back(get(vpm,target(he,*pMesh)));
        }
      }
      all_points.shrink_to_fit();
      CGAL::convex_hull_3(
            all_points.begin(),
            all_points.end(),
            *pConvex_hull);
    }
    else if ( sm_item ){
      SMesh* pMesh = sm_item->polyhedron();
      CGAL::convex_hull_3(*pMesh,
                          *pConvex_hull);

    }
    else{
      if (pts_item)
        CGAL::convex_hull_3(pts_item->point_set()->points().begin(),
                            pts_item->point_set()->points().end(),
                            *pConvex_hull);
      else{
        std::size_t nb_points=0;
        for(std::list<std::vector<Kernel::Point_3> >::const_iterator it = lines_item->polylines.begin();
            it != lines_item->polylines.end();
            ++it)  nb_points+=it->size();

        std::vector<Kernel::Point_3> all_points;
        all_points.reserve( nb_points );

        for(std::list<std::vector<Kernel::Point_3> >::const_iterator it = lines_item->polylines.begin();
            it != lines_item->polylines.end();
            ++it)  std::copy(it->begin(), it->end(),std::back_inserter( all_points ) );

        CGAL::convex_hull_3(all_points.begin(),all_points.end(),*pConvex_hull);
      }
    }
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(pConvex_hull);
    new_item->setName(tr("%1 (convex hull)").arg(scene->item(index)->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(FlatPlusEdges);
    scene->addItem(new_item);

    // default cursor
    QApplication::restoreOverrideCursor();
  }


}

#include "Convex_hull_plugin.moc"
