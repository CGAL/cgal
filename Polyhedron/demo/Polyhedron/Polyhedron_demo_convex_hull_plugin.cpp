#include <QTime>
#include <QApplication>
#include <QAction>
#include <QStringList>

#include "opengl_tools.h"
#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/convex_hull_3.h>
#include <boost/iterator/transform_iterator.hpp>

class Polyhedron_demo_convex_hull_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionConvexHull";
  }

  bool applicable(QAction*) const {
    return 
      qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionConvexHull_triggered();

}; // end Polyhedron_demo_convex_hull_plugin

// for transform iterator
struct Get_point {
  typedef const Polyhedron::Point_3& result_type;
  result_type operator()(const Polyhedron::Vertex_handle v) const
  { return v->point(); }
};

void Polyhedron_demo_convex_hull_plugin::on_actionConvexHull_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
  
  Scene_polylines_item* lines_item = 
    qobject_cast<Scene_polylines_item*>(scene->item(index));
  
  Scene_polyhedron_selection_item* selection_item = 
    qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

  if(poly_item || pts_item || lines_item || selection_item)
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    QTime time;
    time.start();
    std::cout << "Convex hull...";

    // add convex hull as new polyhedron
    Polyhedron *pConvex_hull = new Polyhedron;
    if(selection_item) {
      CGAL::convex_hull_3(
        boost::make_transform_iterator(selection_item->selected_vertices.begin(), Get_point()),
        boost::make_transform_iterator(selection_item->selected_vertices.end(), Get_point()),
        *pConvex_hull);
    }
    else if ( poly_item ){
      Polyhedron* pMesh = poly_item->polyhedron();  
      CGAL::convex_hull_3(pMesh->points_begin(),pMesh->points_end(),*pConvex_hull);
    }
    else{
      if (pts_item)
        CGAL::convex_hull_3(pts_item->point_set()->begin(),pts_item->point_set()->end(),*pConvex_hull);
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

    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pConvex_hull);
    new_item->setName(tr("%1 (convex hull)").arg(scene->item(index)->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(FlatPlusEdges);
    scene->addItem(new_item);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Polyhedron_demo_convex_hull_plugin.moc"
