#include <QTime>
#include <QApplication>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/convex_hull_3.h>

class Polyhedron_demo_convex_hull_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionConvexHull";
  }

  bool applicable() const {
    return 
      qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public slots:
  void on_actionConvexHull_triggered();

}; // end Polyhedron_demo_convex_hull_plugin


void Polyhedron_demo_convex_hull_plugin::on_actionConvexHull_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
  
  Scene_polylines_item* lines_item = 
    qobject_cast<Scene_polylines_item*>(scene->item(index));
  
  if(poly_item || pts_item || lines_item)
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    QTime time;
    time.start();
    std::cout << "Convex hull...";

    // add convex hull as new polyhedron
    Polyhedron *pConvex_hull = new Polyhedron;
    if ( poly_item ){
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

Q_EXPORT_PLUGIN2(Polyhedron_demo_convex_hull_plugin, Polyhedron_demo_convex_hull_plugin)

#include "Polyhedron_demo_convex_hull_plugin.moc"
