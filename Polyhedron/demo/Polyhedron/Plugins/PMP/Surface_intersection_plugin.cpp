#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>

#ifdef USE_SURFACE_MESH
#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#else
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#endif
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polylines_item.h"

#include <boost/foreach.hpp>

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_face_graph_item;
#else
typedef Scene_polyhedron_item Scene_face_graph_item;
#endif


using namespace CGAL::Three;
namespace PMP = CGAL::Polygon_mesh_processing;

class Polyhedron_demo_intersection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  bool applicable(QAction*) const {
    return scene->selectionIndices().size() == 2 &&
      qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().first())) &&
      qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().last()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionPolyhedronIntersection_3;
  }

  void init(QMainWindow* mw, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    this->scene = scene_interface;
    actionPolyhedronIntersection_3 = new QAction("Surface Intersection", mw);
    actionPolyhedronIntersection_3->setProperty("subMenuName", "Polygon Mesh Processing");
    if(actionPolyhedronIntersection_3) {
      connect(actionPolyhedronIntersection_3, SIGNAL(triggered()),
              this, SLOT(intersection()));
    }
  }

private:

  QAction*  actionPolyhedronIntersection_3;
  Scene_interface *scene;

public Q_SLOTS:
  void intersection();

}; // end class Polyhedron_demo_intersection_plugin

void Polyhedron_demo_intersection_plugin::intersection()
{
  Scene_face_graph_item* itemA = NULL;
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_face_graph_item* itemB =
      qobject_cast<Scene_face_graph_item*>(scene->item(index));

    if(itemB)
    {
      if (itemA==NULL)
      {
        itemA = itemB;
        continue;
      }

      QApplication::setOverrideCursor(Qt::WaitCursor);

      Scene_polylines_item* new_item = new Scene_polylines_item();
     // perform Boolean operation
      QTime time;
      time.start();

      try{
        PMP::surface_intersection(*itemA->polyhedron(),
                                  *itemB->polyhedron(),
                                  std::back_inserter(new_item->polylines),
                                  true);
      }
      catch(CGAL::Corefinement::Self_intersection_exception)
      {
        QMessageBox::warning((QWidget*)NULL,
          tr("Self-intersections Found"),
          tr("Some self-intersections were found amongst intersecting facets"));
        delete new_item;
        QApplication::restoreOverrideCursor();
        return;
      }

      QString name = tr("%1 intersection %2");

      new_item->setName(name.arg(itemA->name(), itemB->name()));
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      if (new_item->polylines.empty())
        delete new_item;
      else{
        new_item->setColor(Qt::green);
        new_item->setRenderingMode(Wireframe);
        scene->addItem(new_item);
        new_item->invalidateOpenGLBuffers();
      }

      QApplication::restoreOverrideCursor();
    }
  }
}

#include "Surface_intersection_plugin.moc"
