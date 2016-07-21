#include <QApplication>
#include <QAction>
#include <QMessageBox>
#include <QMainWindow>
#include "opengl_tools.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Make_triangle_soup.h>

typedef Kernel::Triangle_3 Triangle;
using namespace CGAL::Three;
class Polyhedron_demo_self_intersection_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  QList<QAction*> actions() const {
    return _actions;
  }

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
      mw = mainWindow;
      scene = scene_interface;
      QAction *actionSelfIntersection = new QAction(tr("Self-&intersection"), mw);
      actionSelfIntersection->setProperty("subMenuName", "Polygon Mesh Processing");
      connect(actionSelfIntersection, SIGNAL(triggered()), this, SLOT(on_actionSelfIntersection_triggered()));
      _actions <<actionSelfIntersection;

  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionSelfIntersection_triggered();
private:
  QList<QAction*> _actions;
  Scene_interface *scene;
  QMainWindow *mw;

}; // end Polyhedron_demo_self_intersection_plugin

void Polyhedron_demo_self_intersection_plugin::on_actionSelfIntersection_triggered()
{
  CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Polyhedron* pMesh = item->polyhedron();

    // compute self-intersections

    typedef Polyhedron::Facet_handle Facet_handle;
    std::vector<std::pair<Facet_handle, Facet_handle> > facets;
    CGAL::Polygon_mesh_processing::self_intersections
      (*pMesh, std::back_inserter(facets),
      CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *pMesh)));

    std::cout << "ok (" << facets.size() << " triangle pair(s))" << std::endl;

    // add intersecting triangles as a new polyhedron, i.e., a triangle soup.
    if(!facets.empty())
    {
      Scene_polyhedron_selection_item* selection_item = new Scene_polyhedron_selection_item(item, mw);
      for(std::vector<std::pair<Facet_handle, Facet_handle> >::iterator fb = facets.begin();
        fb != facets.end(); ++fb) {
        selection_item->selected_facets.insert(fb->first);
        selection_item->selected_facets.insert(fb->second);

        Polyhedron::Halfedge_around_facet_circulator hc = (fb->first)->facet_begin(), cend = hc;
        CGAL_For_all(hc,cend) {
          selection_item->selected_edges.insert(edge(hc, *selection_item->polyhedron()));
        }
        hc = (fb->second)->facet_begin();
        cend = hc;
        CGAL_For_all(hc,cend) {
          selection_item->selected_edges.insert(edge(hc, *selection_item->polyhedron()));
        }
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->setName(tr("%1 (selection) (intersecting triangles)").arg(item->name()));
      item->setRenderingMode(Wireframe);
      scene->addItem(selection_item);
      scene->itemChanged(item);
      scene->itemChanged(selection_item);
    }
    else
      QMessageBox::information(mw, tr("No self intersection"),
                               tr("The polyhedron \"%1\" does not self-intersect.").
                               arg(item->name()));
  }
  QApplication::restoreOverrideCursor();
}

#include "Self_intersection_plugin.moc"
