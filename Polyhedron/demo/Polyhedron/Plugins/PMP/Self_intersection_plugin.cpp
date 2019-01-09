#include <QApplication>
#include <QAction>
#include <QMessageBox>
#include <QMainWindow>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Make_triangle_soup.h>


#include "Kernel_type.h"
#include "Scene_polyhedron_selection_item.h"

#include "Scene_surface_mesh_item.h"
typedef Scene_surface_mesh_item Scene_face_graph_item;

typedef Scene_face_graph_item::Face_graph Face_graph;
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
      QAction *actionSelfIntersection = new QAction(tr("Self-&Intersection Test"), mw);
      actionSelfIntersection->setProperty("subMenuName", "Polygon Mesh Processing");
      connect(actionSelfIntersection, SIGNAL(triggered()), this, SLOT(on_actionSelfIntersection_triggered()));
      _actions <<actionSelfIntersection;

  }

  bool applicable(QAction*) const { 
    return
        qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionSelfIntersection_triggered();
private:
  QList<QAction*> _actions;
  Scene_interface *scene;
  QMainWindow *mw;

}; // end Polyhedron_demo_self_intersection_plugin

//pretty useless for now but could allow a huge factorization when a selection_item is
// available for SM_items
template<class Mesh>
bool selfIntersect(Mesh* mesh, std::vector<std::pair<typename boost::graph_traits<Mesh>::face_descriptor,typename boost::graph_traits<Mesh>::face_descriptor> > &faces)
{

  // compute self-intersections
  CGAL::Polygon_mesh_processing::self_intersections
    (*mesh, std::back_inserter(faces),
    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, *mesh)));

  std::cout << "ok (" << faces.size() << " triangle pair(s))" << std::endl;
  return !faces.empty();
}

void Polyhedron_demo_self_intersection_plugin::on_actionSelfIntersection_triggered()
{
  typedef boost::graph_traits<Face_graph>::face_descriptor Face_descriptor;
  typedef boost::graph_traits<Face_graph>::halfedge_descriptor halfedge_descriptor;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  bool found = false;
  std::vector<Scene_face_graph_item*> selected_polys;
  Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_face_graph_item* poly_item =
        qobject_cast<Scene_face_graph_item*>(scene->item(index));
    if(poly_item)
    {
      selected_polys.push_back(poly_item);
    }
  }
  Q_FOREACH(Scene_face_graph_item* poly_item, selected_polys)
  {
    Face_graph* mesh = poly_item->face_graph();
    std::vector<std::pair<Face_descriptor, Face_descriptor> > faces;
    // add intersecting triangles to a new Surface_mesh.
    if(selfIntersect(mesh, faces))
    {
      //add the faces
      Scene_polyhedron_selection_item* selection_item = new Scene_polyhedron_selection_item(poly_item, mw);
      for(std::vector<std::pair<Face_descriptor, Face_descriptor> >::iterator fb = faces.begin();
          fb != faces.end(); ++fb) {
        selection_item->selected_facets.insert(fb->first);
        selection_item->selected_facets.insert(fb->second);

        //add the edges
        BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fb->first, *mesh), *mesh))
        {
          selection_item->selected_edges.insert(edge(he_circ, *mesh));
        }
        BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fb->second, *mesh), *mesh))
        {
          selection_item->selected_edges.insert(edge(he_circ, *mesh));
        }
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->setName(tr("%1 (selection) (intersecting triangles)").arg(poly_item->name()));
      poly_item->setRenderingMode(Wireframe);
      scene->addItem(selection_item);
      scene->itemChanged(poly_item);
      scene->itemChanged(selection_item);
      found = true;
    }
  }
  QApplication::restoreOverrideCursor();
  if(!found)
    QMessageBox::information(mw, tr("No self intersection"),
                             tr("None of the selected surfaces self-intersect."));
}

#include "Self_intersection_plugin.moc"
