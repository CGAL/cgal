#include <QApplication>
#include <QAction>
#include <QMessageBox>
#include <QMainWindow>
#include "opengl_tools.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
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
    return
        (qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
         qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex())));
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
  typedef Scene_surface_mesh_item::SMesh Surface_mesh;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor Face_descriptor;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
  typedef Surface_mesh::Vertex_index Vertex_index;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  bool found = false;
  std::vector<Scene_polyhedron_item*> selected_polys;
  std::vector<Scene_surface_mesh_item*> selected_sms;
  Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polyhedron_item* poly_item =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(poly_item)
    {
      selected_polys.push_back(poly_item);
    }
    else
    {
      Scene_surface_mesh_item* sm_item =
          qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
      if(sm_item)
      {
        selected_sms.push_back(sm_item);
      }
    }
  }
  Q_FOREACH(Scene_polyhedron_item* poly_item, selected_polys)
  {
    typedef boost::graph_traits<Polyhedron>::face_descriptor Face_descriptor;
    Polyhedron* pMesh = poly_item->polyhedron();
    std::vector<std::pair<Face_descriptor, Face_descriptor> > facets;
    // add intersecting triangles to a selection_item.
    if(selfIntersect(pMesh, facets))
    {
      Scene_polyhedron_selection_item* selection_item = new Scene_polyhedron_selection_item(poly_item, mw);
      for(std::vector<std::pair<Face_descriptor, Face_descriptor> >::iterator fb = facets.begin();
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
      selection_item->setName(tr("%1 (selection) (intersecting triangles)").arg(poly_item->name()));
      poly_item->setRenderingMode(Wireframe);
      scene->addItem(selection_item);
      scene->itemChanged(poly_item);
      scene->itemChanged(selection_item);
      found = true;
    }
  }

  Q_FOREACH(Scene_surface_mesh_item* sm_item, selected_sms)
  {
    Surface_mesh* sMesh = sm_item->polyhedron();
    std::vector<std::pair<Face_descriptor, Face_descriptor> > faces;
    // add intersecting triangles to a new Surface_mesh.
    if(selfIntersect(sMesh, faces))
    {
      //add the faces
      Surface_mesh* new_mesh = new Surface_mesh();
      for(std::vector<std::pair<Face_descriptor, Face_descriptor> >::iterator fb = faces.begin();
          fb != faces.end(); ++fb)
      {
        std::vector<Vertex_index> face_vertices;
        BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fb->first, *sMesh), *sMesh))
        {
          face_vertices.push_back(new_mesh->add_vertex(get(sMesh->points(), target(he_circ, *sMesh))));
        }
        new_mesh->add_face(face_vertices);
        face_vertices.clear();
        BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fb->second, *sMesh), *sMesh))
        {
          face_vertices.push_back(new_mesh->add_vertex(get(sMesh->points(), target(he_circ, *sMesh))));
        }
        new_mesh->add_face(face_vertices);
      }
      //add the edges
      BOOST_FOREACH(Face_descriptor fd, new_mesh->faces())
      {
        BOOST_FOREACH(halfedge_descriptor he_circ, halfedges_around_face( halfedge(fd, *new_mesh), *new_mesh))
        {
          new_mesh->add_edge(source(he_circ, *new_mesh), target(he_circ, *new_mesh));
        }
      }
      Scene_surface_mesh_item* ism_item = new Scene_surface_mesh_item(new_mesh);
      ism_item->invalidateOpenGLBuffers();
      ism_item->setName(tr("%1(intersecting triangles)").arg(sm_item->name()));
      ism_item->setColor(QColor(Qt::gray).darker(200));
      sm_item->setRenderingMode(Wireframe);
      scene->addItem(ism_item);
      scene->itemChanged(sm_item);
      scene->itemChanged(ism_item);
      found = true;
    }
  }
  if(!found)
    QMessageBox::information(mw, tr("No self intersection"),
                             tr("None of the selected surfaces self-intersect."));
  QApplication::restoreOverrideCursor();
}

#include "Self_intersection_plugin.moc"
