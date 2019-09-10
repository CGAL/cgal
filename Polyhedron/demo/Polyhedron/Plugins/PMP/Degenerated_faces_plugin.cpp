#include <QApplication>
#include <QAction>
#include <QMessageBox>
#include <QMainWindow>
#include "opengl_tools.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#ifdef USE_SURFACE_MESH
#include "Scene_surface_mesh_item.h"
#else
#include "Scene_polyhedron_item.h"
#endif
#include "Scene_polyhedron_selection_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Make_triangle_soup.h>
#include <CGAL/Kernel_traits.h>
#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_facegraph_item;
#else
typedef Scene_polyhedron_item Scene_facegraph_item;
#endif
typedef Scene_facegraph_item::Face_graph Face_graph;
using namespace CGAL::Three;
class Degenerated_faces_plugin :
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
    QAction *actionDegenFaces = new QAction(tr("Select Degenerate Faces and Their Boundary Edges"), mw);
    actionDegenFaces->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionDegenFaces, SIGNAL(triggered()), this, SLOT(on_actionDegenFaces_triggered()));

    QAction *actionDegenEdges = new QAction(tr("Select Degenerate Edges and Their Boundary Vertices"), mw);
    actionDegenEdges->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionDegenEdges, SIGNAL(triggered()), this, SLOT(on_actionDegenEdges_triggered()));

    _actions <<actionDegenFaces
            <<actionDegenEdges;

  }

  bool applicable(QAction*) const {
    return
        qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionDegenFaces_triggered();
  void on_actionDegenEdges_triggered();
private:
  QList<QAction*> _actions;
  Scene_interface *scene;
  QMainWindow *mw;

}; // end Degenerated_faces_plugin

//pretty useless for now but could allow a huge factorization when a selection_item is
// available for SM_items
template<class Mesh>
bool isDegen(Mesh* mesh, std::vector<typename boost::graph_traits<Mesh>::face_descriptor> &out_faces)
{
  typedef typename boost::graph_traits<Mesh>::face_descriptor FaceDescriptor;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type Vpm;
  typedef typename boost::property_traits<Vpm>::value_type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  //filter non-triangle_faces
  BOOST_FOREACH(FaceDescriptor f, faces(*mesh))
  {
    if(is_triangle(halfedge(f, *mesh), *mesh)
       && is_degenerate_triangle_face(f, *mesh, get(boost::vertex_point, *mesh), Kernel()) )
      out_faces.push_back(f);
  }
  return !out_faces.empty();
}

template<class Mesh>
bool isDegen(Mesh* mesh, std::vector<typename boost::graph_traits<Mesh>::edge_descriptor> &out_edges)
{
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor HalfedgeDescriptor;
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type Vpm;
  typedef typename boost::property_traits<Vpm>::value_type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  Vpm vpm = get(boost::vertex_point, *mesh);

  BOOST_FOREACH(HalfedgeDescriptor h, halfedges(*mesh))
  {
    Point s(get(vpm, source(h, *mesh))),
        t(get(vpm, target(h, *mesh)));
    CGAL::Vector_3<Kernel> v(s, t);
    if(v.squared_length() == 0)
      out_edges.push_back(edge(h, *mesh));
  }
  return !out_edges.empty();
}
void Degenerated_faces_plugin::on_actionDegenFaces_triggered()
{

  typedef boost::graph_traits<Face_graph>::face_descriptor Face_descriptor;
  typedef boost::graph_traits<Face_graph>::halfedge_descriptor halfedge_descriptor;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  bool found = false;
  std::vector<Scene_facegraph_item*> selected_polys;
  Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_facegraph_item* poly_item =
        qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(poly_item)
    {
      selected_polys.push_back(poly_item);
    }
  }
  Q_FOREACH(Scene_facegraph_item* poly_item, selected_polys)
  {
    Face_graph* pMesh = poly_item->polyhedron();
    std::vector<Face_descriptor> facets;
    // add intersecting triangles to a selection_item.
    if(isDegen(pMesh, facets))
    {
      Scene_polyhedron_selection_item* selection_item = new Scene_polyhedron_selection_item(poly_item, mw);
      for(std::vector<Face_descriptor>::iterator fb = facets.begin();
          fb != facets.end(); ++fb) {
        selection_item->selected_facets.insert(*fb);
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(*fb, *pMesh), *pMesh))
        {
          selection_item->selected_edges.insert(edge(h, *selection_item->polyhedron()));
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
    QMessageBox::information(mw, tr("No degenerate triangle"),
                             tr("None of the selected surfaces has degenerate triangle faces."));
}

void Degenerated_faces_plugin::on_actionDegenEdges_triggered()
{

  typedef boost::graph_traits<Face_graph>::edge_descriptor edge_descriptor;
  QApplication::setOverrideCursor(Qt::WaitCursor);
  bool found = false;
  std::vector<Scene_facegraph_item*> selected_polys;
  Q_FOREACH(Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_facegraph_item* poly_item =
        qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(poly_item)
    {
      selected_polys.push_back(poly_item);
    }
  }
  Q_FOREACH(Scene_facegraph_item* poly_item, selected_polys)
  {
    Face_graph* pMesh = poly_item->polyhedron();
    std::vector<edge_descriptor> edges;
    // add intersecting triangles to a selection_item.
    if(isDegen(pMesh, edges))
    {
      Scene_polyhedron_selection_item* selection_item = new Scene_polyhedron_selection_item(poly_item, mw);
      for(std::vector<edge_descriptor>::iterator ed = edges.begin();
          ed != edges.end(); ++ed) {
        selection_item->selected_edges.insert(*ed);

        selection_item->selected_vertices.insert(source( halfedge(*ed, *selection_item->polyhedron()),
                                                         *selection_item->polyhedron()));
        selection_item->selected_vertices.insert(target( halfedge(*ed, *selection_item->polyhedron()),
                                                         *selection_item->polyhedron()));
      }
      selection_item->invalidateOpenGLBuffers();
      selection_item->setName(tr("%1 (selection) (null length edges)").arg(poly_item->name()));
      poly_item->setRenderingMode(Wireframe);
      scene->addItem(selection_item);
      scene->itemChanged(poly_item);
      scene->itemChanged(selection_item);
      found = true;
    }
  }
  QApplication::restoreOverrideCursor();
  if(!found)
    QMessageBox::information(mw, tr("No degenerate edge"),
                             tr("None of the selected surfaces has degenerate edges."));
}

#include "Degenerated_faces_plugin.moc"
