#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QStringList>
#include <QInputDialog>
#include <QtPlugin>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#ifdef USE_SURFACE_MESH
#include "Scene_surface_mesh_item.h"
#include <CGAL/Mesh_3/properties_Surface_mesh.h>

#else
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Mesh_3/properties_Polyhedron_3.h>
#endif

#include "Polyhedron_demo_detect_sharp_edges.h"

#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef CGAL::Kernel_traits<Scene_surface_mesh_item::Face_graph::Point>::Kernel Kernel;
#else
typedef Scene_polyhedron_item Scene_facegraph_item;
#endif

typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;


using namespace CGAL::Three;
class Polyhedron_demo_detect_sharp_edges_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionSharEdges = new QAction("Detect Sharp Features", mw);
    actionSharEdges->setObjectName("detectSharpFeaturesAction");
    if(actionSharEdges) {
      connect(actionSharEdges, SIGNAL(triggered()),
              this, SLOT(detectSharpEdgesWithInputDialog()));
    }
  }

  bool applicable(QAction*) const {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      Scene_facegraph_item* item =
        qobject_cast<Scene_facegraph_item*>(scene->item(index));
      if (item) return true;
    }
    return false;
  }
  
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSharEdges;
  }

public Q_SLOTS:
void detectSharpEdges(bool input_dialog = false, double angle = 60);
  void detectSharpEdgesWithInputDialog();

protected:
  Kernel::Vector_3 facet_normal(face_descriptor f);
  bool is_sharp(halfedge_descriptor he);

private:
  QAction* actionSharEdges;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
}; // end Polyhedron_demo_detect_sharp_edges_plugin

void Polyhedron_demo_detect_sharp_edges_plugin::detectSharpEdgesWithInputDialog()
{
  detectSharpEdges(true);
}

void Polyhedron_demo_detect_sharp_edges_plugin::detectSharpEdges(bool input_dialog,
                                                                 double angle)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  typedef std::pair<int,FaceGraph*> Poly_tuple;
  
  // Get selected items
  QList<Poly_tuple> polyhedrons;
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(index));
    if(!item)
      return;
    
    FaceGraph* pMesh = item->polyhedron();
    if(!pMesh)
      return;
    item->show_feature_edges(true);
    polyhedrons << std::make_pair(index, pMesh);
  }

  QApplication::restoreOverrideCursor();
  if(input_dialog) {
    bool ok = true;
    angle = QInputDialog::getDouble(NULL, 
                                    tr("Sharp edges max angle"),
                                    tr("Angle in degrees between 0 and 180:"),
                                    angle, // value
                                    0.,          // min 
                                    180., // max
                                    2,          // decimals
                                    &ok);
    if(!ok) return;
  }
  // Detect edges
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QApplication::processEvents();
  Q_FOREACH(Poly_tuple tuple, polyhedrons)
  {
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(tuple.first));
    FaceGraph* pMesh = tuple.second;
    if (!pMesh) continue;

    CGAL::Polygon_mesh_processing::Detect_features_in_polyhedra<FaceGraph,
        int> detect_features;

    // Get sharp features
    detect_features.detect_sharp_edges(*pMesh, angle);
    detect_features.detect_surface_patches(*pMesh);
    detect_features.detect_vertices_incident_patches(*pMesh);

    //update item
    item->setItemIsMulticolor(true);
    item->invalidateOpenGLBuffers();

    // update scene
    scene->itemChanged(tuple.first);
  }

  // default cursor
  QApplication::restoreOverrideCursor();
}

#include "Detect_sharp_edges_plugin.moc"
