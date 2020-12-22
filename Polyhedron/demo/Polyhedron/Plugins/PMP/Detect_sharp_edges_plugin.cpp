#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QStringList>
#include <QInputDialog>
#include <QtPlugin>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"

#include "Polyhedron_demo_detect_sharp_edges.h"

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef CGAL::Kernel_traits<Scene_surface_mesh_item::Face_graph::Point>::Kernel Kernel;

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
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "detect_sharp_edges_plugin.json")

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

namespace PMP = CGAL::Polygon_mesh_processing;
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
  std::size_t first_patch = 1;
  Q_FOREACH(Poly_tuple tuple, polyhedrons)
  {
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(tuple.first));
    FaceGraph* pMesh = tuple.second;
    if (!pMesh)
      continue;

    typedef boost::property_map<FaceGraph,CGAL::face_patch_id_t<int> >::type PatchID;
    typedef boost::property_map<FaceGraph, CGAL::vertex_incident_patches_t<int> >::type VIP;
    boost::property_map<FaceGraph, CGAL::edge_is_feature_t>::type eif
      = get(CGAL::edge_is_feature, *pMesh);
    PatchID pid = get(CGAL::face_patch_id_t<int>(), *pMesh);
    VIP vip = get(CGAL::vertex_incident_patches_t<int>(), *pMesh);

    first_patch+=PMP::sharp_edges_segmentation(*pMesh, angle, eif, pid,
                                               PMP::parameters::first_index(first_patch)
                                               .vertex_incident_patches_map(vip));
    //update item
    item->setItemIsMulticolor(true);
    item->computeItemColorVectorAutomatically(true);
    item->invalidateOpenGLBuffers();

    // update scene
    scene->itemChanged(tuple.first);
  }

  // default cursor
  QApplication::restoreOverrideCursor();
}

#include "Detect_sharp_edges_plugin.moc"
