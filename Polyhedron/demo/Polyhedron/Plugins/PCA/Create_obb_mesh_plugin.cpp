#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/Optimal_bounding_box/obb.h>

using namespace CGAL::Three;

typedef Scene_surface_mesh_item::Face_graph FaceGraph;
typedef Polyhedron::Point_3 Point_3;


class Create_obb_mesh_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const;

  bool applicable(QAction*) const {
    if(scene->mainSelectionIndex() != -1
       && scene->item(scene->mainSelectionIndex())->isFinite())
      return true;
  return false;}

protected:
  template <class FaceGraph_item>
  void gather_mesh_points(std::vector<Point_3>& points);
  void obb();

public Q_SLOTS:
  void createObb() {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    obb();
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionObb;

}; // end Create_obb_mesh_plugin class



void Create_obb_mesh_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionObb = new QAction(tr("Create &Optimal Bbox Mesh"), mainWindow);
  actionObb->setObjectName("createObbMeshAction");
  connect(actionObb, SIGNAL(triggered()), this, SLOT(createObb()));
}

QList<QAction*> Create_obb_mesh_plugin::actions() const {
  return QList<QAction*>() << actionObb;
}

template<class FaceGraph_item>
void Create_obb_mesh_plugin::gather_mesh_points(std::vector<Point_3>& points)
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  FaceGraph_item*  poly_item = qobject_cast<FaceGraph_item*>(scene->item(index));

  if(poly_item != NULL)
  {
    typedef typename FaceGraph_item::Face_graph FaceGraph;
    typedef typename boost::property_map<FaceGraph, CGAL::vertex_point_t>::type PointPMap;
    typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
    FaceGraph& pmesh = *poly_item->polyhedron();
    PointPMap pmap = get(boost::vertex_point, pmesh);
    BOOST_FOREACH(vertex_descriptor v, vertices(pmesh))
      points.push_back(get(pmap, v ));
  }
}

void Create_obb_mesh_plugin::obb()
{

  // gather point coordinates
  std::vector<Point_3> points;
  if(mw->property("is_polyhedron_mode").toBool())
  {
    gather_mesh_points<Scene_polyhedron_item>(points);
  }
  else
  {
    gather_mesh_points<Scene_surface_mesh_item>(points);
  }

  // find obb
  std::vector<Point_3> obb_points(8);
  CGAL::Optimal_bounding_box::find_obb(points, obb_points, true);

  Scene_item* item;
  if(mw->property("is_polyhedorn_mode").toBool()){
    Polyhedron* p = new Polyhedron;
    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                          obb_points[4], obb_points[5], obb_points[6], obb_points[7], *p);
    item = new Scene_polyhedron_item(p);
  } else {
    SMesh* p = new SMesh;
    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
                          obb_points[4], obb_points[5], obb_points[6], obb_points[7], *p);
    item = new Scene_surface_mesh_item(p);
  }

  item->setName("Optimal bbox mesh");
  item->setRenderingMode(Wireframe);
  scene->addItem(item);
}

#include "Create_obb_mesh_plugin.moc"
