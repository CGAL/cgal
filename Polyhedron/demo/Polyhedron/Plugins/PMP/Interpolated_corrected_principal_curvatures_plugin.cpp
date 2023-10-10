#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"

#include <limits>

#include "Scene.h"
#include <QApplication>
#include <QMainWindow>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

using namespace CGAL::Three;
class Polyhedron_demo_interpolated_corrected_principal_curvatures_and_directions_plugin :
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
  void init(QMainWindow* mw,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
      scene = scene_interface;
      QAction *actionEstimateCurvature = new QAction(tr("Interpolated Corrected Principal Curvatures"), mw);
      connect(actionEstimateCurvature, SIGNAL(triggered()), this, SLOT(on_actionEstimateCurvature_triggered()));
      _actions <<actionEstimateCurvature;

  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionEstimateCurvature_triggered();
private :
  Scene_interface *scene;
  QList<QAction*> _actions;
}; // end Polyhedron_demo_interpolated_corrected_principal_curvatures_and_directions_plugin


void compute(SMesh* sMesh,
             Scene_polylines_item* max_curv,
             Scene_polylines_item* min_curv,
             Scene_polylines_item* max_negative_curv,
             Scene_polylines_item* min_negative_curv)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
  typedef Epic_kernel::Point_3 Point;
  typedef Epic_kernel::Point_3 Point;
  typedef Epic_kernel::Vector_3 Vector;
  typedef boost::graph_traits<SMesh>::vertex_descriptor Vertex_descriptor;

  typename boost::property_map<SMesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, *sMesh);

  bool created = false;
  SMesh::Property_map<Vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>> principal_curvatures_and_directions_map;

  boost::tie(principal_curvatures_and_directions_map, created) = sMesh->add_property_map<Vertex_descriptor, PMP::Principal_curvatures_and_directions<Epic_kernel>>
    ("v:principal_curvatures_and_directions_map", { 0, 0,
        Vector(0,0,0),
        Vector(0,0,0) });
  assert(created);


  PMP::interpolated_corrected_curvatures(
    *sMesh,
    CGAL::parameters::ball_radius(0)
                     .vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map)
  );

  double max_curvature_magnitude_on_mesh = 0;
  for (Vertex_descriptor v : vertices(*sMesh))
  {
    const PMP::Principal_curvatures_and_directions<Epic_kernel> pc = principal_curvatures_and_directions_map[v];
    max_curvature_magnitude_on_mesh = std::max(max_curvature_magnitude_on_mesh, std::max(abs(pc.min_curvature), abs(pc.max_curvature)));
  }

  for(Vertex_descriptor v : vertices(*sMesh))
  {
    std::vector<Point> points;

    // pick central point
    const Point& central_point = get(vpmap,v);
    points.push_back(central_point);

    // compute min edge len around central vertex
    // to scale the ribbons used to display the directions

    const std::size_t n = CGAL::edges(*sMesh).size();

    double avg_edge_length = 0;
    if (n > 0) {
      for (auto e : CGAL::edges(*sMesh))
        avg_edge_length += PMP::edge_length(e, *sMesh);
      avg_edge_length /= n;
    }

    const PMP::Principal_curvatures_and_directions<Epic_kernel> pc = principal_curvatures_and_directions_map[v];

    Vector umin = (pc.min_curvature / max_curvature_magnitude_on_mesh) * pc.min_direction * avg_edge_length;
    Vector umax = (pc.max_curvature / max_curvature_magnitude_on_mesh) * pc.max_direction * avg_edge_length;

    Scene_polylines_item::Polyline max_segment(2), min_segment(2);

    const double du = 0.4;

    min_segment[0] = central_point + du * umin;
    min_segment[1] = central_point - du * umin;
    max_segment[0] = central_point + du * umax;
    max_segment[1] = central_point - du * umax;

    (pc.min_curvature > 0 ? min_curv : min_negative_curv)->polylines.push_back(min_segment);
    (pc.max_curvature > 0 ? max_curv : max_negative_curv)->polylines.push_back(max_segment);
  }
}

void Polyhedron_demo_interpolated_corrected_principal_curvatures_and_directions_plugin::on_actionEstimateCurvature_triggered()
{
  // get active polyhedron
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  QString name = scene->item(index)->name();
  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if(! sm_item){
    return;
  }
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);

  // types
  Scene_polylines_item* max_curv = new Scene_polylines_item;
  max_curv->setColor(Qt::red);
  max_curv->setWidth(3);
  max_curv->setName(tr("%1 (max curvatures)").arg(name));

  Scene_polylines_item* min_curv = new Scene_polylines_item;
  min_curv->setColor(QColor(255,210,0));
  min_curv->setWidth(4);
  min_curv->setName(tr("%1 (min curvatures)").arg(name));

  Scene_polylines_item* max_negative_curv = new Scene_polylines_item;
  max_negative_curv->setColor(Qt::cyan);
  max_negative_curv->setWidth(4);
  max_negative_curv->setName(tr("%1 (max negative curvatures)").arg(name));

  Scene_polylines_item* min_negative_curv = new Scene_polylines_item;
  min_negative_curv->setColor(Qt::blue);
  min_negative_curv->setWidth(3);
  min_negative_curv->setName(tr("%1 (min negative curvatures)").arg(name));

  SMesh* pMesh = sm_item->polyhedron();
  compute(pMesh, max_curv, min_curv, max_negative_curv, min_negative_curv);

  scene->addItem(max_curv);
  scene->addItem(min_curv);
  max_curv->invalidateOpenGLBuffers();
  min_curv->invalidateOpenGLBuffers();
  scene->addItem(max_negative_curv);
  scene->addItem(min_negative_curv);
  max_negative_curv->invalidateOpenGLBuffers();
  min_negative_curv->invalidateOpenGLBuffers();

  // default cursor
  QApplication::restoreOverrideCursor();
}

#include "Interpolated_corrected_principal_curvatures_plugin.moc"
