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
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>

using namespace CGAL::Three;
class Polyhedron_demo_interpolated_corrected_principal_curvatures_plugin :
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
}; // end Polyhedron_demo_interpolated_corrected_principal_curvatures_plugin


void compute(SMesh* sMesh,
             Scene_polylines_item* max_curv,
             Scene_polylines_item* min_curv,
             Scene_polylines_item* max_negative_curv,
             Scene_polylines_item* min_negative_curv)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel EpicKernel;
  typedef EpicKernel::Point_3 Point;
  typedef EpicKernel::Point_3 Point;
  typedef EpicKernel::Vector_3 Vector;
  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  typedef std::tuple<
      EpicKernel::FT,
      EpicKernel::FT,
      Vector,
      Vector
  > PrincipalCurvatureTuple;

  typename boost::property_map<SMesh, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, *sMesh);

  bool created = false;
  SMesh::Property_map<vertex_descriptor, PrincipalCurvatureTuple> principal_curvature_map;

  boost::tie(principal_curvature_map, created) = sMesh->add_property_map<vertex_descriptor, PrincipalCurvatureTuple>
      ("v:principal_curvature_map", { 0, 0,
              Vector(0,0,0),
              Vector(0,0,0)});
      assert(created);

  PMP::interpolated_corrected_principal_curvatures(
      *sMesh,
      principal_curvature_map
  );

  typename EpicKernel::FT max_curvature_magnitude_on_mesh = 0;
  for (vertex_descriptor v : vertices(*sMesh))
  {
      const PrincipalCurvatureTuple pc = principal_curvature_map[v];
      max_curvature_magnitude_on_mesh = std::max(max_curvature_magnitude_on_mesh, std::max(abs(get<0>(pc)), get<1>(pc)));
  }

  for(vertex_descriptor v : vertices(*sMesh))
  {
    std::vector<Point> points;

    // pick central point
    const Point& central_point = get(vpmap,v);
    points.push_back(central_point);

    // compute min edge len around central vertex
    // to scale the ribbons used to display the directions

    typedef EPICK::FT FT;

    const std::size_t n = CGAL::edges(*sMesh).size();

    EpicKernel::FT avg_edge_length = 0;
    if (n > 0) {
        for (auto e : CGAL::edges(*sMesh))
            avg_edge_length += PMP::edge_length(e, *sMesh);
        avg_edge_length /= n;
    }

    const PrincipalCurvatureTuple pc = principal_curvature_map[v];

    Vector umin = (std::get<0>(pc)/ max_curvature_magnitude_on_mesh) * std::get<2>(pc) * avg_edge_length;
    Vector umax = (std::get<1>(pc)/ max_curvature_magnitude_on_mesh) * std::get<3>(pc) * avg_edge_length;

    Scene_polylines_item::Polyline max_segment(2), min_segment(2);

    const double du = 0.4;

    min_segment[0] = central_point + du * umin;
    min_segment[1] = central_point - du * umin;
    max_segment[0] = central_point + du * umax;
    max_segment[1] = central_point - du * umax;

    (std::get<0>(pc) > 0 ? min_curv : min_negative_curv)->polylines.push_back(min_segment);
    (std::get<1>(pc) > 0 ? max_curv : max_negative_curv)->polylines.push_back(max_segment);
  }
}

void Polyhedron_demo_interpolated_corrected_principal_curvatures_plugin::on_actionEstimateCurvature_triggered()
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
