#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"

#include <limits>

#include "Scene.h"
#include <QApplication>
#include <QMainWindow>

#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
using namespace CGAL::Three;
class Polyhedron_demo_jet_fitting_plugin :
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
      QAction *actionEstimateCurvature = new QAction(tr("Curvature Estimation"), mw);
      actionEstimateCurvature->setProperty("subMenuName",
                                           "Estimation of Local Differential Properties");
      connect(actionEstimateCurvature, SIGNAL(triggered()), this, SLOT(on_actionEstimateCurvature_triggered()));
      _actions <<actionEstimateCurvature;

  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionEstimateCurvature_triggered();
private :
  Scene_interface *scene;
  QList<QAction*> _actions;
}; // end Polyhedron_demo_jet_fitting_plugin


template <typename Poly>
void compute(Poly* pMesh,
             Scene_polylines_item* min_curv,
             Scene_polylines_item* max_curv)
{

  typedef CGAL::Monge_via_jet_fitting<Kernel> Fitting;
  typedef Fitting::Monge_form Monge_form;

  typedef Kernel::Point_3 Point;

  typename boost::property_map<Poly, CGAL::vertex_point_t>::type vpmap = get(CGAL::vertex_point, *pMesh);

  BOOST_FOREACH(typename boost::graph_traits<Poly>::vertex_descriptor v, vertices(*pMesh))
  {
    std::vector<Point> points;

    // pick central point
    const Point& central_point = get(vpmap,v);
    points.push_back(central_point);

    // compute min edge len around central vertex
    // to scale the ribbons used to display the directions

    typedef Kernel::FT FT;

    FT min_edge_len = std::numeric_limits<FT>::infinity();
    BOOST_FOREACH(typename boost::graph_traits<Poly>::halfedge_descriptor he, halfedges_around_target(v, *pMesh))
    {
      const Point& p = get( vpmap, target(opposite(he, *pMesh ), *pMesh));
      points.push_back(p);
      FT edge_len = std::sqrt(CGAL::squared_distance(central_point,p));
      min_edge_len = edge_len < min_edge_len ? edge_len : min_edge_len; // avoids #undef min
    }

    if(points.size() > 5)
    {
      // estimate curvature by fitting
      Fitting monge_fit;
      const int dim_monge = 2;
      const int dim_fitting = 2;
      Monge_form monge_form = monge_fit(points.begin(),points.end(),dim_fitting,dim_monge);

      // make monge form comply with vertex normal (to get correct
      // orientation)
      typedef Kernel::Vector_3 Vector;
      Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v, *pMesh);
      monge_form.comply_wrt_given_normal(n);

      Vector umin = min_edge_len * monge_form.minimal_principal_direction();
      Vector umax = min_edge_len * monge_form.maximal_principal_direction();

      Scene_polylines_item::Polyline max_segment(2), min_segment(2);

      const double du = 0.2;

      max_segment[0] = central_point + du * umax;
      max_segment[1] = central_point - du * umax;
      min_segment[0] = central_point + du * umin;
      min_segment[1] = central_point - du * umin;

      max_curv->polylines.push_back(max_segment);
      min_curv->polylines.push_back(min_segment);
    }
  }

}

void Polyhedron_demo_jet_fitting_plugin::on_actionEstimateCurvature_triggered()
{
  // get active polyhedron
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  QString name = scene->item(index)->name();
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Scene_surface_mesh_item* sm_item = 
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
  if((!poly_item) && (! sm_item)){
    return;
  }
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);

  // types
  Scene_polylines_item* max_curv = new Scene_polylines_item;
  max_curv->setColor(Qt::red);
  max_curv->setName(tr("%1 (max curvatures)").arg(name));
  Scene_polylines_item* min_curv = new Scene_polylines_item;
  min_curv->setColor(Qt::green);
  min_curv->setName(tr("%1 (min curvatures)").arg(name));

  if(poly_item){
    Polyhedron* pMesh = poly_item->polyhedron();
    compute(pMesh, min_curv, max_curv);
  } else {
    
    SMesh* pMesh = sm_item->polyhedron();
    compute(pMesh, min_curv, max_curv);
  }
  scene->addItem(max_curv);
  scene->addItem(min_curv);
  max_curv->invalidateOpenGLBuffers();
  min_curv->invalidateOpenGLBuffers();
  
  // default cursor
  QApplication::restoreOverrideCursor();
}

#include "Jet_fitting_plugin.moc"
