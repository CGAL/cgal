#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include <limits>

#include "Scene.h"
#include <QApplication>

#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

class Polyhedron_demo_jet_fitting_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionEstimateCurvature";
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionEstimateCurvature_triggered();
}; // end Polyhedron_demo_jet_fitting_plugin

void Polyhedron_demo_jet_fitting_plugin::on_actionEstimateCurvature_triggered()
{
  // get active polyhedron
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if(!poly_item)
    return;

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Polyhedron* pMesh = poly_item->polyhedron();

  // types
  typedef CGAL::Monge_via_jet_fitting<Kernel> Fitting;
  typedef Fitting::Monge_form Monge_form;

  typedef Kernel::Point_3 Point;

  Scene_polylines_item* max_curv = new Scene_polylines_item;
  max_curv->setColor(Qt::red);
  max_curv->setName(tr("%1 (max curvatures)").arg(poly_item->name()));
  Scene_polylines_item* min_curv = new Scene_polylines_item;
  min_curv->setColor(Qt::green);
  min_curv->setName(tr("%1 (min curvatures)").arg(poly_item->name()));

  Polyhedron::Vertex_iterator v;
  for(v = pMesh->vertices_begin();
      v != pMesh->vertices_end();
      v++)
  {
    std::vector<Point> points;

    // pick central point
    const Point& central_point = v->point();
    points.push_back(central_point);

    // compute min edge len around central vertex
    // to scale the ribbons used to display the directions

    typedef Kernel::FT FT;

    FT min_edge_len = std::numeric_limits<FT>::infinity();
    Polyhedron::Halfedge_around_vertex_circulator he = v->vertex_begin();
    Polyhedron::Halfedge_around_vertex_circulator end = he;
    CGAL_For_all(he,end)
    {
      const Point& p = he->opposite()->vertex()->point();
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

  scene->addItem(max_curv);
  scene->addItem(min_curv);
  max_curv->invalidate_buffers();
  min_curv->invalidate_buffers();
  
  // default cursor
  QApplication::restoreOverrideCursor();
}

#include "Polyhedron_demo_jet_fitting_plugin.moc"
