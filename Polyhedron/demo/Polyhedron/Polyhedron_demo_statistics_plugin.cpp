#include <QtCore/qglobal.h>

#include "Scene_item.h"
#include "Scene_interface.h"
#include <CGAL/gl.h>

#include <QAction>
#include <QMainWindow>
#include <QDialog>
#include <QString>
#include "Viewer_interface.h"
#include "Polyhedron_demo_plugin_helper.h"

#include <Polyhedron_type.h>
#include <Scene_polyhedron_item.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/property_map/property_map.hpp>
#include <map>

#include <boost/foreach.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include "ui_Polyhedron_demo_statistics_on_polyhedron_plugin.h"



namespace PMP = CGAL::Polygon_mesh_processing;
using namespace boost::accumulators;


class Polyhedron_demo_statistics_plugin :
        public QObject,
        public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionStatistics = new QAction(tr("Statistics"), mainWindow);
    connect(actionStatistics, SIGNAL(triggered()), this, SLOT(statistics()));
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionStatistics;
  }

  bool applicable(QAction*) const
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if (poly_item)
      return true;

    return false;
  }

  void statistics_on_polyhedron(Polyhedron* poly);
  //todo : statistics_on_polygon_soup
  //todo : statistics_on_c3t3
  //...
  void angles(Polyhedron* poly, double& mini, double& maxi, double& ave);

public Q_SLOTS:
  void statistics();

private:
  QAction* actionStatistics;

}; // end Polyhedron_demo_statistics_plugin


void Polyhedron_demo_statistics_plugin::statistics()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* poly_item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if (poly_item)
  {
    statistics_on_polyhedron(poly_item->polyhedron());
    return;
  }
}

void Polyhedron_demo_statistics_plugin::statistics_on_polyhedron(Polyhedron* poly)
{
  QDialog dlg(mw);
  Ui::Polyhedron_demo_statistics_on_polyhedron_dialog ui;
  ui.setupUi(&dlg);
  connect(ui.okButtonBox, SIGNAL(accepted()), &dlg, SLOT(accept()));

  ui.label_nbvertices->setText(QString::number(poly->size_of_vertices()));
  ui.label_nbfacets->setText(QString::number(poly->size_of_facets()));
  ui.label_nbborderedges->setText(QString::number(poly->size_of_border_edges()));
  
  typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
  int i = 0;
  BOOST_FOREACH(face_descriptor f, faces(*poly)){
    f->id() = i++;
  }

  boost::vector_property_map<int,
    boost::property_map<Polyhedron, boost::face_index_t>::type>
      fccmap(get(boost::face_index, *poly));
  ui.label_nbconnectedcomponents->setText(
    QString::number(PMP::connected_components(*poly, fccmap)));

  double mini, maxi, ave;
  angles(poly, mini, maxi, ave);
  ui.label_minangle->setText(QString::number(mini));
  ui.label_maxangle->setText(QString::number(maxi));
  ui.label_averageangle->setText(QString::number(ave));

  int x = dlg.exec();
//  if (i == QDialog::Rejected) { return; }
}

void Polyhedron_demo_statistics_plugin::angles(Polyhedron* poly,
  double& mini, double& maxi, double& ave)
{
  double rad_to_deg = 180. / CGAL_PI;

  accumulator_set< double,
                   features< tag::min, tag::max, tag::mean > > acc;

  boost::property_map<Polyhedron, CGAL::vertex_point_t>::type
    vpmap = get(CGAL::vertex_point, *poly);
  typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(*poly))
  {
    Kernel::Point_3 a = get(vpmap, source(h, *poly));
    Kernel::Point_3 b = get(vpmap, target(h, *poly));
    Kernel::Point_3 c = get(vpmap, target(next(h, *poly), *poly));

    Kernel::Vector_3 ba(b, a);
    Kernel::Vector_3 bc(b, c);
    double cos_angle = (ba * bc)
      / std::sqrt(ba.squared_length() * bc.squared_length());

    acc(std::acos(cos_angle) * rad_to_deg);
  }

  mini = extract_result< tag::min >(acc);
  maxi = extract_result< tag::max >(acc);
  ave  = extract_result< tag::mean >(acc);
}

#include "Polyhedron_demo_statistics_plugin.moc"
