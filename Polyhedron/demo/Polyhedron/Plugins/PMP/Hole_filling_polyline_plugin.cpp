#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"
#include <CGAL/Three/Scene_interface.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Polyhedron_type.h"

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Timer.h>
#include <CGAL/array.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QInputDialog>
#include <QMessageBox>

#include <vector>
#include <algorithm>

#include <boost/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>


struct Face : public CGAL::cpp11::array<int,3>
{
  Face(int i, int j, int k)
  {
    (*this)[0] = i;
    (*this)[1] = j;
    (*this)[2] = k;
  } 
};
namespace PMP = CGAL::Polygon_mesh_processing;

using namespace CGAL::Three;
class Polyhedron_demo_hole_filling_polyline_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction *) const { return qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionHoleFillingPolyline; }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m){
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionHoleFillingPolyline = new QAction(tr("Polyline Hole Filling"), mw);
    connect(actionHoleFillingPolyline, SIGNAL(triggered()),
      this, SLOT(hole_filling_polyline_action()));
  }
private:
  struct Nop_functor {
    template<class T>
    void operator()(const T & /*t*/) const {}
  };
  typedef boost::function_output_iterator<Nop_functor> Nop_out;

  struct Get_handle {
    typedef Polyhedron::Facet_handle result_type;
    result_type operator()(Polyhedron::Facet& f) const
    { return f.halfedge()->facet(); }
  };
  
public Q_SLOTS:
  void hole_filling_polyline_action() {
    Scene_polylines_item* polylines_item = qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex()));
    if(!polylines_item) {
      print_message("Error: there is no selected polyline item!");
      return;
    } 

    bool also_refine;
    const double density_control_factor = 
      QInputDialog::getDouble(mw, tr("Density Control Factor"),
      tr("Density Control Factor (Cancel for not Refine): "), 1.41, 0.0, 100.0, 2, &also_refine);

    bool use_DT = 
      QMessageBox::Yes == QMessageBox::question(
      NULL, "Use Delaunay Triangulation", "Use Delaunay Triangulation ?", QMessageBox::Yes|QMessageBox::No);

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    std::size_t counter = 0;
    for(Scene_polylines_item::Polylines_container::iterator it = polylines_item->polylines.begin();
      it != polylines_item->polylines.end(); ++it, ++counter) 
    {
      if(it->front() != it->back()) { //not closed, skip it
        print_message("Warning: skipping not closed polyline!");
        continue; 
      } 
      if(it->size() < 4) { // no triangle, skip it (needs at least 3 + 1 repeat)
        print_message("Warning: skipping polyline which has less than 4 points!");
        continue; 
      }

      CGAL::Timer timer; timer.start();
      std::vector<Face> patch;
      CGAL::Polygon_mesh_processing::triangulate_hole_polyline(*it,
        std::back_inserter(patch),
        PMP::parameters::use_delaunay_triangulation(use_DT));
      print_message(QString("Triangulated in %1 sec.").arg(timer.time()));

      if(patch.empty()) {
        print_message("Warning: generating patch is not successful, please try it without 'Delaunay Triangulation'!");
        continue;
      }

      if(mw->property("is_polyhedron_mode").toBool()){
        Polyhedron* poly = new Polyhedron;
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(*it,
                                                                    patch,
                                                                    *poly);

        if(also_refine) {
          timer.reset();
          CGAL::Polygon_mesh_processing::refine(*poly, faces(*poly),
                                                Nop_out(), Nop_out(),
                                                CGAL::Polygon_mesh_processing::parameters::density_control_factor(density_control_factor));
          print_message(QString("Refined in %1 sec.").arg(timer.time()));
        }

        Scene_polyhedron_item* poly_item = new Scene_polyhedron_item(poly);
        poly_item->setName(tr("%1-filled-%2").arg(polylines_item->name()).arg(counter));
        poly_item->setRenderingMode(FlatPlusEdges);
        scene->setSelectedItem(scene->addItem(poly_item));
      } else {
        SMesh* poly = new SMesh;
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(*it,
                                                                    patch,
                                                                    *poly);

        if(also_refine) {
          timer.reset();
          CGAL::Polygon_mesh_processing::refine(*poly, faces(*poly),
                                                Nop_out(), Nop_out(),
                                                CGAL::Polygon_mesh_processing::parameters::density_control_factor(density_control_factor));
          print_message(QString("Refined in %1 sec.").arg(timer.time()));
        }

        Scene_surface_mesh_item* poly_item = new Scene_surface_mesh_item(poly);
        poly_item->setName(tr("%1-filled-%2").arg(polylines_item->name()).arg(counter));
        poly_item->setRenderingMode(FlatPlusEdges);
        scene->setSelectedItem(scene->addItem(poly_item));
      }
    }
    QApplication::restoreOverrideCursor();
    }

private:
  QMainWindow* mw;
  CGAL::Three::Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionHoleFillingPolyline;

}; // end Polyhedron_demo_hole_filling_polyline_plugin

// Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_polyline_plugin, Polyhedron_demo_hole_filling_polyline_plugin)

#include "Hole_filling_polyline_plugin.moc"
