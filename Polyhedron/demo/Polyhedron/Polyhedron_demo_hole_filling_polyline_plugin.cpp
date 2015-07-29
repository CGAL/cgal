#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"
#include "Scene_interface.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_type.h"

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Timer.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QInputDialog>
#include <QMessageBox>

#include <vector>
#include <algorithm>

#include <boost/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

template<class HDS>
class Polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
  Polyhedron_builder(std::vector<CGAL::Triple<int, int, int> >* triangles, 
    Scene_polylines_item::Polyline* polyline) 
    : triangles(triangles), polyline(polyline) 
  { }

  void operator()(HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface(polyline->size() -1, triangles->size());

    for(Scene_polylines_item::Polyline::iterator it = polyline->begin();
      it != --polyline->end(); ++it) {
      B.add_vertex(*it);
    }

    for(std::vector<CGAL::Triple<int, int, int> >::iterator it = triangles->begin();
      it != triangles->end(); ++it) {
      B.begin_facet();
      B.add_vertex_to_facet(it->first);
      B.add_vertex_to_facet(it->second);
      B.add_vertex_to_facet(it->third);
      B.end_facet();
    }

    B.end_surface();
  }

private:
  std::vector<CGAL::Triple<int, int, int> >* triangles;
  Scene_polylines_item::Polyline* polyline;
};

class Polyhedron_demo_hole_filling_polyline_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction *) const { return qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionHoleFillingPolyline; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m){
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
      std::vector<CGAL::Triple<int, int, int> > patch;
      CGAL::Polygon_mesh_processing::triangulate_hole_polyline(*it, std::back_inserter(patch), use_DT);
      print_message(QString("Triangulated in %1 sec.").arg(timer.time()));

      if(patch.empty()) {
        print_message("Warning: generating patch is not successful, please try it without 'Delaunay Triangulation'!");
        continue;
      }
      Polyhedron* poly = new Polyhedron;
      Polyhedron_builder<Polyhedron::HalfedgeDS> patch_builder(&patch, &(*it));
      poly->delegate(patch_builder);

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
    }
    QApplication::restoreOverrideCursor();
  }

private:
  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionHoleFillingPolyline;

}; // end Polyhedron_demo_hole_filling_polyline_plugin

// Q_EXPORT_PLUGIN2(Polyhedron_demo_hole_filling_polyline_plugin, Polyhedron_demo_hole_filling_polyline_plugin)

#include "Polyhedron_demo_hole_filling_polyline_plugin.moc"
