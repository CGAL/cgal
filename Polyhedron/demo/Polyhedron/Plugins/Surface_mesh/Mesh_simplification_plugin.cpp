#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QElapsedTimer>
#include <QAction>
#include <QMessageBox>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>

#include "ui_Mesh_simplification_dialog.h"

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;

class Custom_stop_predicate
{
  bool m_and;
  CGAL::Surface_mesh_simplification::Count_stop_predicate<FaceGraph> m_count_stop;
  CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double> m_length_stop;

public:

  Custom_stop_predicate (bool use_and, std::size_t nb_edges, double edge_length)
    : m_and (use_and), m_count_stop (nb_edges), m_length_stop (edge_length)
  {
    std::cerr << "Simplifying until:" << std::endl
              << " * Number of edges = " << nb_edges << std::endl
              << (use_and ? " AND " : " OR ") << std::endl
              << " * Minimum edge length = " << edge_length << std::endl;
  }

  template <typename Profile>
  bool operator() (const double& current_cost, const Profile& edge_profile,
                   std::size_t initial_count, std::size_t current_count) const
  {
    if (m_and)
      return (m_count_stop(current_cost, edge_profile, initial_count, current_count)
              && m_length_stop(current_cost, edge_profile, initial_count, current_count));
    else
      return (m_count_stop(current_cost, edge_profile, initial_count, current_count)
              || m_length_stop(current_cost, edge_profile, initial_count, current_count));
  }

};


using namespace CGAL::Three;
class Polyhedron_demo_mesh_simplification_plugin :
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

  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
      mw = mainWindow;
      scene = scene_interface;
     QAction *actionSimplify = new QAction("Simplification", mw);
      actionSimplify->setProperty("subMenuName",
                                  "Triangulated Surface Mesh Simplification");
      connect(actionSimplify, SIGNAL(triggered()), this, SLOT(on_actionSimplify_triggered()));
      _actions <<actionSimplify;

  }
  bool applicable(QAction*) const {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()))
      || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
  }
public Q_SLOTS:
  void on_actionSimplify_triggered();
private :
  Scene_interface *scene;
  QMainWindow *mw;
  QList<QAction*> _actions;

}; // end Polyhedron_demo_mesh_simplification_plugin

void Polyhedron_demo_mesh_simplification_plugin::on_actionSimplify_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_facegraph_item* poly_item =
    qobject_cast<Scene_facegraph_item*>(scene->item(index));

  Scene_polyhedron_selection_item* selection_item =
    qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
  if(selection_item && selection_item->selected_edges.empty())
  {
    QMessageBox::warning(mw, "Empty Edges", "There are no selected edges. Aborting.");
    return;
  }
  if (poly_item || selection_item)
  {
    FaceGraph& pmesh = (poly_item != NULL)
      ? *poly_item->polyhedron()
      : *selection_item->polyhedron();

    // get option
    QDialog dialog(mw);
    Ui::Mesh_simplification_dialog ui;
    ui.setupUi(&dialog);
    connect(ui.buttonBox, SIGNAL(accepted()),
            &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()),
            &dialog, SLOT(reject()));

    Scene_interface::Bbox bbox = poly_item != NULL ? poly_item->bbox()
      : (selection_item != NULL ? selection_item->bbox()
        : scene->bbox());

    double diago_length = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
                                     + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin()) +
                                     (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));

    ui.m_nb_edges->setValue ((int)(num_halfedges(pmesh) / 4));
    ui.m_nb_edges->setMaximum ((int)(num_halfedges(pmesh)));
    ui.m_edge_length->setValue (diago_length * 0.05);

    // check user cancellation
    if(dialog.exec() == QDialog::Rejected)
      return;

    // simplify
    QElapsedTimer time;
    time.start();
    std::cout << "Simplify...";
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    Custom_stop_predicate stop ((ui.m_combinatorial->currentIndex() == 0)
                                && !(ui.m_use_nb_edges->isChecked())
                                && !(ui.m_use_edge_length->isChecked()),
                                (ui.m_use_nb_edges->isChecked()
                                 ? ui.m_nb_edges->value()
                                 : 0),
                                (ui.m_use_edge_length->isChecked()
                                 ? ui.m_edge_length->value()
                                 : (std::numeric_limits<double>::max)()));

    if (selection_item)
      {
        CGAL::Surface_mesh_simplification::Constrained_placement
          <CGAL::Surface_mesh_simplification::Bounded_normal_change_placement
           <CGAL::Surface_mesh_simplification::LindstromTurk_placement
            <FaceGraph> >,
           Scene_polyhedron_selection_item::Is_constrained_map
           <Scene_polyhedron_selection_item::Selection_set_edge> >
          placement (selection_item->constrained_edges_pmap());

        CGAL::Surface_mesh_simplification::edge_collapse
          (pmesh, stop,
           CGAL::parameters::edge_is_constrained_map(selection_item->constrained_edges_pmap())
           .get_placement(placement));
      }
    else
      {
        CGAL::Surface_mesh_simplification::Bounded_normal_change_placement
          <CGAL::Surface_mesh_simplification::LindstromTurk_placement
           <FaceGraph> > placement;

        CGAL::Surface_mesh_simplification::edge_collapse
          (pmesh, stop,
           CGAL::parameters::vertex_index_map(get(boost::vertex_index, pmesh)).get_placement(placement));
      }

    std::cout << "ok (" << time.elapsed() << " ms, "
      << num_halfedges(pmesh) / 2 << " edges)" << std::endl;

    // update scene
    if (poly_item != NULL)
    {
      poly_item->invalidateOpenGLBuffers();
      poly_item->polyhedron()->collect_garbage();
    }
    else
      {
      selection_item->polyhedron_item()->polyhedron()->collect_garbage();
      selection_item->poly_item_changed();
      selection_item->changed_with_poly_item();
      selection_item->invalidateOpenGLBuffers();
      }

    scene->itemChanged(index);
    QApplication::restoreOverrideCursor();
  }
}

#include "Mesh_simplification_plugin.moc"
