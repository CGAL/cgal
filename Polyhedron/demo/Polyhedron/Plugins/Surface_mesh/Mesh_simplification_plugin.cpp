#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

#include "ui_Mesh_simplification_dialog.h"

class Custom_stop_predicate
{
  bool m_and;
    CGAL::Surface_mesh_simplification::Count_stop_predicate<Polyhedron> m_count_stop;
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
      return (m_count_stop(current_count, edge_profile, initial_count, current_count)
              && m_length_stop(current_count, edge_profile, initial_count, current_count));
    else
      return (m_count_stop(current_count, edge_profile, initial_count, current_count)
              || m_length_stop(current_count, edge_profile, initial_count, current_count));
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
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
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
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    // get option
    QDialog dialog(mw);
    Ui::Mesh_simplification_dialog ui;
    ui.setupUi(&dialog);
    connect(ui.buttonBox, SIGNAL(accepted()),
            &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()),
            &dialog, SLOT(reject()));

    Scene_interface::Bbox bbox = item->bbox();
    double diago_length = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
                                     + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin()) +
                                     (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
    
    ui.m_nb_edges->setValue ((int)(pMesh->size_of_halfedges () / 4));
    ui.m_nb_edges->setMaximum ((int)(pMesh->size_of_halfedges ()));
    ui.m_edge_length->setValue (diago_length * 0.05);

    // check user cancellation
    if(dialog.exec() == QDialog::Rejected)
      return;

    // simplify
    QTime time;
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
                                 : std::numeric_limits<double>::max()));
    
    CGAL::Surface_mesh_simplification::edge_collapse( *pMesh, stop,
                        CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,*pMesh))
                                         .halfedge_index_map(get(CGAL::halfedge_external_index,*pMesh)));
    std::cout << "ok (" << time.elapsed() << " ms, " 
      << pMesh->size_of_halfedges() / 2 << " edges)" << std::endl;

    // update scene
    item->invalidateOpenGLBuffers();
    scene->itemChanged(index);
    QApplication::restoreOverrideCursor();
  }
}

#include "Mesh_simplification_plugin.moc"
