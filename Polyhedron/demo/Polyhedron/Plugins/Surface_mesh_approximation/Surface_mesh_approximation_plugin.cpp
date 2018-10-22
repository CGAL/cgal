#include <QApplication>
#include <QMainWindow>
#include <QTime>
#include <QAction>
#include <QObject>
#include <QDockWidget>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include "ui_Surface_mesh_approximation_dockwidget.h"
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include "VSA_approximation_wrapper.h"

using namespace CGAL::Three;

typedef VSA_approximation_wrapper<Polyhedron, Kernel> Approximation_wrapper;
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
typedef Approximation_wrapper::L21_proxy_wrapper L21_proxy_wrapper;
#endif
typedef Approximation_wrapper::Indexed_triangle Indexed_triangle;

typedef Kernel::FT FT;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Vertex_handle Vertex_handle;

class Polyhedron_demo_surface_mesh_approximation_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow *main_window,
    Scene_interface *scene_interface,
    Messages_interface *message_interface) {
    mw = main_window;
    scene = scene_interface;
    mi = message_interface;

    QAction *actionSurfaceMeshApproximation = new QAction("Surface Mesh Approximation", mw);
    actionSurfaceMeshApproximation->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionSurfaceMeshApproximation, SIGNAL(triggered()), this, SLOT(on_actionSurfaceMeshApproximation_triggered()));

    // set metric menu
    // QAction *actionL21 = new QAction("L21 metric", mw);
    actionL21 = new QAction("L21 metric", mw);
    actionL21->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionL21, SIGNAL(triggered()), this, SLOT(on_actionL21_triggered()));

    // QAction *actionL2 = new QAction("L2 metric", mw);
    actionL2 = new QAction("L2 metric", mw);
    actionL2->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionL2, SIGNAL(triggered()), this, SLOT(on_actionL2_triggered()));

    // QAction *actionCompact = new QAction("Compact metric", mw);
    actionCompact = new QAction("Compact metric", mw);
    actionCompact->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionCompact, SIGNAL(triggered()), this, SLOT(on_actionL21_triggered()));

    // operations menu
    // QAction *actionApproximation = new QAction("Approximation", mw);
    actionApproximation = new QAction("Approximation", mw);
    actionApproximation->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionApproximation, SIGNAL(triggered()), this, SLOT(on_actionApproximation_triggered()));

    // QAction *actionSeeding = new QAction("Seeding", mw);
    actionSeeding = new QAction("Seeding", mw);
    actionSeeding->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionSeeding, SIGNAL(triggered()), this, SLOT(on_actionSeeding_triggered()));

    // QAction *actionFit = new QAction("Fit", mw);
    actionFit = new QAction("Fit", mw);
    actionFit->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionFit, SIGNAL(triggered()), this, SLOT(on_actionFit_triggered()));

    // QAction *actionMeshing = new QAction("Meshing", mw);
    actionMeshing = new QAction("Meshing", mw);
    actionMeshing->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionMeshing, SIGNAL(triggered()), this, SLOT(on_actionMeshing_triggered()));

    // QAction *actionAdd = new QAction("Add proxy", mw);
    actionAdd = new QAction("Add proxy", mw);
    actionAdd->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionAdd, SIGNAL(triggered()), this, SLOT(on_actionAdd_triggered()));

    // QAction *actionTeleport = new QAction("Teleport proxy", mw);
    actionTeleport = new QAction("Teleport proxy", mw);
    actionTeleport->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionTeleport, SIGNAL(triggered()), this, SLOT(on_actionTeleport_triggered()));

    // QAction *actionSplit = new QAction("Split proxy", mw);
    actionSplit = new QAction("Split proxy", mw);
    actionSplit->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionSplit, SIGNAL(triggered()), this, SLOT(on_actionSplit_triggered()));

    // view menu
    // QAction *actionViewPolyhedron = new QAction("View polyhedron", mw);
    actionViewPolyhedron = new QAction("View polyhedron", mw);
    actionViewPolyhedron->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionViewPolyhedron, SIGNAL(triggered()), this, SLOT(on_actionViewPolyhedron_triggered()));

    _actions << actionSurfaceMeshApproximation
             // << NULL
             << actionL21
             << actionL2
             << actionCompact
             // << NULL
             << actionApproximation
             << actionSeeding
             << actionFit
             << actionMeshing
             << actionAdd
             << actionTeleport
             << actionSplit
             // << NULL
             << actionViewPolyhedron;

    dock_widget = new QDockWidget("Mesh approximation parameters", mw);
    dock_widget->setVisible(true);
    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);
  }

  void closure() {
    // dock_widget->hide();
  }

  QList<QAction *> actions() const { return _actions; }

  bool applicable(QAction *) const {
    return 
      qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionSurfaceMeshApproximation_triggered();

  // settings
  // void quit();
  // void readSettings();
  // void writeSettings();

  // set metric menu
  void on_actionL21_triggered();
  void on_actionL2_triggered();
  void on_actionCompact_triggered();

  // operations menu
  void on_actionApproximation_triggered();
  void on_actionSeeding_triggered();
  void on_actionFit_triggered();
  void on_actionMeshing_triggered();
  void on_actionAdd_triggered();
  void on_actionTeleport_triggered();
  void on_actionSplit_triggered();

  // view menu
  void on_actionViewPolyhedron_triggered() {
    m_view_polyhedron = !m_view_polyhedron;
  }
  void on_actionViewWireframe_triggered() {
    m_view_wireframe = !m_view_wireframe;
  }
  void on_actionViewBoundary_triggered() {
    m_view_boundary = !m_view_boundary;
  }
  void on_actionViewProxies_triggered() {
    m_view_proxies = !m_view_proxies;
  }
  void on_actionViewAnchors_triggered() {
    m_view_anchors = !m_view_anchors;
  }
  void on_actionViewApproximation_triggered() {
    m_view_approximation = !m_view_approximation;
  }

private:
  // set metric algorithms
  void set_metric(const int init);

  // pseudorandom number for proxy color mapping
  std::size_t rand_0_255() {
    return static_cast<std::size_t>(std::rand() % 255);
  }

private:
  QAction *actionL21;
  QAction *actionL2;
  QAction *actionCompact;
  QAction *actionApproximation;
  QAction *actionSeeding;
  QAction *actionFit;
  QAction *actionMeshing;
  QAction *actionAdd;
  QAction *actionTeleport;
  QAction *actionSplit;
  QAction *actionViewPolyhedron;
  QAction *actionViewWireframe;
  QAction *actionViewBoundary;
  QAction *actionViewProxies;
  QAction *actionViewAnchors;
  QAction *actionViewApproximation;
  QList<QAction *> _actions;

  Ui::Surface_mesh_approximation ui_widget;
  QDockWidget *dock_widget;

  QMainWindow *mw;
  Scene_interface *scene;
  Messages_interface *mi;

  Polyhedron *m_pmesh;

  // property-map for segment-idx
  std::map<Facet_handle, std::size_t> m_fidx_map;
  boost::associative_property_map<std::map<Facet_handle, std::size_t> > m_fidx_pmap;

  // algorithm instance
  Approximation_wrapper m_approx;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::vector<L21_proxy_wrapper> m_proxies;
#endif
  std::vector<std::size_t> m_px_color;
  std::vector<Point_3> m_anchor_pos;
  std::vector<Polyhedron::Vertex_handle> m_anchor_vtx;
  std::vector<std::vector<std::size_t> > m_bdrs; // anchor borders
  std::vector<Indexed_triangle> m_tris;

  // view options
  bool m_view_polyhedron;
  bool m_view_wireframe;
  bool m_view_boundary;
  bool m_view_proxies;
  bool m_view_anchors;
  bool m_view_approximation;
}; // end Polyhedron_demo_surface_mesh_approximation_plugin

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionSurfaceMeshApproximation_triggered()
{
  dock_widget->show();
  return;

  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item = 
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if (poly_item || sm_item) {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    
    QTime time;
    time.start();
    std::cout << "Surface mesh approximation...";

    // add surface mesh approximation output as new polyhedron
    SMesh *pResult  = new SMesh;
    if (poly_item) {
      Polyhedron *pMesh = poly_item->polyhedron();
      VSA::approximate_triangle_mesh(*pMesh, CGAL::parameters::max_number_of_proxies(30).
        number_of_iterations(20).
        subdivision_ratio(3).
        relative_to_chord(false));
    }
    else {
      SMesh *pMesh = sm_item->polyhedron();
      VSA::approximate_triangle_mesh(*pMesh, CGAL::parameters::max_number_of_proxies(30).
        number_of_iterations(20).
        subdivision_ratio(3).
        relative_to_chord(false));
    }
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    if(mw->property("is_polyhedron_mode").toBool()){
      Polyhedron *poly = new Polyhedron;
      CGAL::copy_face_graph(*pResult,*poly);
      delete pResult;

      Scene_polyhedron_item* new_item = new Scene_polyhedron_item(poly);
      new_item->setName(tr("%1 (surface mesh approximation)").arg(scene->item(index)->name()));
      new_item->setColor(Qt::magenta);
      new_item->setRenderingMode(FlatPlusEdges);
      scene->addItem(new_item);
    } else {
       Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(pResult);
       new_item->setName(tr("%1 (surface mesh approximation)").arg(scene->item(index)->name()));
       new_item->setColor(Qt::magenta);
       new_item->setRenderingMode(FlatPlusEdges);
       scene->addItem(new_item);
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionL21_triggered() {
  actionL2->setChecked(false);
  actionCompact->setChecked(false);

  actionViewPolyhedron->setChecked(true);
  actionViewWireframe->setChecked(false);
  actionViewBoundary->setChecked(false);
  actionViewProxies->setChecked(false);
  actionViewAnchors->setChecked(false);
  actionViewApproximation->setChecked(false);

  set_metric(0);
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionL2_triggered() {
  actionL21->setChecked(false);
  actionCompact->setChecked(false);

  actionViewPolyhedron->setChecked(true);
  actionViewWireframe->setChecked(false);
  actionViewBoundary->setChecked(false);
  actionViewProxies->setChecked(false);
  actionViewAnchors->setChecked(false);
  actionViewApproximation->setChecked(false);

  set_metric(1);
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionCompact_triggered() {
  actionL21->setChecked(false);
  actionL2->setChecked(false);

  actionViewPolyhedron->setChecked(true);
  actionViewWireframe->setChecked(false);
  actionViewBoundary->setChecked(false);
  actionViewProxies->setChecked(false);
  actionViewAnchors->setChecked(false);
  actionViewApproximation->setChecked(false);

  set_metric(2);
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionApproximation_triggered() {
  ui_widget.seeding->setEnabled(true);
  ui_widget.mesh_extraction->setEnabled(true);
  QApplication::setOverrideCursor(Qt::WaitCursor);

  const VSA::Seeding_method method = ui_widget.method_random->isChecked() ? VSA::RANDOM : (
    ui_widget.method_incremental->isChecked() ? VSA::INCREMENTAL : VSA::HIERARCHICAL);
  m_approx.initialize_seeds(method,
    (ui_widget.cb_nb_proxies->isChecked() ? boost::optional<std::size_t>(ui_widget.nb_proxies->value()) : boost::none),
    (ui_widget.cb_error_drop->isChecked() ? boost::optional<FT>(ui_widget.error_drop->value()) : boost::none),
    ui_widget.nb_relaxations->value());
  m_approx.run(ui_widget.nb_iterations->value());

  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif
  // generate proxy color map
  m_px_color.clear();
  for (std::size_t i = 0; i < m_approx.number_of_proxies(); i++)
    m_px_color.push_back(rand_0_255());

  m_tris.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();
  m_approx.extract_mesh(ui_widget.chord_error->value(),
    ui_widget.is_relative_to_chord->isChecked(),
    ui_widget.with_dihedral_angle->isChecked(),
    ui_widget.if_optimize_anchor_location->isChecked(),
    ui_widget.pca_plane->isChecked());
  m_approx.indexed_triangles(std::back_inserter(m_tris));
  m_approx.anchor_points(std::back_inserter(m_anchor_pos));
  m_approx.anchor_vertices(std::back_inserter(m_anchor_vtx));
  m_approx.indexed_boundary_polygons(std::back_inserter(m_bdrs));

  // update display options
  m_view_boundary = true;
  m_view_anchors = true;
  m_view_approximation = true;
  actionViewBoundary->setChecked(true);
  actionViewAnchors->setChecked(true);
  actionViewApproximation->setChecked(true);

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionSeeding_triggered() {
  ui_widget.seeding->setEnabled(true);
  QApplication::setOverrideCursor(Qt::WaitCursor);

  const VSA::Seeding_method method = ui_widget.method_random->isChecked() ? VSA::RANDOM : (
    ui_widget.method_incremental->isChecked() ? VSA::INCREMENTAL : VSA::HIERARCHICAL);
  m_approx.initialize_seeds(method,
    (ui_widget.cb_nb_proxies->isChecked() ? boost::optional<std::size_t>(ui_widget.nb_proxies->value()) : boost::none),
    (ui_widget.cb_error_drop->isChecked() ? boost::optional<FT>(ui_widget.error_drop->value()) : boost::none),
    ui_widget.nb_relaxations->value());
  m_approx.run(ui_widget.nb_iterations->value());

  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif
  // generate proxy color map
  m_px_color.clear();
  for (std::size_t i = 0; i < m_approx.number_of_proxies(); i++)
    m_px_color.push_back(rand_0_255());

  m_view_polyhedron = true;
  actionViewBoundary->setChecked(true);

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionFit_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);

  m_approx.run(1);
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif
  
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionMeshing_triggered() {
  ui_widget.mesh_extraction->setEnabled(true);
  QApplication::setOverrideCursor(Qt::WaitCursor);

  m_tris.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();

  m_approx.extract_mesh(ui_widget.chord_error->value(),
    ui_widget.is_relative_to_chord->isChecked(),
    ui_widget.with_dihedral_angle->isChecked(),
    ui_widget.if_optimize_anchor_location->isChecked(),
    ui_widget.pca_plane->isChecked());

  m_approx.indexed_triangles(std::back_inserter(m_tris));
  m_approx.anchor_points(std::back_inserter(m_anchor_pos));
  m_approx.anchor_vertices(std::back_inserter(m_anchor_vtx));
  m_approx.indexed_boundary_polygons(std::back_inserter(m_bdrs));

  m_view_anchors = true;
  m_view_approximation = true;
  actionViewAnchors->setChecked(true);
  actionViewApproximation->setChecked(true);

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionAdd_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (m_approx.add_one_proxy() == 0)
    return;

  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  // add one proxy color
  m_px_color.push_back(rand_0_255());

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionTeleport_triggered() {
  QApplication::setOverrideCursor(Qt::WaitCursor);

  m_approx.teleport_one_proxy();
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionSplit_triggered() {
  ui_widget.operations->setEnabled(true);
  ui_widget.operations->setCurrentIndex(0);
  QApplication::setOverrideCursor(Qt::WaitCursor);

  if (m_approx.split(ui_widget.split_proxy_idx->value(),
    ui_widget.split_nb_sections->value(),
    ui_widget.split_nb_relaxations->value())) {
    // add colors
    for (std::size_t i = 1; i < ui_widget.split_nb_sections->value(); ++i)
      m_px_color.push_back(rand_0_255());
  }
  else
    std::cout << "split failed" << std::endl;

  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::set_metric(const int m) {
  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    m_fidx_map[fitr] = 0;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
#endif
  m_px_color.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();
  m_tris.clear();

  m_view_polyhedron = true;
  m_view_wireframe = false;
  m_view_boundary = false;
  m_view_proxies = false;
  m_view_anchors = false;
  m_view_approximation = false;

  switch (m) {
    case 0: return m_approx.set_metric(Approximation_wrapper::L21);
    case 1: return m_approx.set_metric(Approximation_wrapper::L2);
    case 2: return m_approx.set_metric(Approximation_wrapper::Compact);
  }
}

#include "Surface_mesh_approximation_plugin.moc"
