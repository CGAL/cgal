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

    actionSurfaceMeshApproximation = new QAction("Surface Mesh Approximation", mw);
    actionSurfaceMeshApproximation->setProperty("subMenuName", "Triangulated Surface Mesh Approximation");
    connect(actionSurfaceMeshApproximation, SIGNAL(triggered()), this, SLOT(on_actionSurfaceMeshApproximation_triggered()));

    dock_widget = new QDockWidget("Mesh approximation parameters", mw);
    dock_widget->setVisible(true);
    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    // connect ui actions
    connect(ui_widget.comboMetric, SIGNAL(currentIndexChanged(int)), this, SLOT(on_comboMetric_currentIndexChanged(const int)));
    connect(ui_widget.buttonSeeding, SIGNAL(clicked()), this, SLOT(on_buttonSeeding_clicked()));
    connect(ui_widget.buttonAdd, SIGNAL(clicked()), this, SLOT(on_buttonAdd_clicked()));
    connect(ui_widget.buttonFit, SIGNAL(clicked()), this, SLOT(on_buttonFit_clicked()));
    connect(ui_widget.buttonMeshing, SIGNAL(clicked()), this, SLOT(on_buttonMeshing_clicked()));
    connect(ui_widget.buttonTeleport, SIGNAL(clicked()), this, SLOT(on_buttonTeleport_clicked()));
    connect(ui_widget.buttonSplit, SIGNAL(clicked()), this, SLOT(on_buttonSplit_clicked()));
    // connect(actionApproximation, SIGNAL(clicked()), this, SLOT(on_actionApproximation_clicked()));
  }

  void closure() {
    // dock_widget->hide();
  }

  QList<QAction *> actions() const {
    return QList<QAction*>() << actionSurfaceMeshApproximation;
  }

  bool applicable(QAction *) const {
    return 
      qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionSurfaceMeshApproximation_triggered();

  void on_comboMetric_currentIndexChanged(const int init);
  void on_buttonSeeding_clicked();
  void on_buttonAdd_clicked();
  void on_buttonFit_clicked();
  void on_buttonMeshing_clicked();
  void on_buttonTeleport_clicked();
  void on_buttonSplit_clicked();
  void on_actionApproximation_clicked();

  // settings
  // void quit();
  // void readSettings();
  // void writeSettings();

private:
  // pseudorandom number for proxy color mapping
  std::size_t rand_0_255() {
    return static_cast<std::size_t>(std::rand() % 255);
  }

private:
  QAction *actionSurfaceMeshApproximation;
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

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionApproximation_clicked() {
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

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonSeeding_clicked() {
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

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonFit_clicked() {
  QApplication::setOverrideCursor(Qt::WaitCursor);

  m_approx.run(1);
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonMeshing_clicked() {
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

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonAdd_clicked() {
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

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonTeleport_clicked() {
  QApplication::setOverrideCursor(Qt::WaitCursor);

  m_approx.teleport_one_proxy();
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_surface_mesh_approximation_plugin::on_buttonSplit_clicked() {
  ui_widget.operations->setEnabled(true);
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

void Polyhedron_demo_surface_mesh_approximation_plugin::on_comboMetric_currentIndexChanged(const int m) {
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

  switch (m) {
    case 0: return m_approx.set_metric(Approximation_wrapper::L21);
    case 1: return m_approx.set_metric(Approximation_wrapper::L2);
    case 2: return m_approx.set_metric(Approximation_wrapper::Compact);
  }
}

#include "Surface_mesh_approximation_plugin.moc"
