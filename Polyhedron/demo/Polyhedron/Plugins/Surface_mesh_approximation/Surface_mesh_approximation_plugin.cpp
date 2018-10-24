#include <QApplication>
#include <QMainWindow>
#include <QTime>
#include <QAction>
#include <QObject>
#include <QDockWidget>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"
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
#include "Color_cheat_sheet.h"

using namespace CGAL::Three;

class Polyhedron_demo_surface_mesh_approximation_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  typedef VSA_approximation_wrapper<Polyhedron, Kernel> Approximation_wrapper;
  #ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  typedef Approximation_wrapper::L21_proxy_wrapper L21_proxy_wrapper;
  #endif
  typedef Approximation_wrapper::Indexed_triangle Indexed_triangle;

  // typedef Kernel::Point_3 Point_3;
  typedef Kernel::FT FT;
  typedef Polyhedron::Facet_iterator Facet_iterator;
  typedef Polyhedron::Facet_handle Facet_handle;
  typedef Polyhedron::Vertex_handle Vertex_handle;

public:
  Polyhedron_demo_surface_mesh_approximation_plugin() :
    m_fidx_pmap(m_fidx_map) {}

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
      qobject_cast<Scene_polyhedron_item *>(scene->item(scene->mainSelectionIndex()))
        || qobject_cast<Scene_surface_mesh_item *>(scene->item(scene->mainSelectionIndex()));
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
  std::vector<Polyhedron::Vertex::Point> m_anchor_pos;
  std::vector<Polyhedron::Vertex_handle> m_anchor_vtx;
  std::vector<std::vector<std::size_t> > m_bdrs; // anchor borders
  std::vector<Indexed_triangle> m_tris;
}; // end Polyhedron_demo_surface_mesh_approximation_plugin

void Polyhedron_demo_surface_mesh_approximation_plugin::on_actionSurfaceMeshApproximation_triggered()
{
  dock_widget->show();
  return;
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

  // set mesh
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item *poly_item = qobject_cast<Scene_polyhedron_item *>(scene->item(index));
  if (!poly_item) {
    std::cerr << "The selected should be a polyhedron item" << std::endl;
    return;
  }

  m_pmesh = poly_item->polyhedron();
  m_fidx_map.clear();
  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    m_fidx_map.insert(std::pair<Facet_handle, std::size_t>(fitr, 0));

  m_approx.set_mesh(*m_pmesh);
  m_approx.set_metric(Approximation_wrapper::L21);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
#endif
  m_px_color.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();
  m_tris.clear();

  QTime time;
  time.start();
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

  std::cout << "========" << std::endl;

  // generate proxy color map
  for (std::size_t i = 0; i < m_approx.number_of_proxies(); i++)
    m_px_color.push_back(rand_0_255());

  std::cout << "#proxies " << m_approx.number_of_proxies() << std::endl;
  std::cout << "Done. (" << time.elapsed() << " ms)" << std::endl;

  // face color
  m_approx.proxy_map(m_fidx_pmap);

  // set face ids
  std::size_t face_id = 0;
  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    fitr->id() = face_id++;

  std::cout << "is multi color " << poly_item->isItemMulticolor() << std::endl;

  poly_item->set_color_vector_read_only(true);

  typedef typename boost::property_map<Polyhedron, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
  Patch_id_pmap pidmap = get(CGAL::face_patch_id_t<int>(), *poly_item->face_graph());
  for(Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    put(pidmap, fitr, int(m_fidx_pmap[fitr]));
  std::vector<QColor> &color_vector = poly_item->color_vector();
  std::cout << "#color_vector " << color_vector.size() << std::endl;
  color_vector.clear();
  BOOST_FOREACH(const std::size_t &cidx, m_px_color)
    color_vector.push_back(QColor::fromRgb(
      Color_cheat_sheet::r(cidx), Color_cheat_sheet::g(cidx), Color_cheat_sheet::b(cidx)));
  // std::vector<QColor> &color_vector = poly_item->color_vector();
  // std::cout << "#color_vector " << color_vector.size() << std::endl;
  // color_vector.clear();
  // for(Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr) {
  //   put(pidmap, fitr, int(m_fidx_pmap[fitr]));
  //   const std::size_t cidx = m_px_color[m_fidx_pmap[fitr]];
  //   color_vector.push_back(QColor::fromRgb(
  //     Color_cheat_sheet::r(cidx), Color_cheat_sheet::g(cidx), Color_cheat_sheet::b(cidx)));
  // }
  std::cout << "color done." << std::endl;
  std::cout << "#color_vector " << color_vector.size() << std::endl;

  poly_item->setItemIsMulticolor(true);
  poly_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(poly_item));

  std::cout << "is multi color " << poly_item->isItemMulticolor() << std::endl;
  std::cout << "========" << std::endl;

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

  // add triangle soup item
  Scene_polygon_soup_item *psoup_item = new Scene_polygon_soup_item();
  std::vector<std::vector<std::size_t> > polygons;
  BOOST_FOREACH(const Indexed_triangle &t, m_tris) {
    std::vector<std::size_t> polygon;
    polygon.push_back(t[0]);
    polygon.push_back(t[1]);
    polygon.push_back(t[2]);
    polygons.push_back(polygon);
  }
  psoup_item->load(m_anchor_pos, polygons);
  psoup_item->setName(tr("Approximated triangle soup of %1").arg(
    scene->item(scene->mainSelectionIndex())->name()));
  psoup_item->setColor(Qt::lightGray);
  psoup_item->setRenderingMode(FlatPlusEdges);
  scene->addItem(psoup_item);

  // add patch border item
  Scene_polylines_item *borders_item = new Scene_polylines_item();
  BOOST_FOREACH(const std::vector<std::size_t> &border, m_bdrs) {
    std::vector<Kernel::Point_3> polyline;
    for (std::size_t i = 0; i <= border.size(); ++i)
      polyline.push_back(m_anchor_pos[border[i % border.size()]]);
    borders_item->polylines.push_back(polyline);
  }
  borders_item->setName(tr("Patch boundary of %1").arg(
    scene->item(scene->mainSelectionIndex())->name()));
  borders_item->setColor(Qt::red);
  borders_item->invalidateOpenGLBuffers();
  scene->addItem(borders_item);

  // add anchor vertex and anchor point correspondence item
  Scene_polylines_item *corres_item = new Scene_polylines_item();
  for (std::size_t i = 0; i < m_anchor_vtx.size(); ++i) {
    std::vector<Kernel::Point_3> polyline;
    polyline.push_back(m_anchor_vtx[i]->point());
    polyline.push_back(m_anchor_pos[i]);
    corres_item->polylines.push_back(polyline);
  }
  corres_item->setName(tr("Anchor vertex correspondence of %1").arg(
    scene->item(scene->mainSelectionIndex())->name()));
  corres_item->setColor(Qt::blue);
  corres_item->invalidateOpenGLBuffers();
  scene->addItem(corres_item);

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
