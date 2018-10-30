#include <QApplication>
#include <QMainWindow>
#include <QTime>
#include <QAction>
#include <QObject>
#include <QDockWidget>

#include "SMesh_type.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include "ui_Surface_mesh_approximation_dockwidget.h"
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

  typedef VSA_approximation_wrapper<SMesh, EPICK> Approximation_wrapper;
  #ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  typedef Approximation_wrapper::L21_proxy_wrapper L21_proxy_wrapper;
  #endif
  typedef Approximation_wrapper::Indexed_triangle Indexed_triangle;

  typedef EPICK::FT FT;

public:
  Polyhedron_demo_surface_mesh_approximation_plugin() {}

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
    return qobject_cast<Scene_surface_mesh_item *>(scene->item(scene->mainSelectionIndex()));
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

  SMesh *m_pmesh;

  // property-map for segment-idx
  SMesh::Property_map<face_descriptor, std::size_t> m_fidx_pmap;

  // algorithm instance
  Approximation_wrapper m_approx;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::vector<L21_proxy_wrapper> m_proxies;
#endif
  std::vector<std::size_t> m_px_color;
  std::vector<Point_3> m_anchor_pos;
  std::vector<vertex_descriptor> m_anchor_vtx;
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

  // set mesh
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item *>(scene->item(index));
  if (!sm_item) {
    std::cerr << "No surface mesh item selected." << std::endl;
    return;
  }

  ui_widget.seeding->setEnabled(true);
  QApplication::setOverrideCursor(Qt::WaitCursor);

  m_pmesh = sm_item->polyhedron();
  m_fidx_pmap = m_pmesh->add_property_map<face_descriptor, std::size_t>("m_fidx_pmap", 0).first;

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

  std::cout << "is multi color " << sm_item->isItemMulticolor() << std::endl;

  sm_item->setItemIsMulticolor(true);

  typedef typename boost::property_map<SMesh, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
  Patch_id_pmap pidmap = get(CGAL::face_patch_id_t<int>(), *sm_item->face_graph());
  BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
    put(pidmap, f, int(m_fidx_pmap[f]));
  std::vector<QColor> &color_vector = sm_item->color_vector();
  std::cout << "#color_vector " << color_vector.size() << std::endl;
  color_vector.clear();
  BOOST_FOREACH(const std::size_t &cidx, m_px_color)
    color_vector.push_back(QColor::fromRgb(
      Color_cheat_sheet::r(cidx), Color_cheat_sheet::g(cidx), Color_cheat_sheet::b(cidx)));
  // std::vector<QColor> &color_vector = sm_item->color_vector();
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

  sm_item->setItemIsMulticolor(true);
  sm_item->invalidateOpenGLBuffers();
  scene->itemChanged(scene->item_id(sm_item));

  std::cout << "is multi color " << sm_item->isItemMulticolor() << std::endl;
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
    std::vector<Point_3> polyline;
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
    std::vector<Point_3> polyline;
    polyline.push_back(m_pmesh->point(m_anchor_vtx[i]));
    polyline.push_back(m_anchor_pos[i]);
    corres_item->polylines.push_back(polyline);
  }
  corres_item->setName(tr("Anchor vertex correspondence of %1").arg(
    scene->item(scene->mainSelectionIndex())->name()));
  corres_item->setColor(Qt::blue);
  corres_item->invalidateOpenGLBuffers();
  scene->addItem(corres_item);

  // add patch convex hull item
  std::vector<std::vector<EPICK::Triangle_3> > patch_triangles(m_approx.number_of_proxies());\
  BOOST_FOREACH(face_descriptor f, faces(*m_pmesh)) {
    halfedge_descriptor h = m_pmesh->halfedge(f);
    patch_triangles[m_fidx_pmap[f]].push_back(EPICK::Triangle_3(
      m_pmesh->point(source(h, *m_pmesh)),
      m_pmesh->point(target(h, *m_pmesh)),
      m_pmesh->point(target(next(h, *m_pmesh), *m_pmesh))));
  }
  std::vector<EPICK::Plane_3> patch_planes;
  BOOST_FOREACH(const std::vector<EPICK::Triangle_3> &tris, patch_triangles) {
    EPICK::Plane_3 fit_plane;
    CGAL::linear_least_squares_fitting_3(
      tris.begin(), tris.end(), fit_plane, CGAL::Dimension_tag<2>());
    patch_planes.push_back(fit_plane);
  }
  std::vector<std::vector<Point_3> > patch_points(m_approx.number_of_proxies());
  BOOST_FOREACH(vertex_descriptor v, vertices(*m_pmesh)) {
    BOOST_FOREACH(halfedge_descriptor h, CGAL::halfedges_around_target(v, *m_pmesh)) {
      if (!CGAL::is_border(h, *m_pmesh)) {
        const std::size_t fidx = m_fidx_pmap[ face(h, *m_pmesh) ];
        patch_points[fidx].push_back(m_pmesh->point(v));
      }
    }
  }
  std::vector<Point_3> cvx_hull_points;
  std::vector<std::vector<std::size_t> > cvx_hulls;
  for (std::size_t i = 0; i < m_approx.number_of_proxies(); ++i) {
    const std::vector<Point_3> &pts = patch_points[i];
    const EPICK::Plane_3 plane = patch_planes[i];
    const Point_3 origin = plane.projection(pts.front());

    EPICK::Vector_3 base1 = plane.base1();
    EPICK::Vector_3 base2 = plane.base2();
    base1 = base1 / std::sqrt(base1.squared_length());
    base2 = base2 / std::sqrt(base2.squared_length());

    EPICK::Line_3 base_linex(origin, base1);
    EPICK::Line_3 base_liney(origin, base2);

    std::vector<EPICK::Point_2> pts_2d;
    BOOST_FOREACH(const Point_3 &p, pts) {
      const Point_3 point = plane.projection(p);
      EPICK::Vector_3 vecx(origin, base_linex.projection(point));
      EPICK::Vector_3 vecy(origin, base_liney.projection(point));
      double x = std::sqrt(vecx.squared_length());
      double y = std::sqrt(vecy.squared_length());
      x = vecx * base1 < 0 ? -x : x;
      y = vecy * base2 < 0 ? -y : y;
      pts_2d.push_back(EPICK::Point_2(x, y));
    }

    std::vector<EPICK::Point_2> cvx_hull_2d;
    CGAL::convex_hull_2(pts_2d.begin(), pts_2d.end(), std::back_inserter(cvx_hull_2d));

    cvx_hulls.push_back(std::vector<std::size_t>());
    BOOST_FOREACH(const EPICK::Point_2 &p, cvx_hull_2d) {
      cvx_hulls.back().push_back(cvx_hull_points.size());
      cvx_hull_points.push_back(origin + p.x() * base1 + p.y() * base2);
    }
  }
  Scene_polygon_soup_item *cvx_hull_item = new Scene_polygon_soup_item();
  cvx_hull_item->load(cvx_hull_points, cvx_hulls);
  cvx_hull_item->setName(tr("Patch convex hull of %1").arg(
    scene->item(scene->mainSelectionIndex())->name()));
  cvx_hull_item->setColor(Qt::yellow);
  cvx_hull_item->setRenderingMode(FlatPlusEdges);
  scene->addItem(cvx_hull_item);

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
  BOOST_FOREACH(face_descriptor f, faces(*m_pmesh))
    m_fidx_pmap[f] = 0;

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
