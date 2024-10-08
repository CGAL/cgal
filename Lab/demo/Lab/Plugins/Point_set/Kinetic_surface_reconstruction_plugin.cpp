
#include "ui_Kinetic_surface_reconstruction_plugin.h"

#include "Color_map.h"
#include "Color_ramp.h"
#include "id_printing.h"
#include "Messages_interface.h"

#include "Scene.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Kinetic_surface_reconstruction_3.h>

#include <QAbstractItemView>
#include <QAction>
#include <QApplication>
#include <QColor>
#include <QColorDialog>
#include <QInputDialog>
#include <QMainWindow>
#include <QMessageBox>
#include <QSlider>
#include <QObject>
#include <QPalette>
#include <QStyleFactory>

#include <boost/algorithm/clamp.hpp>
#include <boost/range/value_type.hpp>

#include <type_traits>
#include <unordered_map>
#include <vector>


#include <typeinfo>

using namespace CGAL::Three;

using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using Segment_3 = typename Kernel::Segment_3;

using Point_map = typename Point_set::Point_map;
using Normal_map = typename Point_set::Vector_map;

using KSR = CGAL::Kinetic_surface_reconstruction_3<CGAL::Epick, Point_set, Point_map, Normal_map>;

Viewer_interface* (&getActiveViewer)() = Three::activeViewer;

class DockWidget
  : public QDockWidget,
  public Ui::KineticSurfaceReconstructionWidget
{
public:
  DockWidget(const QString& name, QWidget* parent)
    : QDockWidget(name, parent)
  {
    setupUi(this);
  }
};

class Kinetic_surface_reconstruction_plugin
  : public QObject,
  public CGAL_Lab_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")
private:
  QAction* actionKineticSurfaceReconstruction;

  DockWidget* dock_widget;

  Scene_points_with_normal_item* m_pwn_item = nullptr;

  // A dock widget allows switching between items which makes managing parameters and intermediate results a mess.
  // I could create a map from item to KSR and also reference intermediate results (in case they are recalculated with changed parameters)?

  KSR *m_ksr = nullptr;
  bool m_known_file = false;

public:
  virtual ~Kinetic_surface_reconstruction_plugin() {
    if (m_ksr)
      delete m_ksr;
    m_ksr = nullptr;
  }

  bool applicable(QAction*) const override
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if (!item)
      return false;

    Scene_points_with_normal_item* pwn_item = qobject_cast<Scene_points_with_normal_item*>(item);

    if (pwn_item != nullptr && pwn_item->has_normals())
      return true;
    else
      return false;
  }

  QList<QAction*> actions() const override
  {
    return QList<QAction*>() << actionKineticSurfaceReconstruction;
  }

  void init(QMainWindow* mw,
    Scene_interface* sc,
    Messages_interface*) override
  {
    this->scene = sc;
    this->mw = mw;

    // Main action
    actionKineticSurfaceReconstruction = new QAction(QString("Kinetic Surface Reconstruction"), mw);
    actionKineticSurfaceReconstruction->setObjectName("actionKineticSurfaceReconstruction");

    connect(actionKineticSurfaceReconstruction, SIGNAL(triggered()),
      this, SLOT(openDialog()));

    Scene* scene_item = static_cast<Scene*>(scene);
    connect(scene_item, SIGNAL(itemIndexSelected(int)),
      this, SLOT(onItemIndexSelected(int)));

    // Dock Widget
    dock_widget = new DockWidget("Kinetic Surface Reconstruction", mw);
    addDockWidget(dock_widget);

    dock_widget->setVisible(false);
    //dock_widget->setEnabled(false);

    connect(dock_widget->sdRunButton, SIGNAL(clicked(bool)), this, SLOT(run_detection()));
    connect(dock_widget->partRunButton, SIGNAL(clicked(bool)), this, SLOT(run_partition()));
    connect(dock_widget->recRunButton, SIGNAL(clicked(bool)), this, SLOT(run_reconstruction()));
    connect(dock_widget->partSubdivisionCheck, SIGNAL(stateChanged(int)), this, SLOT(onSubdivisionStateChanged(int)));
    connect(dock_widget, SIGNAL(visibilityChanged(bool)), this, SLOT(onVisibilityChanged(bool)));
  }
private Q_SLOTS:
  void openDialog()
  {
    if (!dock_widget->isVisible())
      dock_widget->show();
    dock_widget->raise();
  }

  void closure() override
  {
    dock_widget->hide();
    if (m_ksr) {
      delete m_ksr;
      m_ksr = nullptr;
    }
  }

  void onItemIndexSelected(int item_index) {
    if (!dock_widget->isVisible())
      return;

    Scene_points_with_normal_item *selection = qobject_cast<Scene_points_with_normal_item*>(scene->item(item_index));
    if (selection == nullptr) {
      // Keep old reference if no new point cloud has been selected.
      if (m_pwn_item == nullptr)
        dock_widget->setEnabled(false);

      return;
    }

    if (m_pwn_item == selection) {
      // The point cloud may have normals added after loading.
      if (m_pwn_item->has_normals()) {
        enable_detection(true);
        enable_regularization(true);
      }
      return;
    }

    QObject* scene_obj = dynamic_cast<QObject*>(scene);

    if (m_pwn_item) {
      disconnect(m_pwn_item, nullptr, this, SLOT(onItemChanged()));

      if (scene_obj)
        disconnect(scene_obj, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(onItemDestroyed(CGAL::Three::Scene_item*)));
    }

    m_pwn_item = selection;
    connect(m_pwn_item, SIGNAL(itemChanged()), this, SLOT(onItemChanged()));

    if (scene_obj)
      connect(scene_obj, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(onItemDestroyed(CGAL::Three::Scene_item*)));

    if (m_ksr) {
      delete m_ksr;
      m_ksr = nullptr;
    }

    dock_widget->setEnabled(true);

    if (m_pwn_item && m_pwn_item->has_normals()) {
      enable_detection(true);
      enable_regularization(true);
    }
    else {
      enable_detection(false);
      enable_regularization(false);
    }

    enable_partition(false);
    enable_reconstruction(false);

    prefill_parameters();
  }

  void onItemChanged() {
    // Enable detection if the point set item has normals now
    assert(m_pwn_item);

    if (m_pwn_item->has_normals())
      enable_detection(true);
  }

  void onItemDestroyed(CGAL::Three::Scene_item* item) {
    if (m_pwn_item == item) {
      m_pwn_item = nullptr;
      if (m_ksr) {
        delete m_ksr;
        m_ksr = nullptr;
      }

      enable_detection(false);
      enable_regularization(false);
      enable_partition(false);
      enable_reconstruction(false);
      dock_widget->hide();
    }
  }

  void run_detection() {
    assert(m_pwn_item);

    if (!m_pwn_item->has_normals())
      return;

    m_ksr->detect_planar_shapes(CGAL::parameters::maximum_distance(dock_widget->sdMaxDistanceBox->value())
      .maximum_angle(dock_widget->sdMaxAngleBox->value())
      .k_neighbors(dock_widget->sdkNeighborsBox->value())
      .minimum_region_size(dock_widget->sdMinRegionSizeBox->value())
      .regularize_parallelism(dock_widget->srParallelismCheck->isChecked())
      .regularize_coplanarity(dock_widget->srCoplanarityCheck->isChecked())
      .regularize_orthogonality(dock_widget->srOrthogonalityCheck->isChecked())
      .regularize_axis_symmetry(false)
      .angle_tolerance(dock_widget->srAngleToleranceBox->value())
      .maximum_offset(dock_widget->srMaxOffsetBox->value()));

    CGAL::Three::Three::information(QString::number(m_ksr->detected_planar_shapes().size()) + " regularized planar shapes detected");

    if (m_ksr->detected_planar_shapes().empty()) {
      enable_partition(false);
      enable_reconstruction(false);

      return;
    }

    const std::vector<CGAL::Epick::Plane_3> &planes = m_ksr->detected_planar_shapes();
    const std::vector<std::vector<std::size_t>> &regions = m_ksr->detected_planar_shape_indices();

    SMesh* mesh = new SMesh();
    std::vector<std::vector<CGAL::Epick::Point_3> > polys;
    for (std::size_t i = 0; i < regions.size(); i++)
      convex_hull(regions[i], planes[i], polys);

    for (std::size_t i = 0; i < polys.size(); i++) {
      std::vector<typename SMesh::Vertex_index> vtx(polys[i].size());
      for (std::size_t j = 0; j < polys[i].size(); j++)
        vtx[j] = mesh->add_vertex(polys[i][j]);
      mesh->add_face(vtx);
    }

    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(mesh);
    new_item->setName(tr("%1 convex shapes d%2 a%3").arg(m_pwn_item->name()).arg(dock_widget->sdMaxDistanceBox->value()).arg(dock_widget->sdMaxAngleBox->value()));
    new_item->setColor(Qt::darkCyan);
    scene->addItem(new_item);

    if (!m_known_file) {
      std::size_t max_depth = m_ksr->estimate_max_subdivision_depth();
      dock_widget->partReorientCheck->setChecked(true);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(max_depth);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(2);
    }

    enable_partition(true);
    enable_reconstruction(false);
  }

  void run_partition() {
    assert(m_pwn_item);
    assert(m_ksr);
    assert(!m_ksr->detected_planar_shapes().empty());

    m_ksr->initialize_partition(CGAL::parameters::reorient_bbox(dock_widget->partReorientCheck->checkState() == Qt::Checked)
        .max_octree_depth(dock_widget->partSubdivisionCheck->checkState() == Qt::Checked ? dock_widget->partMaxDepthBox->value() : 0)
        .max_octree_node_size(dock_widget->partPolygonsPerNodeBox->value()));

    m_ksr->partition(dock_widget->partKBox->value());

    CGAL::Three::Three::information(QString::number(m_ksr->kinetic_partition().number_of_volumes()) + " volumes in partition");

    if (m_ksr->kinetic_partition().number_of_volumes() != 0)
      enable_reconstruction(true);
    else std::cout << "kinetic partition is empty!" << std::endl;
  }

  void run_reconstruction() {
    assert(m_pwn_item);
    assert(m_ksr);
    assert(m_ksr->kinetic_partition().number_of_volumes() != 0);

    std::vector<Point_3> vtx;
    std::vector<std::vector<std::size_t> > polylist;

    std::map<typename KSR::KSP::Face_support, bool> external_nodes;

    if (dock_widget->recGroundCheck->checkState() == Qt::Checked) {
      external_nodes[KSR::KSP::Face_support::ZMIN] = false;
      m_ksr->reconstruct_with_ground(dock_widget->recLambdaBox->value(), std::back_inserter(vtx), std::back_inserter(polylist));
    }
    else {
      m_ksr->reconstruct(dock_widget->recLambdaBox->value(), external_nodes, std::back_inserter(vtx), std::back_inserter(polylist));
    }

    if (!polylist.empty())
    {
      // Add polygon mesh to scene
      Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item();
      new_item->load(vtx, polylist);
      new_item->setName(tr("%1 ksr lambda %2").arg(m_pwn_item->name()).arg(dock_widget->recLambdaBox->value()));
      new_item->setColor(Qt::darkGray);
      new_item->invalidateOpenGLBuffers();
      scene->addItem(new_item);
    }
  }

  void onSubdivisionStateChanged(int state) {
    dock_widget->partMaxDepthBox->setEnabled(state != 0);
    dock_widget->partPolygonsPerNodeBox->setEnabled(state != 0);
  }

  void onVisibilityChanged(bool) {
    if (!dock_widget->isVisible())
      return;
    std::cout << "in visibility changed" << std::endl;
  }

private:
  void convex_hull(const std::vector<std::size_t>& region, const CGAL::Epick::Plane_3& plane, std::vector<std::vector<CGAL::Epick::Point_3> > &polys) {
    if (m_pwn_item == nullptr)
      return;

    Point_set* points = m_pwn_item->point_set();

    std::vector<CGAL::Epick::Point_2> pts2d;
    pts2d.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < points->size());
      const auto& p = points->point(idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      pts2d.push_back(point);
    }
    CGAL_assertion(pts2d.size() == region.size());

    std::vector<CGAL::Epick::Point_2> ch;
    CGAL::convex_hull_2(pts2d.begin(), pts2d.end(), std::back_inserter(ch));

    std::vector<CGAL::Epick::Point_3> polygon;
    for (const auto& p : ch) {
      const auto point = plane.to_3d(p);
      polygon.push_back(point);
    }

    polys.push_back(polygon);
  }
  void enable_detection(bool enable) {
    dock_widget->shapeDetectionGroup->setEnabled(enable);
    dock_widget->sdRunButton->setEnabled(enable);
  }

  void enable_regularization(bool enable) {
    dock_widget->shapeRegularizationGroup->setEnabled(enable);
  }

  void enable_partition(bool enable) {
    dock_widget->kineticPartitionGroup->setEnabled(enable);
    dock_widget->partRunButton->setEnabled(enable);
  }

  void enable_reconstruction(bool enable) {
    dock_widget->kineticReconstructionGroup->setEnabled(enable);
    dock_widget->recRunButton->setEnabled(enable);
  }

  void prefill_parameters() {
    if (m_pwn_item == nullptr)
      return;

    std::string filename = m_pwn_item->name().toStdString();

    m_known_file = true;

    if (m_ksr) {
      delete m_ksr;
      m_ksr = nullptr;
    }

    Point_set* points = m_pwn_item->point_set();
    m_ksr = new KSR(*points);

    if (filename == "foam_box" || filename == "foam_box_new") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(0.05);
      dock_widget->sdMaxAngleBox->setValue(15);
      dock_widget->sdMinRegionSizeBox->setValue(250);
      dock_widget->sdkNeighborsBox->setValue(12);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(true);
      dock_widget->srOrthogonalityCheck->setChecked(false);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srAngleToleranceBox->setValue(10);
      dock_widget->srMaxOffsetBox->setValue(0.01);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(2);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(false);
      dock_widget->recLambdaBox->setValue(0.7);
    }
    else if (filename == "lans") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(0.15);
      dock_widget->sdMaxAngleBox->setValue(20);
      dock_widget->sdMinRegionSizeBox->setValue(300);
      dock_widget->sdkNeighborsBox->setValue(12);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(true);
      dock_widget->srOrthogonalityCheck->setChecked(false);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srAngleToleranceBox->setValue(8);
      dock_widget->srMaxOffsetBox->setValue(0.08);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(2);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(false);
      dock_widget->recLambdaBox->setValue(0.7);
    }
    else if (filename == "building_c_1M") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(1.1);
      dock_widget->sdMaxAngleBox->setValue(26);
      dock_widget->sdMinRegionSizeBox->setValue(500);
      dock_widget->sdkNeighborsBox->setValue(15);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(true);
      dock_widget->srOrthogonalityCheck->setChecked(false);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srAngleToleranceBox->setValue(3);
      dock_widget->srMaxOffsetBox->setValue(0.5);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(2);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(true);
      dock_widget->recLambdaBox->setValue(0.77);
    }
    else if (filename == "dragon") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(0.7);
      dock_widget->sdMaxAngleBox->setValue(26);
      dock_widget->sdMinRegionSizeBox->setValue(150);
      dock_widget->sdkNeighborsBox->setValue(10);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(false);
      dock_widget->srOrthogonalityCheck->setChecked(false);
      dock_widget->srCoplanarityCheck->setChecked(false);
      dock_widget->srAngleToleranceBox->setValue(0);
      dock_widget->srMaxOffsetBox->setValue(0);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(1);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(false);
      dock_widget->recLambdaBox->setValue(0.75);
    }
    else if (filename == "full_thing_pds_1M") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(0.3);
      dock_widget->sdMaxAngleBox->setValue(36);
      dock_widget->sdMinRegionSizeBox->setValue(30);
      dock_widget->sdkNeighborsBox->setValue(12);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(true);
      dock_widget->srOrthogonalityCheck->setChecked(false);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srAngleToleranceBox->setValue(3);
      dock_widget->srMaxOffsetBox->setValue(0.05);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(3);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(false);
      dock_widget->recLambdaBox->setValue(0.5);
    }
    else if (filename == "hilbert_cube2_pds_100k") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(0.3);
      dock_widget->sdMaxAngleBox->setValue(10);
      dock_widget->sdMinRegionSizeBox->setValue(10);
      dock_widget->sdkNeighborsBox->setValue(12);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(true);
      dock_widget->srOrthogonalityCheck->setChecked(true);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srAngleToleranceBox->setValue(5);
      dock_widget->srMaxOffsetBox->setValue(0.03);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(4);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(false);
      dock_widget->recLambdaBox->setValue(0.5);
    }
    else if (filename == "Meetingroom_3M") {
      // Shape detection parameters
      dock_widget->sdMaxDistanceBox->setValue(0.03);
      dock_widget->sdMaxAngleBox->setValue(19);
      dock_widget->sdMinRegionSizeBox->setValue(100);
      dock_widget->sdkNeighborsBox->setValue(15);
      // Shape regularization parameters
      dock_widget->srParallelismCheck->setChecked(true);
      dock_widget->srOrthogonalityCheck->setChecked(true);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srAngleToleranceBox->setValue(10);
      dock_widget->srMaxOffsetBox->setValue(0.03);
      // Partition parameters
      dock_widget->partReorientCheck->setChecked(false);
      dock_widget->partSubdivisionCheck->setChecked(true);
      dock_widget->partMaxDepthBox->setValue(3);
      dock_widget->partPolygonsPerNodeBox->setValue(40);
      dock_widget->partKBox->setValue(3);
      // Reconstruction parameters
      dock_widget->recGroundCheck->setChecked(false);
      dock_widget->recLambdaBox->setValue(0.5);
    }
    else {
      m_known_file = false;
      FT max_distance, max_angle;
      std::size_t min_region_size;
      m_ksr->estimate_detection_parameters(max_distance, max_angle, min_region_size);
      dock_widget->sdMaxDistanceBox->setValue(max_distance);
      dock_widget->sdMaxAngleBox->setValue(max_angle);
      dock_widget->sdMinRegionSizeBox->setValue(min_region_size);
      dock_widget->sdkNeighborsBox->setValue(12);
      dock_widget->srCoplanarityCheck->setChecked(true);
      dock_widget->srMaxOffsetBox->setValue(max_distance * 0.5);
    }
  }
};

#include "Kinetic_surface_reconstruction_plugin.moc"
