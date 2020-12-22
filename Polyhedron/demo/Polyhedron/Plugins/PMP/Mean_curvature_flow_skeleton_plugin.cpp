#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include "ui_Mean_curvature_flow_skeleton_plugin.h"

#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include "Scene_mcf_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene.h"
#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QElapsedTimer>
#include <QMessageBox>

#include <Eigen/Sparse>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Facet_with_id_pmap.h>
#include <queue>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef Scene_surface_mesh_item Scene_face_graph_item;
namespace CGAL {

template<>
void set_halfedgeds_items_id (Scene_face_graph_item::Face_graph&)
{}

} // namespace CGAL

typedef Scene_face_graph_item::Face_graph Face_graph;

typedef boost::graph_traits<Face_graph>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Face_graph>::vertex_iterator            vertex_iterator;
typedef boost::graph_traits<Face_graph>::halfedge_descriptor        halfedge_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Face_graph>      Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton Skeleton;

typedef Kernel::Point_3            Point;

struct Polyline_visitor
{
  typedef std::vector<Point> Polyline;
  typedef std::vector<std::size_t> Polyline_of_ids;

  std::list<Polyline>& polylines;
  Skeleton& skeleton;

  Polyline_visitor(std::list<Polyline>& lines, Skeleton& skeleton)
    : polylines(lines),
      skeleton(skeleton)
  {}

  void start_new_polyline()
  {
    Polyline V;
    polylines.push_back(V);
  }

  void add_node(boost::graph_traits<Skeleton>::vertex_descriptor vd)
  {
    Polyline& polyline = polylines.back();
    polyline.push_back(skeleton[vd].point);
  }

  void end_polyline(){}
};

using namespace CGAL::Three;
class Polyhedron_demo_mean_curvature_flow_skeleton_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  QAction* actionMCFSkeleton;
  QAction* actionConvert_to_medial_skeleton;

public:

  ~Polyhedron_demo_mean_curvature_flow_skeleton_plugin()
  {
    delete ui;
  }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {

    this->mw = mainWindow;
    this->scene = scene_interface;

    dockWidget = NULL;
    ui = NULL;

    actionMCFSkeleton = new QAction(tr(
                                      "Mean Curvature Skeleton (Advanced)"
                                      ), mainWindow);
    actionMCFSkeleton->setProperty("subMenuName", "Triangulated Surface Mesh Skeletonization");
    actionMCFSkeleton->setObjectName("actionMCFSkeleton");

    actionConvert_to_medial_skeleton = new QAction(tr("Extract Medial Skeleton"), mainWindow);
    actionConvert_to_medial_skeleton->setProperty("subMenuName", "Triangulated Surface Mesh Skeletonization");
    actionConvert_to_medial_skeleton->setObjectName("actionConvert_to_medial_skeleton");

    dockWidget = new QDockWidget(tr(
                                   "Mean Curvature Skeleton"
                                   ),mw);
    dockWidget->setVisible(false);
    ui = new Ui::Mean_curvature_flow_skeleton_plugin();
    ui->setupUi(dockWidget);
    dockWidget->setFeatures(QDockWidget::DockWidgetMovable
                          | QDockWidget::DockWidgetFloatable
                          | QDockWidget::DockWidgetClosable);
    dockWidget->setWindowTitle(tr(
                               "Mean Curvature Skeleton"
                                 ));
    addDockWidget(dockWidget);

    connect(ui->pushButton_contract, SIGNAL(clicked()),
            this, SLOT(on_actionContract()));
    connect(ui->pushButton_collapse, SIGNAL(clicked()),
            this, SLOT(on_actionCollapse()));
    connect(ui->pushButton_split, SIGNAL(clicked()),
            this, SLOT(on_actionSplit()));
    connect(ui->pushButton_degeneracy, SIGNAL(clicked()),
            this, SLOT(on_actionDegeneracy()));
    connect(ui->pushButton_run, SIGNAL(clicked()),
            this, SLOT(on_actionRun()));
    connect(ui->pushButton_skeletonize, SIGNAL(clicked()),
            this, SLOT(on_actionSkeletonize()));
    connect(ui->pushButton_converge, SIGNAL(clicked()),
            this, SLOT(on_actionConverge()));
    connect(dynamic_cast<Scene*>(scene), SIGNAL(updated_bbox(bool)),
            this, SLOT(on_actionUpdateBBox(bool)));
    connect(ui->pushButton_segment, SIGNAL(clicked()),
            this, SLOT(on_actionSegment()));

    autoConnectActions();
    QObject* scene_object = dynamic_cast<QObject*>(scene);
    connect(scene_object, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
            this, SLOT(on_actionItemAboutToBeDestroyed(CGAL::Three::Scene_item*)));
  }

  virtual void closure()
  {
    dockWidget->hide();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMCFSkeleton << actionConvert_to_medial_skeleton;
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init_ui(double diag) {
    on_checkbox_toggled(false);
    connect(ui->is_medially_centered, SIGNAL(toggled(bool)),
            this, SLOT(on_checkbox_toggled(bool)));
    connect(ui->helpButton, &QPushButton::clicked,
            [this]{QMessageBox::about(mw, QString("Help"),
                                    QString("This widget gives access to the low level steps of the mean curvature flow sketonization algorithm. "
                                            "The algorithm is iterative. Each iteration consist in calls to Contract, Collapse, Split, "
                                            "and Degeneracy (repectively mesh contraction, edge collapse, edge split, and degenerate edge"
                                            "removal). The skeleton extraction can be called at any time but for a better result it should be"
                                            "called when the iterations are converging. A segmentation of the surface can be extracted using"
                                            "the distance of the mesh to the skeleton computed.\n"
                                             "All operations can be applied to a polyhedron item or "
                                            "to a surface mesh item. The generated mcf group must be selected in "
                                            "order to continue an on-going set of operations. "));});
    ui->omega_H->setValue(0.1);
    ui->omega_P->setValue(0.2);
    ui->min_edge_length->setValue(0.002 * diag);
    ui->delta_area->setValue(1e-4);
    ui->is_medially_centered->setChecked(false);

    ui->label_omega_H->setToolTip(QString("omega_H controls the velocity of movement and approximation quality"));
    ui->label_omega_P->setToolTip(QString("omega_P controls the smoothness of the medial approximation"));
    ui->pushButton_contract->setToolTip(QString("contract mesh based on mean curvature flow"));
    ui->pushButton_collapse->setToolTip(QString("collapse short edges"));
    ui->pushButton_split->setToolTip(QString("split obtuse triangles"));
    ui->pushButton_degeneracy->setToolTip(QString("fix degenerate points"));
    ui->pushButton_skeletonize->setToolTip(QString("Turn mesh to a skeleton curve"));
    ui->pushButton_run->setToolTip(QString("run one iteration of contract, collapse, split, detect degeneracy"));
    ui->pushButton_converge->setToolTip(QString("iteratively contract the mesh until convergence"));
  }

  bool check_item_index(int index) {
    if (index < 0)
    {
      QMessageBox msgBox;
      msgBox.setText("Please select an item first");
      msgBox.exec();
      return false;
    }
    return true;
  }

  /// \todo move this function into an include
  bool is_mesh_valid(Face_graph *pMesh) {
    if (! CGAL::is_closed(*pMesh))
    {
      QMessageBox msgBox;
      msgBox.setText("The mesh is not closed.");
      msgBox.exec();
      return false;
    }
    if (! CGAL::is_triangle_mesh(*pMesh))
    {
      QMessageBox msgBox;
      msgBox.setText("The mesh is not a pure triangle mesh.");
      msgBox.exec();
      return false;
    }

    // the algorithm is only applicable on a mesh
    // that has only one connected component

    boost::unordered_map<boost::graph_traits<Face_graph>::face_descriptor,int> cc(num_faces(*pMesh));
    std::size_t num_component = PMP::connected_components(*pMesh, boost::make_assoc_property_map(cc));

    if (num_component != 1)
    {
      QMessageBox msgBox;
      QString str = QString("The mesh is not a single closed mesh.\n It has %1 components.").arg(num_component);
      msgBox.setText(str);
      msgBox.exec();
      return false;
    }
    return true;
  }

  /// \todo remove duplicated code
  // check if the Mean_curvature_skeleton exists
  // or has the same polyheron item
  // check if the mesh is a watertight triangle mesh
  bool check_mesh(Scene_mcf_item* item) {
    Face_graph *pMesh = item->input_triangle_mesh;

    if (item->mcs == NULL)
    {
      if (!is_mesh_valid(pMesh))
      {
        return false;
      }
      createContractedItem(item);

      item->fixedPointsItemIndex = -1;
      item->nonFixedPointsItemIndex = -1;
      item->poleLinesItemIndex = -1;
    }
    else
    {
      item->mcs->set_quality_speed_tradeoff(ui->omega_H->value());
      item->mcs->set_medially_centered_speed_tradeoff(ui->omega_P->value());
      item->mcs->set_min_edge_length(ui->min_edge_length->value());
      item->mcs->set_area_variation_factor(ui->delta_area->value());
      item->mcs->set_is_medially_centered(ui->is_medially_centered->isChecked());
    }
    return true;
  }

  void update_meso_skeleton(Scene_mcf_item* item)
  {
    clear(*item->meso_skeleton);
    copy_face_graph(item->mcs->meso_skeleton(), *item->meso_skeleton);
    scene->item(item->contractedItemIndex)->invalidateOpenGLBuffers();
    scene->itemChanged(item->contractedItemIndex);
  }

  void update_parameters(Mean_curvature_skeleton* mcs)
  {
    double omega_H = ui->omega_H->value();
    double omega_P = ui->omega_P->value();
    double min_edge_length = ui->min_edge_length->value();
    double delta_area = ui->delta_area->value();
    bool is_medially_centered = ui->is_medially_centered->isChecked();

    mcs->set_quality_speed_tradeoff(omega_H);
    mcs->set_medially_centered_speed_tradeoff(omega_P);
    mcs->set_min_edge_length(min_edge_length);
    mcs->set_area_variation_factor(delta_area);
    mcs->set_is_medially_centered(is_medially_centered);
  }

public Q_SLOTS:
  void on_actionMCFSkeleton_triggered();
  void on_actionConvert_to_medial_skeleton_triggered();
  void on_actionContract();
  void on_actionCollapse();
  void on_actionSplit();
  void on_actionDegeneracy();
  void on_actionRun();
  void on_actionSkeletonize();
  void on_actionConverge();
  void on_actionUpdateBBox(bool);
  void on_checkbox_toggled(bool);
  void on_actionSegment();
  void on_actionItemAboutToBeDestroyed(CGAL::Three::Scene_item*);

private:
  Scene_mcf_item *getMCFItem();
  void createContractedItem(Scene_mcf_item* item);
  QDockWidget* dockWidget;
  Ui::Mean_curvature_flow_skeleton_plugin* ui;

}; // end Polyhedron_demo_mean_curvature_flow_skeleton_plugin

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionMCFSkeleton_triggered()
{
  dockWidget->show();
  dockWidget->raise();
  double diag = scene->len_diagonal();
  init_ui(diag);
  getMCFItem();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_checkbox_toggled(bool b)
{
  ui->omega_P->setEnabled(b);
  ui->label_omega_P->setEnabled(b);
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionUpdateBBox(bool)
{
  double diag = scene->len_diagonal();
  ui->min_edge_length->setValue(0.002 * diag);
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSegment()
{
  Scene_mcf_item* item = getMCFItem();

  if(!item)
  {
    return;
  }

  if (num_vertices(item->skeleton_curve)==0 ) on_actionSkeletonize();
  if (num_vertices(item->skeleton_curve)==0 ) { QApplication::restoreOverrideCursor(); return;}
  QApplication::setOverrideCursor(Qt::WaitCursor);

  QElapsedTimer time;
  time.start();

    // init the polyhedron simplex indices
  CGAL::set_halfedgeds_items_id(*item->input_triangle_mesh);
  boost::property_map<Face_graph, boost::vertex_index_t>::type
    vimap = get(boost::vertex_index, *item->input_triangle_mesh);

  //for each input vertex compute its distance to the skeleton
  std::vector<double> distances(num_vertices(*item->input_triangle_mesh));

  Face_graph *smesh = item->input_triangle_mesh;

  boost::property_map<Face_graph,CGAL::vertex_point_t>::type vpm
    = get(CGAL::vertex_point,*smesh);

  for(boost::graph_traits<Skeleton>::vertex_descriptor v : CGAL::make_range(vertices(item->skeleton_curve)) )
  {
    const Point& skel_pt = item->skeleton_curve[v].point;
    for(vertex_descriptor mesh_v : item->skeleton_curve[v].vertices)
    {
      const Point& mesh_pt = get(vpm,mesh_v);
      distances[get(vimap,mesh_v)] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
    }
  }

  // create a property-map for sdf values
  std::vector<double> sdf_values( num_faces(*item->input_triangle_mesh) );
  Facet_with_id_pmap<Face_graph,double> sdf_property_map(*item->input_triangle_mesh, sdf_values);

  // compute sdf values with skeleton
  for(boost::graph_traits<Face_graph>::face_descriptor f : faces(*item->input_triangle_mesh))
  {
    double dist = 0;
    for(boost::graph_traits<Face_graph>::halfedge_descriptor hd : halfedges_around_face(halfedge(f, *item->input_triangle_mesh), *item->input_triangle_mesh))
      dist+=distances[get(vimap,target(hd, *item->input_triangle_mesh))];
    sdf_property_map[f] = dist / 3.;
  }

  // post-process the sdf values
  CGAL::sdf_values_postprocessing(*item->input_triangle_mesh, sdf_property_map);

  // create a property-map for segment-ids (it is an adaptor for this case)
  std::vector<std::size_t> segment_ids( num_faces(*item->input_triangle_mesh) );
  Facet_with_id_pmap<Face_graph,std::size_t> segment_property_map(*item->input_triangle_mesh, segment_ids);

  // segment the mesh using default parameters
  std::cout << "Number of segments: "
            << CGAL::segmentation_from_sdf_values(*item->input_triangle_mesh, sdf_property_map, segment_property_map) <<"\n";

  Face_graph* segmented_polyhedron = new Face_graph(*item->input_triangle_mesh);

  Scene_face_graph_item* item_segmentation = new Scene_face_graph_item(segmented_polyhedron);
  int i=0;
  typedef boost::property_map<Face_graph, CGAL::face_patch_id_t<int> >::type Fpim;
  Fpim fpim = get(CGAL::face_patch_id_t<int>(), *segmented_polyhedron);
  int nb_segment=0;
  for(boost::graph_traits<Face_graph>::face_descriptor fd : faces(*segmented_polyhedron))
  {
    int segment = static_cast<int>(segment_ids[i++]);
    if(segment > nb_segment)
      nb_segment = segment + 1;
    put(fpim, fd, segment);

  }
  item_segmentation->setItemIsMulticolor(true);
  item_segmentation->computeItemColorVectorAutomatically(true);
  item_segmentation->setProperty("NbPatchIds", nb_segment); //for join_and_split plugin
  item_segmentation->invalidateOpenGLBuffers();
  scene->addItem(item_segmentation);
  item_segmentation->setName(QString("segmentation"));
  scene->changeGroup(item_segmentation, item);
  Scene_item* parent = scene->item(item->InputMeshItemIndex);
  if(parent)
    parent->setVisible(false);
  scene->setSelectedItem(scene->item_id(item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConvert_to_medial_skeleton_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_face_graph_item* item =
    qobject_cast<Scene_face_graph_item*>(scene->item(index));

  if(item)
  {
    Face_graph* pMesh = item->polyhedron();

    if ( !is_mesh_valid(pMesh) ) return;

    QElapsedTimer time;
    time.start();
    QApplication::setOverrideCursor(Qt::WaitCursor);

    Skeleton skeleton;
    CGAL::extract_mean_curvature_flow_skeleton(*pMesh, skeleton);

    std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

    //create the polylines representing the skeleton
    Scene_polylines_item* skeleton_item = new Scene_polylines_item();
    skeleton_item->setColor(QColor(175, 0, 255));

    Polyline_visitor polyline_visitor(skeleton_item->polylines, skeleton);
    CGAL::split_graph_into_polylines( skeleton,
                                      polyline_visitor,
                                      CGAL::internal::IsTerminalDefault() );

    skeleton_item->setName(QString("Medial skeleton curve of %1").arg(item->name()));
    scene->setSelectedItem(-1);
    scene->addItem(skeleton_item);
    skeleton_item->invalidateOpenGLBuffers();

    item->setPointsMode();

    QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionContract()
{
  Scene_mcf_item* item = getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  std::cout << "Contract...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(item->mcs);
  item->mcs->contract_geometry();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  update_meso_skeleton(item);
  QApplication::restoreOverrideCursor();
  scene->setSelectedItem(scene->item_id(item));
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionCollapse()
{
  Scene_mcf_item* item =
    getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  std::cout << "Collapse...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(item->mcs);
  std::size_t num_collapses = item->mcs->collapse_edges();
  std::cout << "collapsed " << num_collapses << " edges.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  update_meso_skeleton(item);
  QApplication::restoreOverrideCursor();
  scene->setSelectedItem(scene->item_id(item));
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSplit()
{
  Scene_mcf_item* item =
    getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  std::cout << "Split...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(item->mcs);
  std::size_t num_split = item->mcs->split_faces();
  std::cout << "split " << num_split << " triangles.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  update_meso_skeleton(item);
  QApplication::restoreOverrideCursor();
  scene->setSelectedItem(scene->item_id(item));
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionDegeneracy()
{
  Scene_mcf_item* item = getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  std::cout << "Detect degeneracy...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(item->mcs);
  item->mcs->detect_degeneracies();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));
  std::vector<Point> fixedPoints;
  item->mcs->fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); ++i)
  {
    Kernel::Point_3 point (fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->insert(point);
  }
  ps->select_all ();

  if (item->fixedPointsItemIndex == -1)
  {
    item->fixedPointsItemIndex = scene->addItem(fixedPointsItem);
    scene->changeGroup(fixedPointsItem, item);
    item->lockChild(fixedPointsItem);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(item->fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }
  // update scene
  update_meso_skeleton(item);
  scene->item(item->fixedPointsItemIndex)->invalidateOpenGLBuffers();
  scene->itemChanged(item->fixedPointsItemIndex);
  scene->setSelectedItem(scene->item_id(item));
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionRun()
{
  Scene_mcf_item* item = getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cout << "Run one iteration...\n";
  Scene_face_graph_item* contracted_item = NULL;
if(item->contractedItemIndex != -1)
    contracted_item = qobject_cast<Scene_face_graph_item*>(scene->item(item->contractedItemIndex));
scene->setSelectedItem(scene->item_id(item));
//todo : create a new contracted item
if(!contracted_item)
{
  createContractedItem(item);
  contracted_item = qobject_cast<Scene_face_graph_item*>(scene->item(item->contractedItemIndex));
}

  update_parameters(item->mcs);
  item->mcs->contract();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(contracted_item->name()));
  std::vector<Point> fixedPoints;
  item->mcs->fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); ++i)
  {
    Kernel::Point_3 point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->insert(point);
  }
  ps->select_all();

  if (item->fixedPointsItemIndex == -1)
  {
    item->fixedPointsItemIndex = scene->addItem(fixedPointsItem);
    scene->changeGroup(fixedPointsItem, item);
    item->lockChild(fixedPointsItem);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(item->fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }

//#define DRAW_NON_FIXED_POINTS
#ifdef DRAW_NON_FIXED_POINTS
  // draw non-fixed points
  Scene_points_with_normal_item* nonFixedPointsItem = new Scene_points_with_normal_item;
  nonFixedPointsItem->setName("non-fixed points");
  nonFixedPointsItem->setColor(QColor(0, 255, 0));
  std::vector<Point> nonFixedPoints;
  mcs->non_fixed_points(nonFixedPoints);
  ps = nonFixedPointsItem->point_set();
  for (size_t i = 0; i < nonFixedPoints.size(); ++i)
  {
    UI_point_3<Kernel> point(nonFixedPoints[i].x(), nonFixedPoints[i].y(), nonFixedPoints[i].z());
    ps->push_back(point);
  }
  if (nonFixedPointsItemIndex == -1)
  {
    nonFixedPointsItemIndex = scene->addItem(nonFixedPointsItem);
    scene->changeGroup(nonFixedPointsItem, item);
    item->lockChild(nonFixedPointsItem);
  }
  else
  {
    scene->replaceItem(nonFixedPointsItemIndex, nonFixedPointsItem, false);
  }
  scene->itemChanged(nonFixedPointsItemIndex);
#endif

//#define DRAW_POLE_LINE
#ifdef DRAW_POLE_LINE
  // draw lines connecting surface points and their correspondent poles
  Scene_polylines_item* poleLinesItem = new Scene_polylines_item();
  Face_graph* pMesh = item->input_triangle_mesh;
  std::vector<Point> pole_points;
  item->mcs->poles(pole_points);
  vertex_iterator vb, ve;
  int id = 0;
  for (boost::tie(vb, ve) = vertices(*pMesh); vb != ve; ++vb)
  {
    std::vector<Point> line;
    line.clear();

    vertex_descriptor v = *vb;
    Point s = v->point();
    Point t = pole_points[id++];

    line.push_back(s);
    line.push_back(t);
    poleLinesItem->polylines.push_back(line);
  }

  if (item->poleLinesItemIndex == -1)
  {
    item->poleLinesItemIndex = scene->addItem(poleLinesItem);
    scene->changeGroup(poleLinesItem, item);
    item->lockChild(poleLinesItem);
  }
  else
  {
    scene->replaceItem(poleLinesItemIndex, poleLinesItem, false);
  }
#endif

  update_meso_skeleton(item);
  scene->setSelectedItem(scene->item_id(item));
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSkeletonize()
{
  Scene_mcf_item* item = getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(item->mcs);

  item->mcs->convert_to_skeleton(item->skeleton_curve);


  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  //create the polylines representing the skeleton
  Scene_polylines_item* skeleton = new Scene_polylines_item();
  skeleton->setColor(QColor(175, 0, 255));

  Polyline_visitor polyline_visitor(skeleton->polylines, item->skeleton_curve);
  CGAL::split_graph_into_polylines( item->skeleton_curve,
                                    polyline_visitor,
                                    CGAL::internal::IsTerminalDefault() );

  skeleton->setName(QString("skeleton curve of %1").arg(item->name()));
  skeleton->invalidateOpenGLBuffers();
  if(item->skeletonItemIndex == -1)
  {
    item->skeletonItemIndex = scene->addItem(skeleton);
    scene->changeGroup(skeleton, item);
    item->lockChild(skeleton);
  }
  else
  {
    scene->replaceItem(item->skeletonItemIndex, skeleton, false);
  }

  // set the fixed points and contracted mesh as invisible
  if (item->fixedPointsItemIndex >= 0)
  {
    scene->item(item->fixedPointsItemIndex)->setVisible(false);
    scene->itemChanged(item->fixedPointsItemIndex);
  }
  scene->item(item->contractedItemIndex)->setVisible(false);
  scene->itemChanged(item->contractedItemIndex);
  if (item->InputMeshItemIndex >= 0)
  {
    scene->item(item->InputMeshItemIndex)->setVisible(true);
    scene->item(item->InputMeshItemIndex)->setPointsMode();
    scene->itemChanged(item->InputMeshItemIndex);
  }
  scene->setSelectedItem(scene->item_id(item));
  // update scene
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConverge()
{
  Scene_mcf_item* item = getMCFItem();

  if (!item || !check_mesh(item))
  {
    return;
  }

  QElapsedTimer time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  item->mcs->contract_until_convergence();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  item->mcs->fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); ++i)
  {
    Kernel::Point_3 point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->insert(point);
  }
  ps->select_all();

  if (item->fixedPointsItemIndex == -1)
  {
    item->fixedPointsItemIndex = scene->addItem(fixedPointsItem);
    scene->changeGroup(fixedPointsItem, item);
    item->lockChild(fixedPointsItem);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(item->fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }

  scene->item(item->fixedPointsItemIndex)->invalidateOpenGLBuffers();
  scene->itemChanged(item->fixedPointsItemIndex);
  update_meso_skeleton(item);
  scene->setSelectedItem(scene->item_id(item));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionItemAboutToBeDestroyed(CGAL::Three::Scene_item* corpse )
{
  Scene_mcf_item *mcf= qobject_cast<Scene_mcf_item*>(corpse);

  if(mcf)
  {
    mcf->mcs = NULL;
    mcf->meso_skeleton = NULL;
    mcf->input_triangle_mesh = NULL;
    mcf->fixedPointsItemIndex = -1;
    mcf->nonFixedPointsItemIndex = -1;
    mcf->poleLinesItemIndex = -1;
    mcf->contractedItemIndex = -1;
    mcf->InputMeshItemIndex = -1;
    mcf->meso_skeleton = NULL;
    mcf->input_triangle_mesh = NULL;
  }
}

void
Polyhedron_demo_mean_curvature_flow_skeleton_plugin::createContractedItem(Scene_mcf_item* item)
{
  if(!item)
    return;
  if(item->mcs != NULL)
    delete item->mcs;
  double omega_H = ui->omega_H->value();
  double omega_P = ui->omega_P->value();
  double min_edge_length = ui->min_edge_length->value();
  double delta_area = ui->delta_area->value();
  bool is_medially_centered = ui->is_medially_centered->isChecked();

  item->mcs = new Mean_curvature_skeleton(*item->input_triangle_mesh);
  item->meso_skeleton = new Face_graph(*item->input_triangle_mesh);
  //set algorithm parameters
  item->mcs->set_quality_speed_tradeoff(omega_H);
  item->mcs->set_medially_centered_speed_tradeoff(omega_P);
  item->mcs->set_min_edge_length(min_edge_length);
  item->mcs->set_is_medially_centered(is_medially_centered);
  item->mcs->set_area_variation_factor(delta_area);

  Scene_face_graph_item* contracted_item = new Scene_face_graph_item(item->meso_skeleton);
  contracted_item->setName(QString("contracted mesh of %1").arg(item->name()));
  contracted_item->setItemIsMulticolor(false); //avoids segfault if item was a multicolor surface_mesh

  item->contractedItemIndex = scene->addItem(contracted_item);
  scene->changeGroup(contracted_item, item);
  item->lockChild(contracted_item);
  scene->setSelectedItem(scene->item_id(item));

}

Scene_mcf_item*
Polyhedron_demo_mean_curvature_flow_skeleton_plugin::getMCFItem()
{
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_mcf_item* mcf = qobject_cast<Scene_mcf_item*>(scene->item(index));
    if(mcf)
      return mcf;
  }

  //if the selected item is not an MCF but is a face_graph_item,
  //then create and add a new MCF
  if(scene->mainSelectionIndex() != -1)
  {
    Scene_face_graph_item* item =
        qobject_cast<Scene_face_graph_item*>(scene->item(
                                               scene->mainSelectionIndex()));
    if(item)
    {
      Face_graph* pMesh = item->face_graph();

      if(!pMesh) return NULL;
      Scene_mcf_item* mcf = new Scene_mcf_item(item->face_graph(),
                                               scene->mainSelectionIndex(),
                                               QString("%1 (mcf)").arg(item->name()));
      connect(item, &Scene_face_graph_item::aboutToBeDestroyed,
              [mcf, this]{
        if(scene->item_id(mcf) != -1){
          scene->erase(scene->item_id(mcf));
      }});
      scene->setSelectedItem(scene->addItem(mcf));
      item->setVisible(false);
      scene->itemChanged(item);
      return mcf;
    }
  }
  return NULL;
}

#include "Mean_curvature_flow_skeleton_plugin.moc"
