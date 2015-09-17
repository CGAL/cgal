#ifdef CGAL_GLEW_ENABLED
#include <GL/glew.h> // tmp hack to make sure gl.his included before glew.h
#endif

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Mean_curvature_flow_skeleton_plugin.h"
#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene.h"

#include "Polyhedron_type.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QMessageBox>

#include <Eigen/Sparse>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Polyhedron_copy_3.h>
#include <queue>


typedef boost::graph_traits<Polyhedron>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator            vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor        halfedge_descriptor;
typedef Polyhedron::Facet_iterator                                  Facet_iterator;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>      Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton Skeleton;

typedef Polyhedron::Traits         Kernel;
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

template<class ValueType>
struct Facet_with_id_pmap
    : public boost::put_get_helper<ValueType&,
             Facet_with_id_pmap<ValueType> >
{
    typedef Polyhedron::Face_handle key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    Facet_with_id_pmap(
      std::vector<ValueType>& internal_vector
    ) : internal_vector(internal_vector) { }

    reference operator[](key_type key) const
    { return internal_vector[key->id()]; }
private:
    std::vector<ValueType>& internal_vector;
};

class Polyhedron_demo_mean_curvature_flow_skeleton_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  QAction* actionMCFSkeleton;
  QAction* actionConvert_to_medial_skeleton;

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionMCFSkeleton" << "actionConvert_to_medial_skeleton";
  }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    mcs = NULL;
    dockWidget = NULL;
    ui = NULL;

    actionMCFSkeleton = new QAction(tr("Mean Curvature Skeleton (Advanced)"), mainWindow);
    actionMCFSkeleton->setObjectName("actionMCFSkeleton");

    actionConvert_to_medial_skeleton = new QAction(tr("Extract Medial Skeleton"), mainWindow);
    actionConvert_to_medial_skeleton->setObjectName("actionConvert_to_medial_skeleton");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);

    dockWidget = new QDockWidget(mw);
    dockWidget->setVisible(false);
    ui = new Ui::Mean_curvature_flow_skeleton_plugin();
    ui->setupUi(dockWidget);
    dockWidget->setFeatures(QDockWidget::DockWidgetMovable
                          | QDockWidget::DockWidgetFloatable
                          | QDockWidget::DockWidgetClosable);
    dockWidget->setWindowTitle("Mean Curvature Flow Skeleton");
    add_dock_widget(dockWidget);

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
    connect(dynamic_cast<Scene*>(scene), SIGNAL(updated_bbox()),
            this, SLOT(on_actionUpdateBBox()));
    connect(ui->pushButton_segment, SIGNAL(clicked()),
            this, SLOT(on_actionSegment()));

    QObject* scene_object = dynamic_cast<QObject*>(scene);
    connect(scene_object, SIGNAL(itemAboutToBeDestroyed(Scene_item*)),
            this, SLOT(on_actionItemAboutToBeDestroyed(Scene_item*)));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMCFSkeleton << actionConvert_to_medial_skeleton;
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init_ui(double diag) {
    ui->omega_H->setValue(0.1);
    ui->omega_H->setSingleStep(0.1);
    ui->omega_H->setDecimals(3);
    ui->omega_P->setValue(0.2);
    ui->omega_P->setSingleStep(0.1);
    ui->omega_P->setDecimals(3);
    ui->min_edge_length->setDecimals(7);
    ui->min_edge_length->setValue(0.002 * diag);
    ui->min_edge_length->setSingleStep(0.0000001);
    ui->delta_area->setDecimals(7);
    ui->delta_area->setValue(1e-4);
    ui->delta_area->setSingleStep(1e-5);
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
  bool is_mesh_valid(Polyhedron *pMesh) {
    if (!pMesh->is_closed())
    {
      QMessageBox msgBox;
      msgBox.setText("The mesh is not closed.");
      msgBox.exec();
      return false;
    }
    if (!pMesh->is_pure_triangle())
    {
      QMessageBox msgBox;
      msgBox.setText("The mesh is not a pure triangle mesh.");
      msgBox.exec();
      return false;
    }

    // the algorithm is only applicable on a mesh
    // that has only one connected component
    std::size_t num_component;
    CGAL::Counting_output_iterator output_it(&num_component);
    CGAL::internal::corefinement::extract_connected_components(*pMesh, output_it);
    ++output_it;
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
  bool check_mesh(Scene_polyhedron_item* item) {
    double omega_H = ui->omega_H->value();
    double omega_P = ui->omega_P->value();
    double min_edge_length = ui->min_edge_length->value();
    double delta_area = ui->delta_area->value();
    bool is_medially_centered = ui->is_medially_centered->isChecked();

    Polyhedron *pMesh = item->polyhedron();

    if (mcs == NULL)
    {
      if (!is_mesh_valid(pMesh))
      {
        return false;
      }

      mcs = new Mean_curvature_skeleton(*pMesh);
      meso_skeleton = new Polyhedron(*pMesh);
      input_triangle_mesh = pMesh;
      //set algorithm parameters
      mcs->set_quality_speed_tradeoff(omega_H);
      mcs->set_medially_centered_speed_tradeoff(omega_P);
      mcs->set_min_edge_length(min_edge_length);
      mcs->set_is_medially_centered(is_medially_centered);
      mcs->set_area_variation_factor(delta_area);

      Scene_polyhedron_item* contracted_item = new Scene_polyhedron_item( meso_skeleton );
      contracted_item->setName(QString("contracted mesh of %1").arg(item->name()));

      InputMeshItemIndex = scene->mainSelectionIndex();

      contractedItemIndex = scene->addItem(contracted_item);

      item->setVisible(false);

      fixedPointsItemIndex = -1;
      nonFixedPointsItemIndex = -1;
      poleLinesItemIndex = -1;
    }
    else
    {
      if (input_triangle_mesh != pMesh)
      {
        if (!is_mesh_valid(pMesh))
        {
          return false;
        }

        delete mcs;

        mcs = new Mean_curvature_skeleton(*pMesh);
        meso_skeleton = new Polyhedron(*pMesh);
        input_triangle_mesh = pMesh;
        //set algorithm parameters
        mcs->set_quality_speed_tradeoff(omega_H);
        mcs->set_medially_centered_speed_tradeoff(omega_P);
        mcs->set_min_edge_length(min_edge_length);
        mcs->set_is_medially_centered(is_medially_centered);
        mcs->set_area_variation_factor(delta_area);

        Scene_polyhedron_item* contracted_item = new Scene_polyhedron_item(meso_skeleton);
        contracted_item->setName(QString("contracted mesh of %1").arg(item->name()));

        InputMeshItemIndex = scene->mainSelectionIndex();

        contractedItemIndex = scene->addItem(contracted_item);

        item->setVisible(false);

        fixedPointsItemIndex = -1;
        nonFixedPointsItemIndex = -1;
        poleLinesItemIndex = -1;
      }
      else
      {
        mcs->set_quality_speed_tradeoff(omega_H);
        mcs->set_medially_centered_speed_tradeoff(omega_P);
        mcs->set_min_edge_length(min_edge_length);
        mcs->set_area_variation_factor(delta_area);
        mcs->set_is_medially_centered(is_medially_centered);
      }
    }
    return true;
  }

  void update_meso_skeleton()
  {
    CGAL::Polyhedron_copy_3<Mean_curvature_skeleton::Meso_skeleton, Polyhedron::HalfedgeDS> modifier(mcs->meso_skeleton());
    meso_skeleton->delegate(modifier);
    scene->item(contractedItemIndex)->invalidate_buffers();
    scene->itemChanged(contractedItemIndex);
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
  void on_actionUpdateBBox();
  void on_actionSegment();
  void on_actionItemAboutToBeDestroyed(Scene_item*);

private:
  Mean_curvature_skeleton* mcs;
  Polyhedron* meso_skeleton; // a copy of the meso_skeleton that is displayed
  Polyhedron* input_triangle_mesh;
  QDockWidget* dockWidget;
  Ui::Mean_curvature_flow_skeleton_plugin* ui;

  int fixedPointsItemIndex;
  int nonFixedPointsItemIndex;
  int poleLinesItemIndex;
  int contractedItemIndex;
  int InputMeshItemIndex;

  Skeleton skeleton_curve;
}; // end Polyhedron_demo_mean_curvature_flow_skeleton_plugin

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionMCFSkeleton_triggered()
{
  dockWidget->show();
  dockWidget->raise();

  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    double diag = scene->len_diagonal();
    init_ui(diag);

    fixedPointsItemIndex = -1;
    nonFixedPointsItemIndex = -1;
    poleLinesItemIndex = -1;
    contractedItemIndex = -1;
    InputMeshItemIndex = -1;
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionUpdateBBox()
{
  double diag = scene->len_diagonal();
  ui->min_edge_length->setValue(0.002 * diag);
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSegment()
{
  if (num_vertices(skeleton_curve)==0 ) on_actionSkeletonize();
  if (num_vertices(skeleton_curve)==0 ) return;

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

    // init the polyhedron simplex indices
  CGAL::set_halfedgeds_items_id(*input_triangle_mesh);

  //for each input vertex compute its distance to the skeleton
  std::vector<double> distances(num_vertices(*input_triangle_mesh));
  BOOST_FOREACH(boost::graph_traits<Skeleton>::vertex_descriptor v, vertices(skeleton_curve) )
  {
    const Point& skel_pt = skeleton_curve[v].point;
    BOOST_FOREACH(vertex_descriptor mesh_v, skeleton_curve[v].vertices)
    {
      const Point& mesh_pt = mesh_v->point();
      distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
    }
  }

  // create a property-map for sdf values
  std::vector<double> sdf_values( num_faces(*input_triangle_mesh) );
  Facet_with_id_pmap<double> sdf_property_map(sdf_values);

  // compute sdf values with skeleton
  BOOST_FOREACH(Polyhedron::Face_handle f, faces(*input_triangle_mesh))
  {
    double dist = 0;
    BOOST_FOREACH(Polyhedron::Halfedge_handle hd, halfedges_around_face(halfedge(f, *input_triangle_mesh), *input_triangle_mesh))
      dist+=distances[target(hd, *input_triangle_mesh)->id()];
    sdf_property_map[f] = dist / 3.;
  }

  // post-process the sdf values
  CGAL::sdf_values_postprocessing(*input_triangle_mesh, sdf_property_map);

  // create a property-map for segment-ids (it is an adaptor for this case)
  std::vector<std::size_t> segment_ids( num_faces(*input_triangle_mesh) );
  Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);

  // segment the mesh using default parameters
  std::cout << "Number of segments: "
            << CGAL::segmentation_from_sdf_values(*input_triangle_mesh, sdf_property_map, segment_property_map) <<"\n";

  Polyhedron* segmented_polyhedron = new Polyhedron(*input_triangle_mesh);

  int i=0;
  BOOST_FOREACH(Polyhedron::Face_handle fd, faces(*segmented_polyhedron))
  {
    fd->set_patch_id( static_cast<int>(segment_ids[i++] ));
  }

  scene->item(InputMeshItemIndex)->setVisible(false);
  Scene_polyhedron_item* item_segmentation = new Scene_polyhedron_item(segmented_polyhedron);
  scene->addItem(item_segmentation);
  item_segmentation->setName(QString("segmentation"));

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConvert_to_medial_skeleton_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if ( !is_mesh_valid(pMesh) ) return;

    QTime time;
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
                                      CGAL::IsTerminalDefault() );

    skeleton_item->setName(QString("Medial skeleton curve of %1").arg(item->name()));
    scene->addItem(skeleton_item);
    skeleton_item->invalidate_buffers();

    item->setPointsMode();

    QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionContract()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Contract...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(mcs);
  mcs->contract_geometry();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  update_meso_skeleton();
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionCollapse()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Collapse...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(mcs);
  std::size_t num_collapses = mcs->collapse_edges();
  std::cout << "collapsed " << num_collapses << " edges.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  update_meso_skeleton();
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSplit()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Split...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(mcs);
  std::size_t num_split = mcs->split_faces();
  std::cout << "split " << num_split << " triangles.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  update_meso_skeleton();
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionDegeneracy()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Detect degeneracy...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(mcs);
  mcs->detect_degeneracies();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); ++i)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }

  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }
  // update scene
  update_meso_skeleton();
  scene->item(fixedPointsItemIndex)->invalidate_buffers();
  scene->itemChanged(fixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionRun()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cout << "Run one iteration...\n";

  update_parameters(mcs);
  mcs->contract();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_interface::Item_id contracted_item_index = scene->mainSelectionIndex();
  Scene_polyhedron_item* contracted_item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(contracted_item_index));

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(contracted_item->name()));

  std::vector<Point> fixedPoints;
  mcs->fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); ++i)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }

  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(fixedPointsItemIndex, fixedPointsItem, false);
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

  Polyhedron* pMesh = item->polyhedron();
  std::vector<Point> pole_points;
  mcs->poles(pole_points);
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

  if (poleLinesItemIndex == -1)
  {
    poleLinesItemIndex = scene->addItem(poleLinesItem);
  }
  else
  {
    scene->replaceItem(poleLinesItemIndex, poleLinesItem, false);
  }
#endif

  update_meso_skeleton();
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSkeletonize()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  update_parameters(mcs);

  mcs->convert_to_skeleton(skeleton_curve);


  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  //create the polylines representing the skeleton
  Scene_polylines_item* skeleton = new Scene_polylines_item();
  skeleton->setColor(QColor(175, 0, 255));

  Polyline_visitor polyline_visitor(skeleton->polylines, skeleton_curve);
  CGAL::split_graph_into_polylines( skeleton_curve,
                                    polyline_visitor,
                                    CGAL::IsTerminalDefault() );

  skeleton->setName(QString("skeleton curve of %1").arg(item->name()));
  scene->addItem(skeleton);
  skeleton->invalidate_buffers();

  // set the fixed points and contracted mesh as invisible
  if (fixedPointsItemIndex >= 0)
  {
    scene->item(fixedPointsItemIndex)->setVisible(false);
    scene->itemChanged(fixedPointsItemIndex);
  }
  scene->item(contractedItemIndex)->setVisible(false);
  scene->itemChanged(contractedItemIndex);
  // display the original mesh in transparent mode
  item->setVisible(false);
  if (InputMeshItemIndex >= 0)
  {
    scene->item(InputMeshItemIndex)->setVisible(true);
    scene->item(InputMeshItemIndex)->setPointsMode();
    scene->itemChanged(InputMeshItemIndex);
  }

  // update scene
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConverge()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->contract_until_convergence();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); ++i)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }
  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }

  scene->item(fixedPointsItemIndex)->invalidate_buffers();
  scene->itemChanged(fixedPointsItemIndex);
  update_meso_skeleton();
  scene->setSelectedItem(index);

  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionItemAboutToBeDestroyed(Scene_item* /* item */)
{
  if (mcs != NULL)
  {
    delete mcs;
    mcs = NULL;
  }
}

#include "Polyhedron_demo_mean_curvature_flow_skeleton_plugin.moc"
