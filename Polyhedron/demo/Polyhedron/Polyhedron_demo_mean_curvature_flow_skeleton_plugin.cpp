#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Mean_curvature_flow_skeleton_plugin.h"

#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_points_with_normal_item.cpp"
#include "Scene_polylines_item.h"
#include "Scene_polylines_item.cpp"

#include "Polyhedron_type.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>

#include <Eigen/Sparse>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Mean_curvature_skeleton.h>

template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;

    reference operator[](key_type key) const { return key->id(); }
};

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<CGAL::Eigen_sparse_matrix<double>::EigenType> > Sparse_linear_solver;

typedef CGAL::Mean_curvature_skeleton<Polyhedron, Sparse_linear_solver, Vertex_index_map, Edge_index_map> Mean_curvature_skeleton;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

typedef Polyhedron::Traits         Kernel;
typedef Kernel::Point_3            Point;

class Polyhedron_demo_mean_curvature_flow_skeleton_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionMCFSkeleton";
  }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    mcs = NULL;
    dockWidget = NULL;
    ui = NULL;

    std::cerr << "init plugin\n";
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable() const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init_ui(double diag) {
    ui->omega_L->setValue(1);
    ui->omega_H->setValue(0.1);
    ui->edgelength_TH->setDecimals(7);
    ui->edgelength_TH->setValue(0.002 * diag);
    ui->alpha->setValue(0.15);
    ui->zero_TH->setDecimals(8);
    ui->zero_TH->setValue(1e-07);
    ui->area_TH->setDecimals(8);
    ui->area_TH->setValue(1e-5);
  }

public slots:
  void on_actionMCFSkeleton_triggered();
  void on_actionContract();
  void on_actionCollapse();
  void on_actionSplit();
  void on_actionDegeneracy();
  void on_actionRun();
  void on_actionSkeletonize();
  void on_actionConverge();

private:
  Mean_curvature_skeleton* mcs;
  QDockWidget* dockWidget;
  Ui::Mean_curvature_flow_skeleton_plugin* ui;
  int fixedPointsItemIndex;
  int nonFixedPointsItemIndex;
}; // end Polyhedron_demo_mean_curvature_flow_skeleton_plugin

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionMCFSkeleton_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    dockWidget = new QDockWidget(mw);
    ui = new Ui::Mean_curvature_flow_skeleton_plugin();
    ui->setupUi(dockWidget);
    dockWidget->setFeatures(QDockWidget::DockWidgetMovable
                          | QDockWidget::DockWidgetFloatable
                          | QDockWidget::DockWidgetClosable);
    dockWidget->setWindowTitle("Mean Curvature Flow Skeleton");
    mw->addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    mw->tabifyDockWidget(static_cast<MainWindow*>(mw)->get_ui()->consoleDockWidget, dockWidget);
    dockWidget->show();
    dockWidget->raise();

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

    double diag = scene->len_diagonal();
    init_ui(diag);

    fixedPointsItemIndex = -1;
    nonFixedPointsItemIndex = -1;
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionContract()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  double omega_L = ui->omega_L->value();
  double omega_H = ui->omega_H->value();
  double edgelength_TH = ui->edgelength_TH->value();
  double alpha = ui->alpha->value();
  double zero_TH = ui->zero_TH->value();
  double area_TH = ui->area_TH->value();
  double diag = scene->len_diagonal();

  if (mcs == NULL)
  {
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }
  else
  {
    Polyhedron* mesh = mcs->get_polyhedron();
    if (mesh != pMesh)
    {
      delete mcs;
      init_ui(diag);
      omega_L = ui->omega_L->value();
      omega_H = ui->omega_H->value();
      edgelength_TH = ui->edgelength_TH->value();
      alpha = ui->alpha->value();
      zero_TH = ui->zero_TH->value();
      area_TH = ui->area_TH->value();
      mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                        omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
    }
    else
    {
      mcs->set_omega_L(omega_L);
      mcs->set_omega_H(omega_H);
      mcs->set_edgelength_TH(edgelength_TH);
      mcs->set_zero_TH(zero_TH);
    }
  }

  QTime time;
  time.start();
  std::cout << "Contract...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->contract_geometry();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  scene->itemChanged(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionCollapse()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (mcs == NULL)
  {
    double omega_L = ui->omega_L->value();
    double omega_H = ui->omega_H->value();
    double edgelength_TH = ui->edgelength_TH->value();
    double alpha = ui->alpha->value();
    double zero_TH = ui->zero_TH->value();
    double area_TH = ui->area_TH->value();
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }

  QTime time;
  time.start();
  std::cout << "Collapse...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cout << "before collapse edges\n";
  int num_collapses = mcs->collapse_short_edges();
  std::cout << "collapse " << num_collapses << " edges.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  item->color();
  // update scene
  scene->itemChanged(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSplit()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (mcs == NULL)
  {
    double omega_L = ui->omega_L->value();
    double omega_H = ui->omega_H->value();
    double edgelength_TH = ui->edgelength_TH->value();
    double alpha = ui->alpha->value();
    double zero_TH = ui->zero_TH->value();
    double area_TH = ui->area_TH->value();
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }

  QTime time;
  time.start();
  std::cout << "Split...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cout << "before split triangles\n";
  int num_split = mcs->iteratively_split_triangles();
  std::cout << "split " << num_split << " triangles.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  scene->itemChanged(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionDegeneracy()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (mcs == NULL)
  {
    double omega_L = ui->omega_L->value();
    double omega_H = ui->omega_H->value();
    double edgelength_TH = ui->edgelength_TH->value();
    double alpha = ui->alpha->value();
    double zero_TH = ui->zero_TH->value();
    double area_TH = ui->area_TH->value();
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }

  QTime time;
  time.start();
  std::cout << "Degeneracy\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->detect_degeneracies();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->get_fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); i++)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }

  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem);
    std::cerr << "add item " << fixedPointsItemIndex << "\n";
  }
  else
  {
    std::cerr << "replace item " << fixedPointsItemIndex << "\n";
    scene->replaceItem(fixedPointsItemIndex, fixedPointsItem);
  }
  // update scene
  scene->itemChanged(index);
  scene->itemChanged(fixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionRun()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  double omega_L = ui->omega_L->value();
  double omega_H = ui->omega_H->value();
  double edgelength_TH = ui->edgelength_TH->value();
  double alpha = ui->alpha->value();
  double zero_TH = ui->zero_TH->value();
  double area_TH = ui->area_TH->value();
  double diag = scene->len_diagonal();

  if (mcs == NULL)
  {
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }
  else
  {
    Polyhedron* mesh = mcs->get_polyhedron();
    if (mesh != pMesh)
    {
      delete mcs;
      init_ui(diag);
      omega_L = ui->omega_L->value();
      omega_H = ui->omega_H->value();
      edgelength_TH = ui->edgelength_TH->value();
      alpha = ui->alpha->value();
      zero_TH = ui->zero_TH->value();
      area_TH = ui->area_TH->value();
      mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                        omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
    }
    else
    {
      mcs->set_omega_L(omega_L);
      mcs->set_omega_H(omega_H);
      mcs->set_edgelength_TH(edgelength_TH);
      mcs->set_zero_TH(zero_TH);
    }
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->contract();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->get_fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); i++)
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
    scene->replaceItem(fixedPointsItemIndex, fixedPointsItem);
  }

//  Scene_points_with_normal_item* nonFixedPointsItem = new Scene_points_with_normal_item;
//  nonFixedPointsItem->setName("non-fixed points");
//  nonFixedPointsItem->setColor(QColor(0, 255, 0));
//  std::vector<Point> nonFixedPoints;
//  mcs->get_non_fixed_points(nonFixedPoints);
//  ps = nonFixedPointsItem->point_set();
//  for (size_t i = 0; i < nonFixedPoints.size(); i++)
//  {
//    UI_point_3<Kernel> point(nonFixedPoints[i].x(), nonFixedPoints[i].y(), nonFixedPoints[i].z());
//    ps->push_back(point);
//  }
//  if (nonFixedPointsItemIndex == -1)
//  {
//    nonFixedPointsItemIndex = scene->addItem(nonFixedPointsItem);
//  }
//  else
//  {
//    scene->replaceItem(nonFixedPointsItemIndex, nonFixedPointsItem);
//  }

  scene->itemChanged(index);
  scene->itemChanged(fixedPointsItemIndex);
//  scene->itemChanged(nonFixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSkeletonize()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (mcs == NULL)
  {
    double omega_L = ui->omega_L->value();
    double omega_H = ui->omega_H->value();
    double edgelength_TH = ui->edgelength_TH->value();
    double alpha = ui->alpha->value();
    double zero_TH = ui->zero_TH->value();
    double area_TH = ui->area_TH->value();
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Graph g;
  std::vector<Point> points;

  mcs->convert_to_skeleton();
  mcs->get_skeleton(g, points);

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_polylines_item* skeleton = new Scene_polylines_item();

  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    std::vector<Point> line;
    line.clear();
    Point s = points[boost::source(*ei, g)];
    Point t = points[boost::target(*ei, g)];
    line.push_back(s);
    line.push_back(t);
    skeleton->polylines.push_back(line);
  }
  skeleton->setName(QString("skeleton curve of %1").arg(item->name()));
  scene->addItem(skeleton);

  // update scene
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConverge()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  double omega_L = ui->omega_L->value();
  double omega_H = ui->omega_H->value();
  double edgelength_TH = ui->edgelength_TH->value();
  double alpha = ui->alpha->value();
  double zero_TH = ui->zero_TH->value();
  double area_TH = ui->area_TH->value();
  double diag = scene->len_diagonal();

  if (mcs == NULL)
  {
    mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
  }
  else
  {
    Polyhedron* mesh = mcs->get_polyhedron();
    if (mesh != pMesh)
    {
      delete mcs;
      init_ui(diag);
      omega_L = ui->omega_L->value();
      omega_H = ui->omega_H->value();
      edgelength_TH = ui->edgelength_TH->value();
      alpha = ui->alpha->value();
      zero_TH = ui->zero_TH->value();
      area_TH = ui->area_TH->value();
      mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                        omega_L, omega_H, edgelength_TH, zero_TH, area_TH);
    }
    else
    {
      mcs->set_omega_L(omega_L);
      mcs->set_omega_H(omega_H);
      mcs->set_edgelength_TH(edgelength_TH);
      mcs->set_zero_TH(zero_TH);
    }
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->run_to_converge();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->get_fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); i++)
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
    scene->replaceItem(fixedPointsItemIndex, fixedPointsItem);
  }

//  Scene_points_with_normal_item* nonFixedPointsItem = new Scene_points_with_normal_item;
//  nonFixedPointsItem->setName("non-fixed points");
//  nonFixedPointsItem->setColor(QColor(0, 255, 0));
//  std::vector<Point> nonFixedPoints;
//  mcs->get_non_fixed_points(nonFixedPoints);
//  ps = nonFixedPointsItem->point_set();
//  for (size_t i = 0; i < nonFixedPoints.size(); i++)
//  {
//    UI_point_3<Kernel> point(nonFixedPoints[i].x(), nonFixedPoints[i].y(), nonFixedPoints[i].z());
//    ps->push_back(point);
//  }
//  if (nonFixedPointsItemIndex == -1)
//  {
//    nonFixedPointsItemIndex = scene->addItem(nonFixedPointsItem);
//  }
//  else
//  {
//    scene->replaceItem(nonFixedPointsItemIndex, nonFixedPointsItem);
//  }

  scene->itemChanged(index);
  scene->itemChanged(fixedPointsItemIndex);
//  scene->itemChanged(nonFixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_mean_curvature_flow_skeleton_plugin, Polyhedron_demo_mean_curvature_flow_skeleton_plugin)

#include "Polyhedron_demo_mean_curvature_flow_skeleton_plugin.moc"
