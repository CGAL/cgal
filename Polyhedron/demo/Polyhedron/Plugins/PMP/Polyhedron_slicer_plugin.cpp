#include <QtCore/qglobal.h>
#include <CGAL/intersections.h>

#include "Messages_interface.h"
#include "Scene_plane_item.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"
#include "Scene.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include "ui_Polyhedron_slicer_widget.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_slicer.h>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
using namespace CGAL::Three;
class Polyhedron_demo_polyhedron_slicer_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void print_message(QString message) { CGAL::Three::Three::information(message);}

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m);
  virtual void closure()
  {
    dock_widget->hide();
  }
  QList<QAction*> actions() const;
  bool get_base_1_2(double bases[6]) {
    bool oks[6];
    bases[0] = ui_widget.Base_1_x->text().toDouble(&oks[0]);
    bases[1] = ui_widget.Base_1_y->text().toDouble(&oks[1]);
    bases[2] = ui_widget.Base_1_z->text().toDouble(&oks[2]);

    bases[3] = ui_widget.Base_2_x->text().toDouble(&oks[3]);
    bases[4] = ui_widget.Base_2_y->text().toDouble(&oks[4]);
    bases[5] = ui_widget.Base_2_z->text().toDouble(&oks[5]);

    bool total_ok = true;
    for(int i = 0; i < 6; ++i && total_ok) { total_ok &= oks[i];}
    return total_ok;
  }
public Q_SLOTS:
  void slicer_widget_action();
  void on_Generate_button_clicked();
  bool on_Update_plane_button_clicked();
  void plane_manipulated_frame_modified();
  void item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item);
  void dock_widget_closed();

protected:
  bool eventFilter(QObject *, QEvent *event) {
    if(event->type() == QEvent::Close) {
      dock_widget_closed();
    }
    return false;
  }

private:
  Messages_interface* messages;
  Scene_plane_item* plane_item;
  QAction* actionSlicerWidget;

  QDockWidget* dock_widget;
  Ui::Polyhedron_slicer ui_widget;

  template <typename TriangleMesh>
  void intersection_of_plane_Polyhedra_3_using_AABB_wrapper(TriangleMesh& mesh,
    const std::vector<Epic_kernel::Plane_3>& planes,
    const std::vector<CGAL::qglviewer::Vec>& plane_positions,
    std::list<std::vector<Epic_kernel::Point_3> >& polylines);

}; // end Polyhedron_demo_polyhedron_slicer_plugin

void Polyhedron_demo_polyhedron_slicer_plugin::init(QMainWindow* mainWindow,
                                      CGAL::Three::Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  mw = mainWindow;
  scene = scene_interface;
  messages = m;
  plane_item = NULL;

  actionSlicerWidget = new QAction(tr("Polyhedron Slicer"), mw);
  actionSlicerWidget->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionSlicerWidget, SIGNAL(triggered()), this, SLOT(slicer_widget_action()));

  dock_widget = new QDockWidget("Polyhedron Slicer", mw);
  dock_widget->setVisible(false);
  dock_widget->installEventFilter(this);
  ui_widget.setupUi(dock_widget);

  addDockWidget(dock_widget);

  connect(ui_widget.Generate_button,  SIGNAL(clicked()), this, SLOT(on_Generate_button_clicked()));
  connect(ui_widget.Update_plane_button,  SIGNAL(clicked()), this, SLOT(on_Update_plane_button_clicked()));
}

QList<QAction*> Polyhedron_demo_polyhedron_slicer_plugin::actions() const {
  return QList<QAction*>() << actionSlicerWidget;
}

void Polyhedron_demo_polyhedron_slicer_plugin::slicer_widget_action(){
  if(dock_widget->isVisible()) { return; }
  dock_widget->show();
  dock_widget->raise();
  ///// from cut plugin /////
  CGAL_assertion(plane_item == NULL);

  plane_item = new Scene_plane_item(scene);
  const CGAL::Three::Scene_interface::Bbox& bbox = scene->bbox();
  plane_item->setPosition((bbox.xmin() + bbox.xmax())/2.f,
    (bbox.ymin()+bbox.ymax())/2.f,
    (bbox.zmin()+bbox.zmax())/2.f);
  plane_item->setNormal(0., 0., 1.);
  plane_item->setManipulatable(true);
  plane_item->setClonable(false);
  plane_item->setColor(Qt::green);
  plane_item->setName(tr("Cutting plane"));
  connect(plane_item->manipulatedFrame(), SIGNAL(modified()),
    this, SLOT(plane_manipulated_frame_modified()));

  if(Scene* scene_casted = dynamic_cast<Scene*>(scene))
  { connect(scene_casted, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this, SLOT(item_about_to_be_destroyed(CGAL::Three::Scene_item*))); }
  scene->addItem(plane_item);

  // set distance_with_planes = bbox_diagona / 30
  double diagonal = std::sqrt(
    CGAL::squared_distanceC3( bbox.xmin(), bbox.ymin(), bbox.zmin(), bbox.xmax(), bbox.ymax(), bbox.zmax()) );
  ui_widget.Distance_with_planes->setText(QString::number(diagonal / 30.0));

  plane_manipulated_frame_modified(); // update text boxes
}

// when manipulated frame of plane is modified, update line-edits
void Polyhedron_demo_polyhedron_slicer_plugin::plane_manipulated_frame_modified() {
  CGAL::qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  const CGAL::qglviewer::Vec& pos = mf->position();
  ui_widget.Center_x->setText(QString::number(pos.x));
  ui_widget.Center_y->setText(QString::number(pos.y));
  ui_widget.Center_z->setText(QString::number(pos.z));

  const CGAL::qglviewer::Vec& base_1 = mf->inverseTransformOf(CGAL::qglviewer::Vec(1., 0., 0.));
  const CGAL::qglviewer::Vec& base_2 = mf->inverseTransformOf(CGAL::qglviewer::Vec(0., 1., 0.));

  ui_widget.Base_1_x->setText(QString::number(base_1.x));
  ui_widget.Base_1_y->setText(QString::number(base_1.y));
  ui_widget.Base_1_z->setText(QString::number(base_1.z));

  ui_widget.Base_2_x->setText(QString::number(base_2.x));
  ui_widget.Base_2_y->setText(QString::number(base_2.y));
  ui_widget.Base_2_z->setText(QString::number(base_2.z));
}

// when Update Plane button is clicked, update manipulated frame of plane with line-edits
bool Polyhedron_demo_polyhedron_slicer_plugin::on_Update_plane_button_clicked() {
  CGAL::qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  // get center
  bool ok_1 = true, ok_2 = true, ok_3 = true;
  double center_x = ui_widget.Center_x->text().toDouble(&ok_1);
  double center_y = ui_widget.Center_y->text().toDouble(&ok_2);
  double center_z = ui_widget.Center_z->text().toDouble(&ok_3);
  if(!ok_1 || !ok_2 || !ok_3)
  { print_message("Error: center coordinates not convertible to double."); return false; }

  // set center
  bool oldState = mf->blockSignals(true); // dont let it signal, it will invoke plane_manipulated_frame_modified otherwise
  mf->setPosition(center_x, center_y, center_z);
  mf->blockSignals(oldState);

  // get base 1 and base 2
  double bases[6];
  if(!get_base_1_2(bases))
  { print_message("Error: Base-1, Base-2 coordinates not convertible to double."); return false; }

  // compute other axis
  CGAL::qglviewer::Vec base_1(bases[0], bases[1], bases[2]);
  CGAL::qglviewer::Vec base_2(bases[3], bases[4], bases[5]);
  CGAL::qglviewer::Vec other = cross(base_1, base_2);
  if(other.norm() == 0.0) { print_message("Error: collinear base vectors are not accepted!"); return false; }

  // set orientation
  CGAL::qglviewer::Quaternion orientation_from_bases;
  orientation_from_bases.setFromRotatedBasis(base_1, base_2, other);

  oldState = mf->blockSignals(true); // dont let it signal, it will invoke plane_manipulated_frame_modified otherwise
  mf->setOrientation(orientation_from_bases);
  mf->blockSignals(oldState);

  scene->itemChanged(plane_item); // redraw
  return true;
}

// generate multiple cuts, until any cut does not intersect with bbox
void Polyhedron_demo_polyhedron_slicer_plugin::on_Generate_button_clicked()
{
  Scene_surface_mesh_item* sm_item = getSelectedItem<Scene_surface_mesh_item>();
  if(! sm_item) {
    print_message("Error: There is no selected Scene_surface_mesh_item!");
    return;
  }
  QString item_name = sm_item->name();

  if(!on_Update_plane_button_clicked()) { return; }
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  QApplication::setOverrideCursor(Qt::WaitCursor);
  // get plane position and normal
  CGAL::qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  const CGAL::qglviewer::Vec& pos = mf->position()-offset;
  // WARNING: due to fp arithmetic (setting quaternion based orientation from base vectors then getting plane normal back from this orientation)
  // for base vectors like: 1,0,0 - 0,1,0 we might not have exact corresponding normal vector.
  // So not using below normal but construct plane directly from bases from text boxes
  const CGAL::qglviewer::Vec& n = mf->inverseTransformOf(CGAL::qglviewer::Vec(0.f, 0.f, 1.f));

  // get bases
  double bases[6];
  get_base_1_2(bases); // no need to check since we call on_Update_plane_button_clicked
  Epic_kernel::Vector_3 base_1(bases[0], bases[1], bases[2]);
  Epic_kernel::Vector_3 base_2(bases[3], bases[4], bases[5]);
  const Epic_kernel::Vector_3 normal = CGAL::cross_product(base_1, base_2);

  // get distance between planes
  bool to_double_ok = true;
  double distance_with_planes = ui_widget.Distance_with_planes->text().toDouble(&to_double_ok);
  if(!to_double_ok) {
    print_message("Error: Set Distance_with_planes text box!");
    QApplication::restoreOverrideCursor();
    return;
  }

  // construct a bbox for selected polyhedron
  const CGAL::Three::Scene_interface::Bbox& bbox = sm_item->bbox();
  CGAL::Bbox_3 cgal_bbox(bbox.xmin(), bbox.ymin(), bbox.zmin(),
    bbox.xmax(), bbox.ymax(), bbox.zmax());
  SMesh* smesh = sm_item->polyhedron();

  // continue generating planes while inside bbox
  std::vector<Epic_kernel::Plane_3> planes;
  std::vector<CGAL::qglviewer::Vec> plane_positions;

  for(int dir = 1, step = 0; /* */ ; ++step)
  {
    double distance_norm = (dir * step) * distance_with_planes;
    CGAL::qglviewer::Vec new_pos = pos + (n*distance_norm);

    //Epic_kernel::Plane_3 plane(n[0], n[1],  n[2], - n * new_pos);
    Epic_kernel::Point_3 new_pos_cgal(new_pos[0], new_pos[1], new_pos[2]);
    Epic_kernel::Plane_3 plane(new_pos_cgal, normal);

    if(!CGAL::do_intersect(cgal_bbox, plane)) {
      if(dir == -1) { break; }
      std::reverse(planes.begin(), planes.end());
      std::reverse(plane_positions.begin(), plane_positions.end());
      dir = -1; // reverse direction
      step = 0; // we should skip the plane itself, and we will when continue cause ++step
      continue;
    }
    planes.push_back(plane);
    plane_positions.push_back(new_pos);
  }
  print_message(QString("Created %1 cuts inside bbox...").arg(planes.size()));

  bool new_polyline_item_for_polylines = ui_widget.newPolylineItemCheckBox->checkState() == Qt::Checked;
  if(!new_polyline_item_for_polylines)
  {
    Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
    QElapsedTimer time; time.start();
    // call algorithm and fill polylines in polylines_item
    intersection_of_plane_Polyhedra_3_using_AABB_wrapper(*smesh, planes, plane_positions, new_polylines_item->polylines);
    // set names etc and print timing
    print_message( QString("Done: processed %1 cuts - generated %2 polylines in %3 ms!").
      arg(planes.size()).arg(new_polylines_item->polylines.size()).arg(time.elapsed()) );

    new_polylines_item->setName(QString("%1 with %2 cuts").
      arg(item_name).arg(planes.size()) );
    new_polylines_item->setColor(Qt::green);
    new_polylines_item->setRenderingMode(Wireframe);
    scene->addItem(new_polylines_item);
  }
  else {
    QElapsedTimer time; time.start();
    std::list<std::vector<Epic_kernel::Point_3> > polylines;
    // call algorithm and fill polylines in polylines_item
    intersection_of_plane_Polyhedra_3_using_AABB_wrapper(*smesh, planes, plane_positions, polylines);
    // set names etc and print timing
    print_message( QString("Done: processed %1 cuts - generated %2 polylines in %3 ms!").
      arg(planes.size()).arg(polylines.size()).arg(time.elapsed()) );

    int counter = 0;
    for(std::list<std::vector<Epic_kernel::Point_3> >::iterator it = polylines.begin(); it != polylines.end(); ++it, ++counter) {
      Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
      new_polylines_item->polylines.push_back(*it);
      new_polylines_item->setName(QString("%1 with %2 cuts %3").
        arg(item_name).arg(planes.size()).arg(counter) );
      new_polylines_item->setColor(Qt::green);
      new_polylines_item->setRenderingMode(Wireframe);
      scene->addItem(new_polylines_item);
      new_polylines_item->invalidateOpenGLBuffers();
    }
  }
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_polyhedron_slicer_plugin::item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item) {
  if(plane_item == NULL) { return; }// which means this plugin erased plane_item
  Scene_plane_item* destroyed_plane = qobject_cast<Scene_plane_item*>(scene_item);
  if(destroyed_plane && destroyed_plane == plane_item) {
    plane_item = NULL;
    dock_widget->hide();
  }
}

void Polyhedron_demo_polyhedron_slicer_plugin::dock_widget_closed() {
  CGAL_assertion(plane_item != NULL);
  CGAL::Three::Scene_interface::Item_id id = scene->item_id(plane_item);
  plane_item = NULL;
  scene->erase(id);
}
// this function assumes 'planes' are parallel
template <typename TriangleMesh>
void Polyhedron_demo_polyhedron_slicer_plugin::intersection_of_plane_Polyhedra_3_using_AABB_wrapper(
  TriangleMesh& poly,
  const std::vector<Epic_kernel::Plane_3>& planes,
  const std::vector<CGAL::qglviewer::Vec>& plane_positions,
  std::list<std::vector<Epic_kernel::Point_3> >& polylines)
{
  CGAL::Polygon_mesh_slicer<TriangleMesh, Epic_kernel> slicer(poly);
  std::vector<CGAL::qglviewer::Vec>::const_iterator plane_position_it = plane_positions.begin();
  for(std::vector<Epic_kernel::Plane_3>::const_iterator plane_it = planes.begin(); plane_it != planes.end(); ++plane_it, ++plane_position_it)
    slicer(*plane_it, std::front_inserter(polylines));

}

#include "Polyhedron_slicer_plugin.moc"
