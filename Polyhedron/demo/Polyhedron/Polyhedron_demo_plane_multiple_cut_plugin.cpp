#undef NDEBUG // there are assertions to notify degenarete cases which are not handled
              // so keep it open in release mode to see whether there is a prob in algo, or a failure in predefined cases
#define CGAL_PROFILE

#include <QtCore/qglobal.h>
#include <CGAL/AABB_intersections.h>

#include "Messages_interface.h"
#include "Scene_item_with_display_list.h"
#include "Scene_plane_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Plane_multiple_cut_widget.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h> 
#include <CGAL/intersection_of_plane_Polyhedra_3_using_AABB.h>

#include "Polyhedron_type.h"

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

class Polyhedron_demo_plane_multiple_cut_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable() const { return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message);}

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m);
  QList<QAction*> actions() const;

public slots:
  void multiple_plane_action();
  void on_Generate_button_clicked();
  void on_Update_plane_button_clicked();
  void plane_manipulated_frame_modified();
  void plane_destroyed();

private:
  Scene_interface* scene;
  Messages_interface* messages;
  Scene_plane_item* plane_item;
  QAction* actionMultiplePlane;

  QDockWidget* dock_widget;
  Ui::Plane_multiple_cut_widget* ui_widget;
}; // end Polyhedron_demo_plane_multiple_cut_plugin

void Polyhedron_demo_plane_multiple_cut_plugin::init(QMainWindow* mw,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  scene = scene_interface;
  messages = m;
  actionMultiplePlane = new QAction(tr("Multiple plane cut"), mw);
  connect(actionMultiplePlane, SIGNAL(triggered()),
          this, SLOT(multiple_plane_action()));

  dock_widget = new QDockWidget("Multiple plane cut parameters", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  ui_widget = new Ui::Plane_multiple_cut_widget();

  QWidget* qw =new QWidget();
  ui_widget->setupUi(qw);
  dock_widget->setWidget(qw);
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget->Generate_button,  SIGNAL(clicked()), this, SLOT(on_Generate_button_clicked()));   
  connect(ui_widget->Update_plane_button,  SIGNAL(clicked()), this, SLOT(on_Update_plane_button_clicked())); 
}

QList<QAction*> Polyhedron_demo_plane_multiple_cut_plugin::actions() const {
  return QList<QAction*>() << actionMultiplePlane;
}

void Polyhedron_demo_plane_multiple_cut_plugin::multiple_plane_action(){
  if(dock_widget != NULL && !dock_widget->isVisible()) { 
    dock_widget->show(); 

    ///// from cut plugin /////
    plane_item = new Scene_plane_item(scene);
    const Scene_interface::Bbox& bbox = scene->bbox();
    plane_item->setPosition((bbox.xmin + bbox.xmax)/2.f,
      (bbox.ymin+bbox.ymax)/2.f,
      (bbox.zmin+bbox.zmax)/2.f);
    plane_item->setNormal(0., 0., 1.);
    plane_item->setManipulatable(true);
    plane_item->setClonable(false);
    plane_item->setColor(Qt::green);
    plane_item->setName(tr("Cutting plane"));
    connect(plane_item->manipulatedFrame(), SIGNAL(modified()),
      this, SLOT(plane_manipulated_frame_modified()));
    connect(plane_item, SIGNAL(destroyed()),
      this, SLOT(plane_destroyed()));
    scene->addItem(plane_item);

    // set distance_with_planes = bbox_diagona / 30
    double diagonal = std::sqrt(
      CGAL::squared_distanceC3( bbox.xmin, bbox.ymin, bbox.zmin, bbox.xmax, bbox.ymax, bbox.zmax) );
    ui_widget->Distance_with_planes->setText(QString::number(diagonal / 30.0));

    plane_manipulated_frame_modified(); // update text boxes
  }
}

// when manipulated frame of plane is modified, update line-edits
void Polyhedron_demo_plane_multiple_cut_plugin::plane_manipulated_frame_modified() {
  qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  const qglviewer::Vec& pos = mf->position();
  ui_widget->Center_x->setText(QString::number(pos.x));
  ui_widget->Center_y->setText(QString::number(pos.y));
  ui_widget->Center_z->setText(QString::number(pos.z));

  const qglviewer::Vec& base_1 = mf->inverseTransformOf(qglviewer::Vec(1., 0., 0.));
  const qglviewer::Vec& base_2 = mf->inverseTransformOf(qglviewer::Vec(0., 1., 0.));

  ui_widget->Base_1_x->setText(QString::number(base_1.x));
  ui_widget->Base_1_y->setText(QString::number(base_1.y));
  ui_widget->Base_1_z->setText(QString::number(base_1.z));

  ui_widget->Base_2_x->setText(QString::number(base_2.x));
  ui_widget->Base_2_y->setText(QString::number(base_2.y));
  ui_widget->Base_2_z->setText(QString::number(base_2.z));
}

// when Update Plane button is clicked, update manipulated frame of plane with line-edits
void Polyhedron_demo_plane_multiple_cut_plugin::on_Update_plane_button_clicked() {
  qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  //todo: ugly code, check is there any work-around for to_double_ok 
  bool to_double_ok = true;
  double center_x = ui_widget->Center_x->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }
  double center_y = ui_widget->Center_y->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }
  double center_z = ui_widget->Center_z->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }

  bool oldState = mf->blockSignals(true); // dont let it signal, it will invoke plane_manipulated_frame_modified otherwise
  mf->setPosition(center_x, center_y, center_z);
  mf->blockSignals(oldState);

  double base_1_x = ui_widget->Base_1_x->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }
  double base_1_y = ui_widget->Base_1_y->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }
  double base_1_z = ui_widget->Base_1_z->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }

  double base_2_x = ui_widget->Base_2_x->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }
  double base_2_y = ui_widget->Base_2_y->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }
  double base_2_z = ui_widget->Base_2_z->text().toDouble(&to_double_ok);
  if(!to_double_ok) { return; }

  qglviewer::Vec base_1(base_1_x, base_1_y, base_1_z);
  qglviewer::Vec base_2(base_2_x, base_2_y, base_2_z);
  
  qglviewer::Quaternion orientation_from_bases;
  orientation_from_bases.setFromRotatedBasis(base_1, base_2, cross(base_1, base_2));

  oldState = mf->blockSignals(true); // dont let it signal, it will invoke plane_manipulated_frame_modified otherwise
  mf->setOrientation(orientation_from_bases);
  mf->blockSignals(oldState);

  scene->itemChanged(plane_item); // redraw
}

// generate multiple cuts, until any cut does not intersect with bbox
void Polyhedron_demo_plane_multiple_cut_plugin::on_Generate_button_clicked() {

  print_message("Generating multiple cuts...");

  Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  if(!item) { 
    print_message("Error: There is no selected Scene_polyhedron_item!");
    return; 
  }

  // get plane position and normal
  qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  const qglviewer::Vec& pos = mf->position();
  const qglviewer::Vec& n = mf->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));

  // get distance between planes
  bool to_double_ok = true;
  double distance_with_planes = ui_widget->Distance_with_planes->text().toDouble(&to_double_ok);
  if(!to_double_ok) { 
    print_message("Error: Set Distance_with_planes text box!");
    return; 
  }

  // construct a bbox for selected polyhedron
  const Scene_interface::Bbox& bbox = item->bbox();
  CGAL::Bbox_3 cgal_bbox(bbox.xmin, bbox.ymin, bbox.zmin,
    bbox.xmax, bbox.ymax, bbox.zmax);
  Polyhedron* poly = item->polyhedron();

  // continue generating planes while inside bbox
  std::vector<Epic_kernel::Plane_3> planes;
  for(int dir = 1, step = 0; /* */ ; ++step) 
  {
    double distance_norm = (dir * step) * distance_with_planes;
    qglviewer::Vec new_pos = pos + (n*distance_norm); 
    Epic_kernel::Plane_3 plane(n[0], n[1],  n[2], - n * new_pos);

    if(!CGAL::do_intersect(cgal_bbox, plane)) { 
      if(dir == -1) { break; }
      std::reverse(planes.begin(), planes.end());
      dir = -1; // reverse direction
      step = 0; // we should skip the plane itself, and we will when continue cause ++step
      continue;
    }
    planes.push_back(plane);
  }
  print_message(QString("Created %1 cuts inside bbox...").arg(planes.size()));
  
  //CGAL::intersection_of_plane_Polyhedra_3_using_AABB(
  //  *poly, Epic_kernel::Plane_3(0, 0,  1, -0.135),
  //  std::back_inserter(new_polylines_item->polylines)
  //  );

  bool new_polyline_item_for_polylines = ui_widget->newPolylineItemCheckBox->checkState() == Qt::Checked;
  
  if(!new_polyline_item_for_polylines) 
  {
    Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
    QTime time; time.start();
    // call algorithm and fill polylines in polylines_item
    CGAL::intersection_of_plane_Polyhedra_3_using_AABB<Epic_kernel>(
      *poly, planes.begin(), planes.end(),
      std::back_inserter(new_polylines_item->polylines ));
    // set names etc and print timing
    print_message( QString("Processed %1 cuts in %2 ms!").arg(planes.size()).arg(time.elapsed()) );
    new_polylines_item->setName(QString("%1 with %2 cuts").
      arg(item->name()).arg(planes.size()) );
    new_polylines_item->setColor(Qt::green);
    new_polylines_item->setRenderingMode(Wireframe);
    scene->addItem(new_polylines_item);
  }
  else {
    QTime time; time.start();
    std::list<std::vector<Epic_kernel::Point_3> > polylines;
    // call algorithm and fill polylines in polylines_item
    CGAL::intersection_of_plane_Polyhedra_3_using_AABB<Epic_kernel>(
      *poly, planes.begin(), planes.end(),
      std::back_inserter(polylines));
    // set names etc and print timing
    print_message( QString("Processed %1 cuts in %2 ms!").arg(planes.size()).arg(time.elapsed()) );

    int counter = 0;
    for(std::list<std::vector<Epic_kernel::Point_3> >::iterator it = polylines.begin(); it != polylines.end(); ++it, ++counter) {
      Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
      new_polylines_item->polylines.push_back(*it);
      new_polylines_item->setName(QString("%1 with %2 cuts %3").
        arg(item->name()).arg(planes.size()).arg(counter) );
      new_polylines_item->setColor(Qt::green);
      new_polylines_item->setRenderingMode(Wireframe);
      scene->addItem(new_polylines_item);
    }
  }
}

void Polyhedron_demo_plane_multiple_cut_plugin::plane_destroyed() {
  dock_widget->hide(); 
}
Q_EXPORT_PLUGIN2(Polyhedron_demo_plane_multiple_cut_plugin, Polyhedron_demo_plane_multiple_cut_plugin)

#include "Polyhedron_demo_plane_multiple_cut_plugin.moc"
