#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>
#include <QtCore/qglobal.h>
#include <QGLViewer/manipulatedCameraFrame.h>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include "Scene_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "ui_Point_set_selection_widget.h"
#include "Point_set_3.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>

#include <map>
#include <fstream>


// Class for visualizing selection 
// provides mouse selection functionality
class Q_DECL_EXPORT Scene_point_set_selection_visualizer : public Scene_item
{
  Q_OBJECT

 private:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;
  typedef K::Point_3 Point_3;
  typedef CGAL::Polygon_2<K> Polygon_2;
  
  bool rectangle;
  std::vector<Point_2> contour_2d;
  Scene_polylines_item* polyline;
  Bbox point_set_bbox;
  CGAL::Bbox_2 domain_rectangle;
  Polygon_2 domain_freeform;
  
public:

  Scene_point_set_selection_visualizer(bool rectangle, const Bbox& point_set_bbox)
    : rectangle (rectangle), point_set_bbox (point_set_bbox)
  {
    polyline = new Scene_polylines_item();
    polyline->setRenderingMode (Wireframe);
    polyline->setVisible (true);
    polyline->polylines.push_back (Scene_polylines_item::Polyline());
  }
  ~Scene_point_set_selection_visualizer() {
  }
  bool isFinite() const { return true; }
  bool isEmpty() const { return poly().empty(); }
  Bbox bbox() const {
    return point_set_bbox;
  }
  Scene_point_set_selection_visualizer* clone() const {
    return 0;
  }
  QString toolTip() const {
    return tr("%1").arg(name());
  }

  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe);
  }
  
  void draw_edges(Viewer_interface* viewer) const {
    viewer->glLineWidth(3.f);
    polyline->setRbgColor(0, 255, 0); 

    polyline->draw_edges(viewer);
  }

  Scene_polylines_item::Polyline& poly() const
  { return polyline->polylines.front(); }
  
  bool update_polyline () const
  {
    if (contour_2d.size() < 2 ||
	(!(poly().empty()) && scene_point (contour_2d.back ()) == poly().back()))
      return false;
    
    if (rectangle)
      {
	poly().clear();
	
	poly().push_back (scene_point (Point_2 (domain_rectangle.xmin(),
						domain_rectangle.ymin())));
	poly().push_back (scene_point (Point_2 (domain_rectangle.xmax(),
						domain_rectangle.ymin())));
	poly().push_back (scene_point (Point_2 (domain_rectangle.xmax(),
						domain_rectangle.ymax())));
	poly().push_back (scene_point (Point_2 (domain_rectangle.xmin(),
						domain_rectangle.ymax())));
	poly().push_back (scene_point (Point_2 (domain_rectangle.xmin(),
						domain_rectangle.ymin())));

      }
    else
      {
	if (!(poly().empty()) && scene_point (contour_2d.back ()) == poly().back())
	  return false;

	poly().clear();

	for (unsigned int i = 0; i < contour_2d.size (); ++ i)
	  poly().push_back (scene_point (contour_2d[i]));
      }
    
    return true;
  }

  Point_3 scene_point (const Point_2& p) const
  {
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();
    qglviewer::Vec vp (p.x(), p.y(), 0.1);
    qglviewer::Vec vsp = camera->unprojectedCoordinatesOf (vp);
    
    return Point_3 (vsp.x, vsp.y, vsp.z);
  }


  
  void sample_mouse_path()
  {
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
    
    if (rectangle && contour_2d.size () == 2)
      {
	contour_2d[1] = Point_2 (p.x (), p.y ());
	domain_rectangle = CGAL::bbox_2 (contour_2d.begin (), contour_2d.end ());
      }
    else
      contour_2d.push_back (Point_2 (p.x (), p.y ()));

    if (update_polyline ())
      {
	Q_EMIT itemChanged();
      }
  }

  void apply_path()
  {
    update_polyline ();
    domain_rectangle = CGAL::bbox_2 (contour_2d.begin (), contour_2d.end ());    
    if (!rectangle)
      domain_freeform = Polygon_2 (contour_2d.begin (), contour_2d.end ());
  }

  bool is_selected (qglviewer::Vec& p)
  {
    if (domain_rectangle.xmin () < p.x &&
	p.x < domain_rectangle.xmax () &&
	domain_rectangle.ymin () < p.y &&
	p.y < domain_rectangle.ymax ())
      {
	if (rectangle)
	  return true;
	
	if (domain_freeform.has_on_bounded_side (Point_2 (p.x, p.y)))
	  return true;
      }
    return false;
  }


}; // end class Scene_point_set_selection_visualizer
///////////////////////////////////////////////////////////////////////////////////////////////////


class Polyhedron_demo_point_set_selection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction*) const { 
      return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointSetSelection; }
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionPointSetSelection = new QAction(tr("Selection"), mw);
    connect(actionPointSetSelection, SIGNAL(triggered()), this, SLOT(selection_action()));

    dock_widget = new QDockWidget("Point Set Selection", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    add_dock_widget(dock_widget);

    connect(ui_widget.Selection_tool_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_tool_combo_box_changed(int)));
    connect(ui_widget.Selection_mode_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_mode_combo_box_changed(int)));
    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Invert_selection_button,  SIGNAL(clicked()), this, SLOT(on_Invert_selection_button_clicked()));
    connect(ui_widget.Erase_selected_points_button,  SIGNAL(clicked()), this, SLOT(on_Erase_selected_points_button_clicked()));
    connect(ui_widget.Create_point_set_item_button, SIGNAL(clicked()), this, SLOT(on_Create_point_set_item_button_clicked()));

    rectangle = true;
    selection_mode = 0;
    visualizer = NULL;
    shift_pressing = false;
    
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    mainWindow->installEventFilter(this);


  }


protected:

  bool eventFilter(QObject *, QEvent *event) {
    if (dock_widget->isHidden() || !(dock_widget->isActiveWindow()))
      return false;
    
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      return false; 
    }
    int item_id = scene->item_id (point_set_item);
      
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)  {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
    }

    // mouse events
    if(shift_pressing && event->type() == QEvent::MouseButtonPress)
      {
	QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
	// Start selection
	if (mouseEvent->button() == Qt::LeftButton && !visualizer)
	  {
	    QApplication::setOverrideCursor(Qt::CrossCursor);
	    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
	    if (viewer->camera()->frame()->isSpinning())
	      viewer->camera()->frame()->stopSpinning();
	
	    visualizer = new Scene_point_set_selection_visualizer(rectangle,
								  point_set_item->bbox());
	    visualizer->setName(tr("Point set selection visualizer"));
	    visualizer->setRenderingMode (Wireframe);
	    visualizer->setVisible (true);

	    // Hack to prevent camera for "jumping" when creating new item
	    scene->addItem(visualizer);
	
	    scene->setSelectedItem(item_id);
	    visualizer->sample_mouse_path();
	    return true;
	  }
	// Cancel selection
	else if (mouseEvent->button() == Qt::RightButton && visualizer)
	  {

	    scene->erase( scene->item_id(visualizer) );
	    scene->setSelectedItem(item_id);
	    visualizer = NULL;
	    QApplication::restoreOverrideCursor();
	    return true;
	  }

      }
    // End selection
    else if (event->type() == QEvent::MouseButtonRelease && visualizer)
      {
	visualizer->apply_path();
	select_points();
	scene->erase( scene->item_id(visualizer) );
	scene->setSelectedItem(item_id);
	visualizer = NULL;
	QApplication::restoreOverrideCursor();
	return true;
      }
    // Update selection
    else if (event->type() == QEvent::MouseMove && visualizer)
      {
	visualizer->sample_mouse_path();
	return true;
      }

    return false;
  }

  void select_points()
  {
    Scene_points_with_normal_item* point_set_item = get_selected_item<Scene_points_with_normal_item>();
    if(!point_set_item)
      {
	print_message("Error: no point set selected!");
	return; 
      }

    if (selection_mode == 0) // New selection
      {
	for(Point_set::iterator it = point_set_item->point_set()->begin ();
	    it != point_set_item->point_set()->end(); ++ it)
	  point_set_item->point_set()->select(&*it,false);
      }
    
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();
    
    for(Point_set::iterator it = point_set_item->point_set()->begin ();
	it != point_set_item->point_set()->end(); ++ it)
      {
	bool already_selected = it->is_selected ();
	    
	qglviewer::Vec vp (it->x (), it->y (), it->z ());
	qglviewer::Vec vsp = camera->projectedCoordinatesOf (vp);
	    
	bool now_selected = visualizer->is_selected (vsp);

	// NEW INTERSECTION or UNION
	//  * Select point if it is now selected
	if (selection_mode < 2 && now_selected)
	  point_set_item->point_set()->select(&*it,true);
	// INTERSECTION
	//  * Unselect point if it was selected and is not anymore
	else if (selection_mode == 2)
	  {
	    if (already_selected && !now_selected)
	      point_set_item->point_set()->select(&*it,false);
	  }
	// DIFFERENCE
	//  * Unselect point if it was selected and is now selected
	else if (selection_mode == 3)
	  {
	    if (already_selected && now_selected)
	      point_set_item->point_set()->select(&*it,false);
	  }

      }

    point_set_item->invalidate_buffers();
  }

  


public Q_SLOTS:
  void selection_action() { 
    dock_widget->show();
    dock_widget->raise();
  }
  
  // Select all
  void on_Select_all_button_clicked() {
    Scene_points_with_normal_item* point_set_item = get_selected_item<Scene_points_with_normal_item>();
    if(!point_set_item)
      {
	print_message("Error: no point set selected!");
	return; 
      }

    point_set_item->selectAll();
  }
  
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }

    point_set_item->resetSelection();
  }

  void on_Erase_selected_points_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }

    point_set_item->deleteSelection();
  }

  void on_Invert_selection_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }

    point_set_item->invertSelection();
  }

  void on_Create_point_set_item_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }
    if(point_set_item->isSelectionEmpty ()) {
      print_message("Error: there is no selected point in point set item!");
      return;
    }
    
    Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item();
    new_item->setName(QString("%1 (selected points)").arg(point_set_item->name()));
    new_item->set_has_normals (point_set_item->has_normals());
    new_item->setColor(point_set_item->color());
    new_item->setRenderingMode(point_set_item->renderingMode());
    new_item->setVisible(point_set_item->visible());

    typedef Point_set_3<Kernel> Point_set;
    for(Point_set::iterator it = point_set_item->point_set()->begin ();
	it != point_set_item->point_set()->end(); ++ it) {
      if (it->is_selected ())
	new_item->point_set()->push_back(*it);
    }
    new_item->resetSelection();
    new_item->invalidate_buffers();

    scene->addItem(new_item);
 }

  void on_Selection_tool_combo_box_changed (int index)
  {
    rectangle = (index == 0);
  }

  void on_Selection_mode_combo_box_changed (int index)
  {
    selection_mode = index;
  }


private:
  Messages_interface* messages;
  QAction* actionPointSetSelection;

  QDockWidget* dock_widget;
  Ui::PointSetSelection ui_widget;
  bool rectangle;
  int selection_mode;
  Scene_point_set_selection_visualizer* visualizer;
  bool shift_pressing;

}; // end Polyhedron_demo_point_set_selection_plugin

//Q_EXPORT_PLUGIN2(Polyhedron_demo_point_set_selection_plugin, Polyhedron_demo_point_set_selection_plugin)

#include "Polyhedron_demo_point_set_selection_plugin.moc"
