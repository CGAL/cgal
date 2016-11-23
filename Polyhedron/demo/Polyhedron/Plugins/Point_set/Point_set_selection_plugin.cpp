#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>
#include <QtCore/qglobal.h>
#include <QGLViewer/manipulatedCameraFrame.h>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
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
class Q_DECL_EXPORT Scene_point_set_selection_visualizer : public CGAL::Three::Scene_item
{
  Q_OBJECT

 private:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;
  typedef K::Point_3 Point_3;
  typedef CGAL::Polygon_2<K> Polygon_2;
  typedef std::vector<Point_2> Polyline_2;
  typedef std::vector<Polyline_2> Polylines;
  
  bool rectangle;
  std::vector<Point_2> contour_2d;
  Polylines* polyline;
  Bbox point_set_bbox;
  CGAL::Bbox_2 domain_rectangle;
  Polygon_2 domain_freeform;
  
public:

  Scene_point_set_selection_visualizer(bool rectangle, const Bbox& point_set_bbox)
    : rectangle (rectangle), point_set_bbox (point_set_bbox)
  {
    polyline = new Polylines(0);
    polyline->push_back(Polyline_2());
  }
  ~Scene_point_set_selection_visualizer() {
  }
  bool isFinite() const { return true; }
  bool isEmpty() const { return poly().empty(); }
  void compute_bbox() const {
    _bbox = point_set_bbox;
  }
  Scene_point_set_selection_visualizer* clone() const {
    return 0;
  }
  QString toolTip() const {
    return tr("%1").arg(name());
  }

  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Points);
  }
  /*
   * We use drawPoints because it is the last drawing function to be called
   * This way, as this item is always the last in the list, it is always the last
   * thing to be drawn. It allow us to safely call glClear(GL_DEPTH_BUFFER_BIT)
   * and have the selection lines always on top of the other items.
   */
  void drawPoints(CGAL::Three::Viewer_interface* viewer) const {

    QPainter *painter = viewer->getPainter();
    QPen pen;
    pen.setColor(QColor(Qt::green));
    pen.setWidth(5);

    painter->setPen(pen);
    glClear(GL_DEPTH_BUFFER_BIT);
    for(std::size_t i=0; i<polyline->size(); ++i)
    {
      Polyline_2 poly = (*polyline)[i];
      if(!poly.empty())
        for(std::size_t j=0; j<poly.size()-1; ++j)
        {
          viewer->getPainter()->drawLine(poly[j].x(), poly[j].y(), poly[j+1].x(), poly[j+1].y());
        }
    }
  }

  Polyline_2& poly() const
  { return polyline->front(); }
  
  bool update_polyline () const
  {
    if (contour_2d.size() < 2 ||
        (!(poly().empty()) && contour_2d.back () == poly().back()))
      return false;
    
    if (rectangle)
      {
	poly().clear();
	
        poly().push_back ( Point_2 (domain_rectangle.xmin(),
                                domain_rectangle.ymin()));
        poly().push_back ( Point_2 (domain_rectangle.xmax(),
                                domain_rectangle.ymin()));
        poly().push_back ( Point_2 (domain_rectangle.xmax(),
                                domain_rectangle.ymax()));
        poly().push_back ( Point_2 (domain_rectangle.xmin(),
                                domain_rectangle.ymax()));
        poly().push_back ( Point_2 (domain_rectangle.xmin(),
                                                domain_rectangle.ymin()));

      }
    else
      {
        if (!(poly().empty()) && contour_2d.back () == poly().back())
	  return false;

	poly().clear();

	for (unsigned int i = 0; i < contour_2d.size (); ++ i)
          poly().push_back (contour_2d[i]);
      }
    return true;
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
/*
 * domain_freeform.has_on_bounded_side() requires the polygon to be simple, which is never the case.
 * However, it works very well even if the polygon is not simple, so we use this instead to avoid
 * the cgal_assertion on is_simple().*/


        if (CGAL::bounded_side_2(domain_freeform.container().begin(),
                                 domain_freeform.container().end(),
                                 Point_2(p.x, p.y),
                                 domain_freeform.traits_member())  == CGAL::ON_BOUNDED_SIDE)
          return true;
      }
    return false;
  }


}; // end class Scene_point_set_selection_visualizer
///////////////////////////////////////////////////////////////////////////////////////////////////

using namespace CGAL::Three;
class Polyhedron_demo_point_set_selection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction*) const { 
      return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointSetSelection; }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionPointSetSelection = new QAction(tr("Selection"), mw);
    connect(actionPointSetSelection, SIGNAL(triggered()), this, SLOT(selection_action()));

    dock_widget = new QDockWidget("Point Set Selection", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);

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
  virtual void closure()
  {
    dock_widget->hide();
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
            visualizer->setRenderingMode (Points);
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
    Scene_points_with_normal_item* point_set_item = getSelectedItem<Scene_points_with_normal_item>();
    if(!point_set_item)
      {
	print_message("Error: no point set selected!");
	return; 
      }

    Point_set* points = point_set_item->point_set();

    if (selection_mode == 0) // New selection
      points->unselect_all();
    
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();

    std::vector<Point_set::Index> unselected, selected;
    
    for(Point_set::iterator it = points->begin ();
	it != points->end(); ++ it)
      {
	bool already_selected = points->is_selected (it);

        const Kernel::Point_3 p = points->point (*it);
	qglviewer::Vec vp (p.x (), p.y (), p.z ());
	qglviewer::Vec vsp = camera->projectedCoordinatesOf (vp);
	    
	bool now_selected = visualizer->is_selected (vsp);

	// NEW INTERSECTION
	if (selection_mode == 0)
	  {
	    if (now_selected)
	      selected.push_back (*it);
	    else
	      unselected.push_back (*it);
	  }
	// UNION
	else if (selection_mode == 1)
	  {
	    if (already_selected || now_selected)
	      selected.push_back (*it);
	    else
	      unselected.push_back (*it);
	  }
	// INTERSECTION
	//  * Unselect point if it was selected and is not anymore
	else if (selection_mode == 2)
	  {
	    if (already_selected && now_selected)
	      selected.push_back (*it);
	    else
	      unselected.push_back (*it);
	  }
	// DIFFERENCE
	//  * Unselect point if it was selected and is now selected
	else if (selection_mode == 3)
	  {
	    if (already_selected && !now_selected)
	      selected.push_back (*it);
	    else
	      unselected.push_back (*it);
	  }

      }

    for (std::size_t i = 0; i < unselected.size(); ++ i)
      *(points->begin() + i) = unselected[i];
    for (std::size_t i = 0; i < selected.size(); ++ i)
      *(points->begin() + (unselected.size() + i)) = selected[i];

    if (selected.empty ())
      {
	points->unselect_all();
      }
    else
      {
	points->set_first_selected
	  (points->begin() + unselected.size());
      } 
    point_set_item->invalidateOpenGLBuffers();
    point_set_item->itemChanged();
  }

  


public Q_SLOTS:
  void selection_action() { 
    dock_widget->show();
    dock_widget->raise();
  }
  
  // Select all
  void on_Select_all_button_clicked() {
    Scene_points_with_normal_item* point_set_item = getSelectedItem<Scene_points_with_normal_item>();
    if(!point_set_item)
      {
	print_message("Error: no point set selected!");
	return; 
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    point_set_item->selectAll();
    QApplication::restoreOverrideCursor();
  }
  
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    point_set_item->resetSelection();
    QApplication::restoreOverrideCursor();
  }

  void on_Erase_selected_points_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    point_set_item->deleteSelection();
    QApplication::restoreOverrideCursor();
  }

  void on_Invert_selection_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    point_set_item->invertSelection();
    QApplication::restoreOverrideCursor();
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
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item();
    new_item->setName(QString("%1 (selected points)").arg(point_set_item->name()));
    if (point_set_item->has_normals())
      new_item->point_set()->add_normal_map();
    if (point_set_item->point_set()->has_colors())
      {
        new_item->point_set()->add_property_map<unsigned char>("red", 0);
        new_item->point_set()->add_property_map<unsigned char>("green", 0);
        new_item->point_set()->add_property_map<unsigned char>("blue", 0);
        new_item->point_set()->check_colors(); 
      }
    
    new_item->setColor(point_set_item->color());
    new_item->setRenderingMode(point_set_item->renderingMode());
    new_item->setVisible(point_set_item->visible());

    typedef Point_set_3<Kernel> Point_set;
    for(Point_set::iterator it = point_set_item->point_set()->begin ();
	it != point_set_item->point_set()->end(); ++ it) {
      if (point_set_item->point_set()->is_selected (it))
        {
          Point_set::iterator new_point =
            new_item->point_set()->insert(point_set_item->point_set()->point(*it));
          if (point_set_item->has_normals())
            new_item->point_set()->normal(*new_point) = point_set_item->point_set()->normal(*it);
          if (point_set_item->point_set()->has_colors())
            {
              new_item->point_set()->red(*new_point) = point_set_item->point_set()->red(*it);
              new_item->point_set()->green(*new_point) = point_set_item->point_set()->green(*it);
              new_item->point_set()->blue(*new_point) = point_set_item->point_set()->blue(*it);
            }
        }
    }
    new_item->resetSelection();
    new_item->invalidateOpenGLBuffers();

    scene->addItem(new_item);
    QApplication::restoreOverrideCursor();
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

#include "Point_set_selection_plugin.moc"
