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
#include <Plugins/PCA/Scene_edit_box_item.h>
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
class Q_DECL_EXPORT Scene_point_set_selection_visualizer
{

 private:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;
  typedef K::Point_3 Point_3;
  typedef CGAL::Polygon_2<K> Polygon_2;
  typedef std::vector<Point_2> Polyline_2;
  typedef std::vector<Polyline_2> Polylines;
  typedef Scene_item::Bbox Bbox;
  
  bool rectangle;
  std::vector<Point_2> contour_2d;
  Polylines* polyline;
  Scene_item::Bbox point_set_bbox;
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

  void render(QImage& image) const {

    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
    QPainter *painter = viewer->getPainter();
    QPen pen;
    pen.setColor(QColor(Qt::green));
    pen.setWidth(5);


    painter->begin(viewer);
    painter->drawImage(QPoint(0,0), image);
    painter->setPen(pen);
    for(std::size_t i=0; i<polyline->size(); ++i)
    {
      Polyline_2 poly = (*polyline)[i];
      if(!poly.empty())
        for(std::size_t j=0; j<poly.size()-1; ++j)
        {
          painter->drawLine(poly[j].x(), poly[j].y(), poly[j+1].x(), poly[j+1].y());
        }
    }
    painter->end();
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

  
  void sample_mouse_path(QImage& image)
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

        render(image);
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

    ui_widget.Add_to_selection_Button->setVisible(false);
    connect(ui_widget.Selection_tool_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_tool_combo_box_changed(int)));
    connect(ui_widget.Selection_mode_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_mode_combo_box_changed(int)));
    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Invert_selection_button,  SIGNAL(clicked()), this, SLOT(on_Invert_selection_button_clicked()));
    connect(ui_widget.Erase_selected_points_button,  SIGNAL(clicked()), this, SLOT(on_Erase_selected_points_button_clicked()));
    connect(ui_widget.Create_point_set_item_button, SIGNAL(clicked()), this, SLOT(on_Create_point_set_item_button_clicked()));
    connect(ui_widget.Add_to_selection_Button, SIGNAL(clicked()), this, SLOT(select_points()));
    rectangle = true;
    selection_mode = 0;
    visualizer = NULL;
    edit_box = NULL;
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
    static QImage background;
    if (dock_widget->isHidden() || !(dock_widget->isActiveWindow()) || ui_widget.Selection_tool_combo_box->currentIndex()==2)
      return false;
    
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      return false; 
    }
      
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
    {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
#if QGLVIEWER_VERSION >= 0x020700
      background = static_cast<CGAL::Three::Viewer_interface*>(*QGLViewer::QGLViewerPool().begin())->grabFramebuffer();
#else
      background = static_cast<CGAL::Three::Viewer_interface*>(*QGLViewer::QGLViewerPool().begin())->grabFrameBuffer();

#endif
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

            visualizer->sample_mouse_path(background);
            return true;
	  }
	// Cancel selection
	else if (mouseEvent->button() == Qt::RightButton && visualizer)
	  {
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
	visualizer = NULL;
	QApplication::restoreOverrideCursor();
	return true;
      }
    // Update selection
    else if (event->type() == QEvent::MouseMove && visualizer)
      {
        visualizer->sample_mouse_path(background);
	return true;
      }

    return false;
  }

protected Q_SLOTS:
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
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(viewer)->offset();

    qglviewer::Camera* camera = viewer->camera();

    std::vector<Point_set::Index> unselected, selected;
    for(Point_set::iterator it = points->begin ();
	it != points->end(); ++ it)
      {
	bool already_selected = points->is_selected (it);

        const Kernel::Point_3 p = points->point (*it);
        qglviewer::Vec vp (p.x (), p.y (), p.z ());
        qglviewer::Vec vsp = camera->projectedCoordinatesOf (vp + offset);
        bool now_selected = false;
        if(ui_widget.Selection_tool_combo_box->currentIndex() != 2)
          now_selected = visualizer->is_selected (vsp);
        else if(edit_box)
        {
         if(p.x() >= edit_box->point(0,0) - offset.x &&
            p.y() >= edit_box->point(0,1) - offset.y &&
            p.z() <= edit_box->point(0,2) - offset.z &&
            p.x() <= edit_box->point(6,0) - offset.x &&
            p.y() <= edit_box->point(6,1) - offset.y &&
            p.z() >= edit_box->point(6,2) - offset.z
            )
         {
           now_selected = true;
         }
        }

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
    Point_set::Byte_map red, green, blue;
    Point_set::Double_map fred, fgreen, fblue;
    if (point_set_item->point_set()->has_colors())
      {
        if (point_set_item->point_set()->has_byte_colors())
          {
            red = new_item->point_set()->add_property_map<unsigned char>("red", 0).first;
            green = new_item->point_set()->add_property_map<unsigned char>("green", 0).first;
            blue = new_item->point_set()->add_property_map<unsigned char>("blue", 0).first;
          }
        else
          {
            fred = new_item->point_set()->add_property_map<double>("red", 0).first;
            fgreen = new_item->point_set()->add_property_map<double>("green", 0).first;
            fblue = new_item->point_set()->add_property_map<double>("blue", 0).first;
          }
        new_item->point_set()->check_colors(); 
      }
    
    new_item->setColor(point_set_item->color());
    new_item->setRenderingMode(point_set_item->renderingMode());
    new_item->setVisible(point_set_item->visible());

    typedef Point_set_3<Kernel> Point_set;
    for(Point_set::iterator it = point_set_item->point_set()->first_selected ();
	it != point_set_item->point_set()->end(); ++ it)
      {
        Point_set::iterator new_point =
          new_item->point_set()->insert(point_set_item->point_set()->point(*it));
        if (point_set_item->has_normals())
          new_item->point_set()->normal(*new_point) = point_set_item->point_set()->normal(*it);
        if (point_set_item->point_set()->has_colors())
          {
            if (point_set_item->point_set()->has_byte_colors())
              {
                red[*new_point] = (unsigned char)(255. * point_set_item->point_set()->red(*it));
                green[*new_point] = (unsigned char)(255. * point_set_item->point_set()->green(*it));
                blue[*new_point] = (unsigned char)(255. * point_set_item->point_set()->blue(*it));
              }
            else
              {
                fred[*new_point] = point_set_item->point_set()->red(*it);
                fgreen[*new_point] = point_set_item->point_set()->green(*it);
                fblue[*new_point] = point_set_item->point_set()->blue(*it);
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
    if(index == 2)
    {
      edit_box = new Scene_edit_box_item(scene);
      edit_box->setRenderingMode(Wireframe);
      edit_box->setName("Selection Box");
      connect(edit_box, &Scene_edit_box_item::aboutToBeDestroyed,
              this, &Polyhedron_demo_point_set_selection_plugin::reset_editbox);
      scene->addItem(edit_box);
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      viewer->installEventFilter(edit_box);
      ui_widget.Add_to_selection_Button->setVisible(true);
    }
    else
    {
      if(edit_box)
        scene->erase(scene->item_id(edit_box));
      edit_box = NULL;
      ui_widget.Add_to_selection_Button->setVisible(false);
    }


  }

  void on_Selection_mode_combo_box_changed (int index)
  {
    selection_mode = index;
  }

  void reset_editbox()
  {
    edit_box = NULL;
    ui_widget.Add_to_selection_Button->setVisible(false);
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
  Scene_edit_box_item* edit_box;

}; // end Polyhedron_demo_point_set_selection_plugin

//Q_EXPORT_PLUGIN2(Polyhedron_demo_point_set_selection_plugin, Polyhedron_demo_point_set_selection_plugin)

#include "Point_set_selection_plugin.moc"
