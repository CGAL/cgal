#include <CGAL/basic.h>

// STL headers
#include <list>

// CGAL kernel headers
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

// Visibility_complex_2 headers
#include <CEP/Visibility_complex/Visibility_complex_point_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>
#include <CEP/Visibility_complex/Visibility_complex_antichain.h>

// CGAL::Qt_widget headers
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_line.h>
#include "Qt_layer_show_points.h"
#include "Qt_layer_show_lines.h"
#include "Show_pseudo_triangulation.h"

// Qt headers
#include <qapplication.h>
#include <qhbox.h>
#include <qtimer.h>
#include <qtooltip.h>

// ********* TYPEDEFS ********** //
// kernel definition
typedef double FT;
typedef CGAL::Simple_cartesian<FT> K1;
typedef CGAL::Filtered_kernel<K1> Kernel;
struct K : public Kernel {};

// needed kernel objects
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

// vivibility complex definition
typedef CGAL::Visibility_complex_point_traits<K> Gt;
typedef CGAL::Visibility_complex_2<Gt> Visibility_complex;
typedef Visibility_complex::Antichain Antichain;
typedef Visibility_complex::Vertex Vertex;

// ********* MAIN WINDOW CLASS  ********** //
class MyWindow : public QMainWindow
{
  Q_OBJECT
private:
  typedef std::list<Point_2> List_of_points;
  typedef std::list<Line_2> List_of_lines;

  List_of_points* list_of_points;
  List_of_lines* list_of_lines;
  CGAL::Qt_widget* points_widget;
  CGAL::Qt_widget* lines_widget;
  Antichain* antichain;

  QTimer* timer;

public:
  MyWindow(): QMainWindow(0, "main window")
    {
      // title
      setCaption("Arrangements");

      // member datas
      list_of_points = new List_of_points();
      list_of_lines = new List_of_lines();
      antichain = new Antichain();
      timer = new QTimer(this, "timer");

      // widgets
      QHBox* hbox = new QHBox(this, "hbox");
      points_widget = new CGAL::Qt_widget(hbox, "Points");
      lines_widget = new CGAL::Qt_widget(hbox, "Lines");
      hbox->setSpacing(5);

      // create std toolbars
      CGAL::Qt_widget_standard_toolbar* points_stdbar =
	new CGAL::Qt_widget_standard_toolbar(points_widget, this);
      CGAL::Qt_widget_standard_toolbar* lines_stdbar =
	new CGAL::Qt_widget_standard_toolbar(lines_widget, this);
      // move them left and right
      this->moveToolBar(points_stdbar->toolbar(), Left);
      this->moveToolBar(lines_stdbar->toolbar(), Right);
      // add tooltips to distinguish them
      QToolTip::add(points_stdbar->toolbar(), "Primal (left)");
      QToolTip::add(lines_stdbar->toolbar(), "Dual (left)");

      // set points_widget layers
      CGAL::Qt_widget_get_point<K>* get_point = 
	new CGAL::Qt_widget_get_point<K>();
      points_widget->attach(get_point);
      connect(points_widget, SIGNAL(new_cgal_object(CGAL::Object)),
	      this, SLOT(get_cgal_object(CGAL::Object)));
      CGAL::Show_pseudo_triangulation<Antichain>* show_ptrig =
	new CGAL::Show_pseudo_triangulation<Antichain>(antichain, false,
						       CGAL::RED,
						       1);
      points_widget->attach(show_ptrig);
      CGAL::Qt_layer_show_points<List_of_points>* show_points =
	new CGAL::Qt_layer_show_points<List_of_points>(list_of_points,
						       CGAL::BLUE, 4);
      points_widget->attach(show_points);

      // set lines_widget layers
      CGAL::Qt_widget_get_line<K>* get_line = 
	new CGAL::Qt_widget_get_line<K>();
      lines_widget->attach(get_line);
      connect(lines_widget, SIGNAL(new_cgal_object(CGAL::Object)),
	      this, SLOT(get_cgal_object(CGAL::Object)));

      CGAL::Qt_layer_show_lines<List_of_lines>* show_lines = 
	new CGAL::Qt_layer_show_lines<List_of_lines>(list_of_lines,
						       CGAL::BLUE, 1);
      lines_widget->attach(show_lines);
      CGAL::Show_pseudo_triangulation<Antichain>* show_ptrig_dual =
	new CGAL::Show_pseudo_triangulation<Antichain>(antichain, true,
						       CGAL::RED,
						       5);
      lines_widget->attach(show_ptrig_dual);

      // set mouse tracking, need for get_line and show_mouse_coordinates
      lines_widget->setMouseTracking(TRUE);
      points_widget->setMouseTracking(TRUE);

      connect(timer, SIGNAL(timeout()),
	      this, SLOT(sweep_all_minimals()));

      // final settings
      setCentralWidget(hbox);
      resize(700,500);
      hbox->show();
      points_widget->set_window(-1., 1., -1., 1.);
      lines_widget->set_window(-1., 1., -1., 1.);
    }
private slots:
  void get_cgal_object(CGAL::Object obj)
    {
      Point_2 p;
      Line_2 l;
      
      if(CGAL::assign(p, obj))
	{
	  list_of_points->push_back(p);

	  // duality p(a,b) -> l: y=ax-b -> l: -ax+y+b=0
	  l = Line_2(-p.x(), 1, p.y()); // WARNING: kernel dependant!

	  list_of_lines->push_back(l);
	}
      else 
	if(CGAL::assign(l, obj))
	   {
	     if (l.b() == 0) return;
	     // duality l: ax+by+c=0 -> l: y = (-a/b)x - (c/b) -> p(-a/b,c/b)
	     p = Point_2( -l.a()/l.b(), l.c()/l.b() );
	     list_of_lines->push_back(l);
	     list_of_points->push_back(p);
	   }
      std::list<Vertex> empty;
      if(list_of_points->size()>1)
	{
	  timer->stop();
	  delete antichain;
	  antichain = new Antichain(list_of_points->begin(),
				    list_of_points->end(),
				    empty.begin(),
				    empty.end());
	  timer->start(300);
	}
      points_widget->redraw();
      lines_widget->redraw();
    };

  void sweep_all_minimals()
    {
      antichain->sweep_all_minimals();
      points_widget->redraw();
      lines_widget->redraw();
    }

  void clear_scene()
    {
      list_of_points->clear();
      list_of_lines->clear();
      delete antichain;
      antichain= new Antichain();
 
      timer->stop();

      points_widget->redraw();
      lines_widget->redraw();
     }
};

#include "arrangments.moc"
//moc_source_file: arrangments.C

int main(int argc, char** argv)
{
  QApplication app( argc, argv );
  MyWindow* W = new MyWindow();
  app.setMainWidget(W);
  W->show();

  return app.exec();
}
