#include <CGAL/basic.h>

// STL headers
#include <list>

// CGAL kernel headers
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

// CGAL utilities headers
#include <CGAL/iterator.h>

// Visibility_complex_2 headers
#include <CEP/Visibility_complex/Visibility_complex_point_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>
#include <CEP/Visibility_complex/Visibility_complex_antichain.h>

// Delaunay triangulations headers
#include <CGAL/Delaunay_triangulation_2.h>

// CGAL::Qt_widget headers
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_line.h>
#include "Qt_layer_show_points.h"
#include "Qt_layer_show_lines.h"
#include "Show_pseudo_triangulation.h"
#include "Qt_layer_show_nearest_vertex.h"

// Qt headers
#include <qapplication.h>
#include <qmainwindow.h>
#include <qsplitter.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qpushbutton.h>
#include <qframe.h>
#include <qhbox.h>
#include <qtimer.h>
#include <qtooltip.h>
#include <qpixmap.h>

// pixmaps
#include "pixmaps/flip.xpm"
#include "pixmaps/flip_parallele.xpm"
#include "pixmaps/cw_flip.xpm"
#include "pixmaps/cw_flip_parallele.xpm"
#include "pixmaps/gmin.xpm"

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

// visibility complex definition
typedef CGAL::Visibility_complex_point_traits<K> Gt;
typedef CGAL::Visibility_complex_2<Gt> Visibility_complex;
typedef Visibility_complex::Antichain Antichain;
typedef Antichain::Cw_traits Cw_traits;
typedef Visibility_complex::Vertex Vertex;

// Delaunay triangulations definition
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef DT::Vertex DT_vertex;

// utility global functions
const Line_2 dual(const Point_2& p)
{
  // duality p(a,b) -> l: y=ax-b -> l: -ax+y+b=0
  return Line_2(-p.x(), 1, p.y()); // WARNING: kernel dependant!
}

const Point_2 dual(const Line_2& l)
{
  // duality l: ax+by+c=0 -> l: y = (-a/b)x - (c/b) -> p(-a/b,c/b)
  return Point_2( -l.a()/l.b(), l.c()/l.b() );
}

// ********* MAIN WINDOW CLASS  ********** //
class MyWindow : public QMainWindow
{
  Q_OBJECT
private:
  typedef std::list<Point_2> List_of_points;
  typedef std::list<Line_2> List_of_lines;

  // typdefs of layers
  typedef CGAL::Qt_layer_show_points<List_of_points,
    List_of_points::const_iterator> Show_points_from_list;
  typedef CGAL::Qt_layer_show_nearest_vertex<DT, Line_2> Nearest_vertex;

  List_of_points* list_of_points;
  List_of_lines* list_of_lines;
  DT dt_of_points;
  DT dt_of_bitangents;
  CGAL::Qt_widget* points_widget;
  CGAL::Qt_widget* lines_widget;
  Antichain* antichain;

  QTimer* timer;

  // show points from a 2d-triangulation
  struct Vertex_to_Point_2 {
    typedef Point_2 result_type;
    typedef DT_vertex argument_type;

    inline Point_2 operator()(const DT_vertex& v)
      { return v.point(); }
  };
  
public:
  MyWindow()
    : QMainWindow(0, "main window"),
      dt_of_points(),
      dt_of_bitangents()
    {
      // title
      setCaption("Arrangements");

      // member datas
      list_of_points = new List_of_points();
      list_of_lines = new List_of_lines();
      antichain = new Antichain();
      timer = new QTimer(this, "timer");

      // widgets
      QSplitter* horizontal_splitter = 
	new QSplitter(this, "hor_splitter");
      //      QHBox* hbox = new QHBox(this, "hbox");
      points_widget = new CGAL::Qt_widget(horizontal_splitter, "Points");
      lines_widget = new CGAL::Qt_widget(horizontal_splitter, "Lines");
      //      hbox->setSpacing(5);
      //      hbox->setMargin(5);

      // create std toolbars
      CGAL::Qt_widget_standard_toolbar* points_stdbar =
	new CGAL::Qt_widget_standard_toolbar(points_widget, this, Left);
      CGAL::Qt_widget_standard_toolbar* lines_stdbar =
	new CGAL::Qt_widget_standard_toolbar(lines_widget, this, Right);
      // add tooltips to distinguish them
      QToolTip::add(points_stdbar->toolbar(), "Primal (left)");
      QToolTip::add(lines_stdbar->toolbar(), "Dual (left)");

      // create the main toolbar
      QToolBar* maintoolbar = new QToolBar("Main toolbar", this);
      QPushButton* pbClear = new QPushButton("Clear",
					     maintoolbar, "pbClear");
      QToolTip::add(pbClear, "Clear the scene");
      connect(pbClear, SIGNAL(clicked()),
	      this, SLOT(clear_scene()));
      new QToolButton(QPixmap((const char**)flip_xpm),
		      "Flip minimal",
		      "Flip one minimal bitangent of "
		      "the pseudo-triangulation",
		      this, SLOT(sweep_one_minimal()),
		      maintoolbar, "tb_flip_minimal");
      new QToolButton(QPixmap((const char**)flip_parallele_xpm),
		      "Flip all minimals",
		      "Flip all minimal bitangents of "
		      "the pseudo-triangulation",
		      this, SLOT(sweep_all_minimals()),
		      maintoolbar, "tb_flip_all_minimals");
      new QToolButton(QPixmap((const char**)cw_flip_xpm),
		      "CW Flip minimal",
		      "Flip clockwise one minimal bitangent of "
		      "the pseudo-triangulation",
		      this, SLOT(sweep_cw_one_minimal()),
		      maintoolbar, "tb_cw_flip_minimal");
      new QToolButton(QPixmap((const char**)cw_flip_parallele_xpm),
		      "CW Flip all minimals",
		      "Flip clockwise all minimal bitangents of "
		      "the pseudo-triangulation",
		      this, SLOT(sweep_cw_all_minimals()),
		      maintoolbar, "tb_cw_flip_all_minimals");
      new QToolButton(QPixmap((const char**)gmin_xpm),
		      "Compute Gmin",
		      "Compute the initial pseudo-triangulation",
		      this, SLOT(compute_gmin()),
		      maintoolbar, "tb_gmin");
      
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
      Show_points_from_list* show_points = 
	new Show_points_from_list(list_of_points,
				  &List_of_points::begin,
				  &List_of_points::end,
				  CGAL::BLUE, 4);
      points_widget->attach(show_points);
      
      Nearest_vertex* nearest_point = 
	new Nearest_vertex(dt_of_points,
			   0,
			   Qt::green,
			   4,
			   CGAL::DISC);
      points_widget->attach(nearest_point);

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

      Nearest_vertex* nearest_bitangent = 
	new Nearest_vertex(dt_of_bitangents,
			   nearest_point,
			   Qt::green,
			   5,
			   CGAL::DISC);
      lines_widget->attach(nearest_bitangent);
      nearest_point->set_twin(nearest_bitangent);

      // set mouse tracking, need for get_line and show_mouse_coordinates
      lines_widget->setMouseTracking(TRUE);
      points_widget->setMouseTracking(TRUE);

      connect(timer, SIGNAL(timeout()),
	      this, SLOT(sweep_all_minimals()));

      // final settings
      setCentralWidget(horizontal_splitter);
      setUsesBigPixmaps(true);
      resize(700,500);
      //      hbox->show();

      points_widget->set_window(-1., 1., -1., 1.);
      lines_widget->set_window(-1., 1., -1., 1.);
    }
private slots:
  void get_cgal_object(CGAL::Object obj)
    {
      Point_2 p;
      Line_2 l;
      
      if(CGAL::assign(p, obj))
	  l = dual(p);
      else 
	if(CGAL::assign(l, obj))
	   {
	     if (l.b() == 0) return;
	     p = dual(l);
	   }
      list_of_lines->push_back(l);
      list_of_points->push_back(p);
      dt_of_points.insert(p);

      if(list_of_points->size()>1)
	{
	  compute_gmin();
	  // the slot call redraw() itself
	}
      else
	{
	  points_widget->redraw();
	  lines_widget->redraw();
	}
    };

  void sweep_one_minimal()
    {
      antichain->sweep(&(*(antichain->minimals_begin()))); 
      //TODO: berk!
      update_dt_of_bitangents();
      points_widget->redraw();
      lines_widget->redraw();
    }

  void sweep_all_minimals()
    {
      antichain->sweep_all_minimals();
      update_dt_of_bitangents();
      points_widget->redraw();
      lines_widget->redraw();
    }

  void sweep_cw_one_minimal()
    {
      antichain->sweep(&(*(antichain->cw_minimals_begin())),
		       Cw_traits()); //TODO: berk!
      update_dt_of_bitangents();
      points_widget->redraw();
      lines_widget->redraw();
    }

  void sweep_cw_all_minimals()
    {
      antichain->sweep_all_minimals(Cw_traits());
      points_widget->redraw();
      lines_widget->redraw();
    }

  void compute_gmin()
    {
      bool active_timer = timer->isActive();
      timer->stop();
      delete antichain;
      
      std::list<Vertex> empty;
      antichain = new Antichain(list_of_points->begin(),
				list_of_points->end(),
				empty.begin(),
				empty.end());
      update_dt_of_bitangents();
      if(active_timer)
	timer->start(300);
      points_widget->redraw();
      lines_widget->redraw();
    }

  void update_dt_of_bitangents()
    {
      typedef Antichain::Vertex_iterator Vertex_iterator;
      
      dt_of_bitangents.clear();
      if(antichain==0) return;
      
      for(Vertex_iterator it = antichain->vertices_begin();
	  it!=antichain->vertices_end();
	  ++it)
	{
	  Line_2 l = it->supporting_line();
	  if(l.b() != 0)
	    dt_of_bitangents.insert(Point_2( -l.a()/l.b(), l.c()/l.b() ));
	}
    }

  void clear_scene()
    {
      list_of_points->clear();
      dt_of_points.clear();
      list_of_lines->clear();
      dt_of_bitangents.clear();
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
