#ifdef CGAL_USE_QT

#include <iostream>
#include <fstream>

#include <qapplication.h>
#include <qmainwindow.h>
#include <qpixmap.h>

//#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_circle.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_Apollonius_site_2.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>


#include "typedefs.h"

#include "qt_file_toolbar.h"
#include "qt_layers_toolbar.h"
#include "qt_layers.h"


//************************************
// global variables
//************************************
AG_2 ag;

//************************************
// conversion functions
//************************************
inline Apollonius_site_2
to_site(const Circle_2 &c)
{
  double r = CGAL_NTS sqrt(CGAL_NTS to_double(c.squared_radius()));
  return  Apollonius_site_2(c.center(), Rep::RT(r));
}

inline Circle_2
to_circle(const Apollonius_site_2 &wp)
{
  return Circle_2(wp.point(), CGAL_NTS square(wp.weight()) );
}


//************************************
// my window
//************************************
class My_Window : public QMainWindow {
  Q_OBJECT

  friend class Layers_toolbar;
private:
  CGAL::Qt_widget *widget;
  Layers_toolbar *layers_toolbar;
  File_toolbar   *file_toolbar;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  CGAL::Qt_widget_get_circle<Rep> get_circle;
  CGAL::Qt_widget_get_point<Rep> get_point;
  bool is_edit_mode;
  bool is_remove_mode;
  bool is_insert_point_mode;

public:
  My_Window(int x, int y)
  {
    is_edit_mode = false;
    is_remove_mode = false;
    is_insert_point_mode = false;

    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);

    *widget << CGAL::BackgroundColor(CGAL::YELLOW);
    resize(x,y);
    widget->set_window(0, x, 0, y);
    widget->show();

    //    setUsesBigPixmaps(TRUE);

    //How to attach the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this,
						    this, FALSE, "");

    file_toolbar = new File_toolbar("File Operations",
				    this, this, FALSE,
				    "File Operations");

    layers_toolbar = new Layers_toolbar(widget, ag,
					"Geometric Operations",
					this, this, FALSE,
					"Geometric Operations");

    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), this,
	    SLOT(get_object(CGAL::Object)));

    connect(layers_toolbar, SIGNAL(inputModeChanged(bool)), this,
	    SLOT(get_input_mode(bool)));

    connect(layers_toolbar, SIGNAL(removeModeChanged(bool)), this,
	    SLOT(get_remove_mode(bool)));

    connect(layers_toolbar, SIGNAL(editModeChanged(bool)), this,
    	    SLOT(get_edit_mode(bool)));

    connect(file_toolbar, SIGNAL(fileToRead(const QString&)), this,
	    SLOT(read_from_file(const QString&)));

    connect(file_toolbar, SIGNAL(fileToWrite(const QString&)), this,
	    SLOT(write_to_file(const QString&)));

    connect(file_toolbar, SIGNAL(printScreen()), this,
	    SLOT(print_screen()));
    connect(file_toolbar, SIGNAL(clearAll()), this,
	    SLOT(remove_all()));

    widget->attach(&get_circle);

    setMouseTracking(true);
    widget->setMouseTracking(true);

    widget->attach(&get_point);

    get_circle.activate();
    get_point.deactivate();

    //    get_circle.setMouseTracking(true);
  }
  ~My_Window(){}

  void set_window(double xmin, double xmax,
		  double ymin, double ymax)
  {
    widget->set_window(xmin, xmax, ymin, ymax);
  }

private slots:
  void get_object(CGAL::Object obj)
  {
    if ( is_edit_mode ) { return; }
    if ( is_remove_mode ) {
      if ( ag.number_of_vertices() == 0 ) { return; }
      Point_2 p;
      if ( CGAL::assign(p, obj) ) {
	Vertex_handle v = ag.nearest_neighbor(p);
#if 0
	AG_2::Vertex_circulator vc = v->incident_vertices();
	AG_2::Vertex_circulator vc_start = vc;
	widget->redraw();
	*widget << CGAL::PURPLE;
	std::cout << "==============================" << std::endl;
	std::cout << v->point() << " " << v->timestamp() << std::endl;
	std::cout << "------------------------------" << std::endl;
	do {
	  Vertex_handle v1(vc);
	  if ( !ag.is_infinite(v1) ) {
	    *widget << to_circle(v1->point());
	    *widget << v1->point().point();
	    std::cout << v1->point() << " " << v1->timestamp() << std::endl;
	  }
	  ++vc;
	} while ( vc != vc_start );
	std::cout << "==============================" << std::endl;
#else
	if ( &(*v) != NULL ) {
	  ag.remove(v);
	  assert( ag.is_valid(false,1) );
	  //	  std::cout << "--------" << std::endl;
	}
#endif
      }
      widget->redraw();
      return;
    }
    Circle_2 c;
    Point_2 p;
    if ( CGAL::assign(c, obj) ) {
      Apollonius_site_2 wp = to_site(c);
      ag.insert(wp);
    } else if ( CGAL::assign(p, obj) ) {
      Apollonius_site_2 wp(p, Weight(0));
      ag.insert(wp);
    }
    assert( ag.is_valid(false, 1) );
    //    std::cout << "--------" << std::endl;
    //      *widget << CGAL::RED << c;
    widget->redraw();
  }

  void get_input_mode(bool b)
  {
    is_insert_point_mode = b;

    if ( !is_remove_mode && !is_edit_mode ) {
      if ( is_insert_point_mode ) {
	get_point.activate();
	get_circle.deactivate();
      } else {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void get_remove_mode(bool b)
  {
    is_remove_mode = b;

    if ( is_remove_mode ) {	
      get_point.activate();
      get_circle.deactivate();
    } else {
      if ( !is_insert_point_mode ) {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void get_edit_mode(bool b)
  {
    is_edit_mode = b;

    if ( is_edit_mode ) {	
      get_point.activate();
      get_circle.deactivate();
    } else {
      if ( !is_insert_point_mode ) {
	get_point.deactivate();
	get_circle.activate();
      }
    }
  }

  void read_from_file(const QString& fileName)
  {
    std::ifstream f(fileName);
    assert( f );

    //      int n;
    //      f >> n;
    Apollonius_site_2 wp;

    int counter = 0;
    std::cout << std::endl;

    while ( f >> wp ) {
      ag.insert(wp);
      counter++;
      if ( counter % 500 == 0 ) {
	std::cout << "\r" << counter
		  << " sites have been inserted..." << std::flush;
      }
    }
    std::cout << "\r" << counter 
	      << " sites have been inserted... Done!"
	      << std::endl;
    assert( ag.is_valid(false, 1) );
    widget->redraw();
  }

  void write_to_file(const QString& fileName)
  {
    std::ofstream f(fileName);
    assert( f );
    f.precision(18);

    for( AG_2::Sites_iterator it = ag.sites_begin();
	 it != ag.sites_end(); it++ ) {
      f << (*it) << std::endl;
    }
  }

  void print_screen()
  {
    widget->print_to_ps();
  }

  void remove_all()
  {
    ag.clear();
    widget->redraw();
  }

};

#include "qt_file_toolbar.moc"
#include "qt_layers_toolbar.moc"
#include "apollonius_graph_2.moc"

int
main(int argc, char* argv[])
{
  int size = 750;

  QApplication app( argc, argv );
  My_Window W(size,size);
  app.setMainWidget( &W );
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W.setIcon(cgal_icon);
  W.show();
  W.set_window(0,size,0,size);
  W.setCaption("Apollonius diagram 2");
  W.setMouseTracking(TRUE);
  return app.exec();
}


// moc_source_file: apollonius_graph_2.C

#else

#include <iostream>

int
main(int argc, char* argv[])
{
  std::cerr << "This demo needs CGAL's Qt_widget installed "
	    << "in order to run..."
	    << std::endl << std::endl;
  return 0;
}


#endif
